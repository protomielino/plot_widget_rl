#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

typedef struct
{
    Rectangle viewport_px;   // area in pixel sullo schermo dedicata al grafico
    double xmin, xmax;       // vista in coordinate mondo (x)
    double ymin, ymax;       // vista in coordinate mondo (y)
    bool dragging;
    Vector2 last_mouse;
} plotWidget;

typedef struct
{
    double x, y;
} pointD;

// utility functions: world <-> screen transformations
static inline double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}
static inline double clamp_double(double v, double a, double b)
{
    return (v < a) ? a : (v > b) ? b : v;
}

static Vector2 world_to_screen(const plotWidget *w, double wx, double wy)
{
    float sx = (float)( (wx - w->xmin) / (w->xmax - w->xmin) * w->viewport_px.width + w->viewport_px.x );
    float sy = (float)( w->viewport_px.y + w->viewport_px.height - (wy - w->ymin) / (w->ymax - w->ymin) * w->viewport_px.height );
    return (Vector2){ sx, sy };
}

static void screen_to_world(const plotWidget *w, double sx, double sy, double *wx, double *wy)
{
    double nx = (sx - w->viewport_px.x) / w->viewport_px.width;
    double ny = 1.0 - (sy - w->viewport_px.y) / w->viewport_px.height;
    *wx = lerp(w->xmin, w->xmax, nx);
    *wy = lerp(w->ymin, w->ymax, ny);
}

// tick spacing utility: ottiene un "nice" step per ticks in base al range e numero minimo di ticks ---
static double nice_tick_step(double range, int target_ticks)
{
    if (range <= 0) return 0.0;
    double raw = range / (double)target_ticks;
    double mag = pow(10.0, floor(log10(raw)));
    double norm = raw / mag;
    double nice;
    if (norm < 1.5) nice = 1.0;
    else if (norm < 3.0) nice = 2.0;
    else if (norm < 7.0) nice = 5.0;
    else nice = 10.0;
    return nice * mag;
}

// threshold value in pixel: subdivide if distance of real point from line is > thresh_px
#define SUBDIV_THRESH_PX 0.6
// max recursion depth
#define SUBDIV_MAX_DEPTH 12

// example function y = f(x). substitute funzione with callback(?).
static double eval_func(double x)
{
    return sin(x) + 0.2 * sin(3.0 * x);
}

// calculate distance point->segment in pixel
static double point_line_distance_px(const plotWidget *w, double x0, double y0, double x1, double y1, double xm, double ym)
{
    Vector2 a = world_to_screen(w, x0, y0);
    Vector2 b = world_to_screen(w, x1, y1);
    Vector2 m = world_to_screen(w, xm, ym);
    // distance of point m from line ab
    float dx = b.x - a.x, dy = b.y - a.y;
    float len2 = dx*dx + dy*dy;
    if (len2 == 0) return Vector2Distance(m, a);
    float t = ((m.x - a.x)*dx + (m.y - a.y)*dy) / len2;
    if (t < 0) return Vector2Distance(m, a);
    if (t > 1) return Vector2Distance(m, b);
    Vector2 proj = { a.x + t*dx, a.y + t*dy };
    return Vector2Distance(m, proj);
}

// recursion for subdivision and draw
static void subdiv_draw(const plotWidget *w, double x0, double y0, double x1, double y1, int depth)
{
    if (depth >= SUBDIV_MAX_DEPTH) {
        Vector2 A = world_to_screen(w, x0, y0);
        Vector2 B = world_to_screen(w, x1, y1);
        DrawLineEx(A, B, 2.0f, (Color){50,120,200,255});
        return;
    }
    double xm = 0.5*(x0 + x1);
    double ym = eval_func(xm);
    double dist_px = point_line_distance_px(w, x0, y0, x1, y1, xm, ym);
    if (dist_px > SUBDIV_THRESH_PX) {
        subdiv_draw(w, x0, y0, xm, ym, depth+1);
        subdiv_draw(w, xm, ym, x1, y1, depth+1);
    } else {
        Vector2 A = world_to_screen(w, x0, y0);
        Vector2 B = world_to_screen(w, x1, y1);
        DrawLineEx(A, B, 2.0f, (Color){50,120,200,255});
    }
}

// adaptive function drawing: partition visible interval in base steps depending on pixel
static void draw_adaptive_curve(const plotWidget *w, double xfunc_min, double xfunc_max)
{
    // clip to viewport
    BeginScissorMode((int)w->viewport_px.x, (int)w->viewport_px.y, (int)w->viewport_px.width, (int)w->viewport_px.height);

    // base step : fraction of pixel in world coordinates
    double dx_pixel = (w->xmax - w->xmin) / w->viewport_px.width;
    double base_step = dx_pixel * 2.0; // ~2 samples per pixel; tweak parameter
    if (base_step <= 0) base_step = (xfunc_max - xfunc_min) / 200.0;

    // check boundaries
    double x0 = xfunc_min;
    while (x0 < xfunc_max) {
        double x1 = x0 + base_step;
        if (x1 > xfunc_max) x1 = xfunc_max;
        double y0 = eval_func(x0);
        double y1 = eval_func(x1);

        // filter: skip if both y are NaN or inf
        if (!isfinite(y0) || !isfinite(y1)) {
            x0 = x1;
            continue;
        }

        subdiv_draw(w, x0, y0, x1, y1, 0);
        x0 = x1;
    }

    EndScissorMode();
}

// draw dynamic grid and ticks
static void draw_grid_and_ticks(const plotWidget *w)
{
    // label's inside margin

    const int target_ticks_x = 8;
    const int target_ticks_y = 6;

    double xrange = w->xmax - w->xmin;
    double yrange = w->ymax - w->ymin;

    double stepx = nice_tick_step(xrange, target_ticks_x);
    double stepy = nice_tick_step(yrange, target_ticks_y);

    // find first tick >= xmin
    double startx = ceil(w->xmin / stepx) * stepx;
    double starty = ceil(w->ymin / stepy) * stepy;

    // colors
    Color grid_col = (Color){200,200,200,80};
    Color axis_col = (Color){180,50,50,200};
    Color tick_col = (Color){80,80,80,200};
    Color label_col = (Color){40,40,40,255};

    // clip to viewport
    BeginScissorMode((int)w->viewport_px.x, (int)w->viewport_px.y, (int)w->viewport_px.width, (int)w->viewport_px.height);

    // vertical grid lines + ticks (x axis)
    for (double tx = startx; tx <= w->xmax + stepx*0.5; tx += stepx) {
        Vector2 p1 = world_to_screen(w, tx, w->ymin);
        Vector2 p2 = world_to_screen(w, tx, w->ymax);
        DrawLineEx(p1, p2, 1.0f, grid_col);

        // tick on bottom
        Vector2 tl = (Vector2){ p1.x, p1.y };
        DrawLine((int)tl.x, (int)(tl.y), (int)tl.x, (int)(tl.y - 6), tick_col);

        // label
        char buf[64];
        // "smart" formatting: use %.2g for big/small values
        if (fabs(tx) < 1e-3 || fabs(tx) > 1e4) {
            snprintf(buf, sizeof(buf), "%.3g", tx);
        } else {
            // scegli precisione in base a stepx
            int prec = (int)ceil(-floor(log10(stepx)));
            if (prec < 0) prec = 0;
            char fmt[16]; snprintf(fmt, sizeof(fmt), "%%.%df", prec);
            snprintf(buf, sizeof(buf), fmt, tx);
        }
        Vector2 lb = (Vector2){ p1.x + 4, p1.y + 4 };
        DrawText(buf, (int)lb.x - 6, (int)lb.y - 20, 10, label_col);
    }

    // horizontal grid lines + ticks (y axis)
    for (double ty = starty; ty <= w->ymax + stepy*0.5; ty += stepy) {
        Vector2 p1 = world_to_screen(w, w->xmin, ty);
        Vector2 p2 = world_to_screen(w, w->xmax, ty);
        DrawLineEx(p1, p2, 1.0f, grid_col);

        // tick on left
        Vector2 tl = (Vector2){ p1.x, p1.y };
        DrawLine((int)(tl.x), (int)tl.y, (int)(tl.x + 6), (int)tl.y, tick_col);

        char buf[64];
        if (fabs(ty) < 1e-3 || fabs(ty) > 1e4) {
            snprintf(buf, sizeof(buf), "%.3g", ty);
        } else {
            int prec = (int)ceil(-floor(log10(stepy)));
            if (prec < 0) prec = 0;
            char fmt[16]; snprintf(fmt, sizeof(fmt), "%%.%df", prec);
            snprintf(buf, sizeof(buf), fmt, ty);
        }
        Vector2 lb = world_to_screen(w, w->xmin, ty);
        DrawText(buf, (int)(lb.x + 10), (int)(lb.y - 6), 10, label_col);
    }

    // draw axis x=0 e y=0, if visible
    if (w->xmin <= 0 && w->xmax >= 0) {
        Vector2 a = world_to_screen(w, 0.0, w->ymin);
        Vector2 b = world_to_screen(w, 0.0, w->ymax);
        DrawLineEx(a, b, 2.0f, axis_col);
    }
    if (w->ymin <= 0 && w->ymax >= 0) {
        Vector2 a = world_to_screen(w, w->xmin, 0.0);
        Vector2 b = world_to_screen(w, w->xmax, 0.0);
        DrawLineEx(a, b, 2.0f, axis_col);
    }

    EndScissorMode();
}

// draw curve data (w\lines) with viewport clipping
static void draw_curve(const plotWidget *w, const pointD *pts, size_t n)
{
    if (n < 2) return;
    BeginScissorMode((int)w->viewport_px.x, (int)w->viewport_px.y, (int)w->viewport_px.width, (int)w->viewport_px.height);
    Color curve_col = (Color){50,120,200,255};
    for (size_t i = 1; i < n; ++i) {
        Vector2 a = world_to_screen(w, pts[i-1].x, pts[i-1].y);
        Vector2 b = world_to_screen(w, pts[i].x, pts[i].y);
        DrawLineEx(a, b, 2.0f, curve_col);
    }
    EndScissorMode();
}

// handle inputs: pan & zoom
static void handle_input(plotWidget *w)
{
    Vector2 mpos = GetMousePosition();

    // start dragging with middle mouse button click
    if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
        if (CheckCollisionPointRec(mpos, w->viewport_px)) {
            w->dragging = true;
            w->last_mouse = mpos;
            SetMouseCursor(MOUSE_CURSOR_DEFAULT);
        }
    }
    if (IsMouseButtonReleased(MOUSE_BUTTON_MIDDLE)) {
        w->dragging = false;
        SetMouseCursor(MOUSE_CURSOR_CROSSHAIR);
    }
    if (w->dragging) {
        Vector2 delta = { mpos.x - w->last_mouse.x, mpos.y - w->last_mouse.y };
        // map delta in world space (inverse traslation)
        double dx_world1, dy_world1, dx_world2, dy_world2;
        screen_to_world(w, w->last_mouse.x, w->last_mouse.y, &dx_world1, &dy_world1);
        screen_to_world(w, mpos.x, mpos.y, &dx_world2, &dy_world2);
        double ddx = dx_world1 - dx_world2;
        double ddy = dy_world1 - dy_world2;
        w->xmin += ddx; w->xmax += ddx;
        w->ymin += ddy; w->ymax += ddy;
        w->last_mouse = mpos;
    }

    // zoom with mouse wheel: zoom centered on mouse position
    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f && CheckCollisionPointRec(mpos, w->viewport_px)) {
        double world_cx, world_cy;
        screen_to_world(w, mpos.x, mpos.y, &world_cx, &world_cy);
        double zoomFactor = (wheel > 0) ? 0.9 : 1.111111; // 10% per notch
        double nxmin = world_cx + (w->xmin - world_cx) * zoomFactor;
        double nxmax = world_cx + (w->xmax - world_cx) * zoomFactor;
        double nymin = world_cy + (w->ymin - world_cy) * zoomFactor;
        double nymax = world_cy + (w->ymax - world_cy) * zoomFactor;

        // optional: limit max/min zoom
        double min_range_x = 1e-12;
        double min_range_y = 1e-12;
        if ((nxmax - nxmin) > min_range_x && (nymax - nymin) > min_range_y) {
            w->xmin = nxmin; w->xmax = nxmax;
            w->ymin = nymin; w->ymax = nymax;
        }
    }
}

// example: generate points for noisy y = sin(x) + ...
static pointD* generate_example_curve(size_t *out_n)
{
    size_t n = 2000;
    pointD *arr = malloc(sizeof(pointD) * n);
    double x0 = -20.0, x1 = 20.0;
    for (size_t i = 0; i < n; ++i) {
        double t = (double)i / (double)(n-1);
        double x = lerp(x0, x1, t);
        double y = sin(x) + 0.2 * sin(3.0 * x) + 0.2 * ((rand() / (double)RAND_MAX) - 0.5);
        arr[i].x = x;
        arr[i].y = y;
    }
    *out_n = n;
    return arr;
}

// draw widget's frame UI (frame, title)
static void draw_widget_frame(const plotWidget *w, const char *title)
{
    // frame
    DrawRectangleLinesEx(w->viewport_px, 2, DARKGRAY);
    // title
    if (title) {
        DrawText(title, (int)w->viewport_px.x + 6, (int)w->viewport_px.y + 6, 10, BLACK);
    }
    // world area info (lower right corner)
    char buf[128];
    snprintf(buf, sizeof(buf), "x:[%.3g, %.3g] y:[%.3g, %.3g]", w->xmin, w->xmax, w->ymin, w->ymax);
    int tw = MeasureText(buf, 10);
    DrawText(buf, (int)(w->viewport_px.x + w->viewport_px.width - tw - 6), (int)(w->viewport_px.y + w->viewport_px.height - 26), 10, DARKGRAY);
}

// main: raylib window and loop
int main(void)
{
    srand((unsigned int)clock());

    const int screenW = 1000;
    const int screenH = 700;
    InitWindow(screenW, screenH, "Plot Widget - raylib");
    SetTargetFPS(30);

    // widget initialization
    plotWidget widget;
    widget.viewport_px = (Rectangle){ 60, 40, screenW - 120.0f, screenH - 100.0f };
    widget.xmin = -10.0; widget.xmax = 10.0;
    widget.ymin = -3.0; widget.ymax = 3.0;
    widget.dragging = false;

    size_t npts;
    pointD *curve = generate_example_curve(&npts);

    // loop
    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(RAYWHITE);

        // UI background per viewport
        DrawRectangleRec(widget.viewport_px, (Color){245,245,245,255});

        handle_input(&widget);

        // Draw grid, ticks, curve, frame
        draw_grid_and_ticks(&widget);
        draw_curve(&widget, curve, npts);
        draw_adaptive_curve(&widget, widget.xmin, widget.xmax);
        draw_widget_frame(&widget, "Example: sin(x) + 0.2 * sin(3.0 * x)");

        // istruzioni interazione
        DrawText("Pan: Middle mouse drag    Zoom: Mouse wheel", 10, 10, 10, DARKGRAY);

        EndDrawing();
    }

    free(curve);

    CloseWindow();
    return 0;
}
