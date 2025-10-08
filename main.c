#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdbool.h>

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

typedef struct
{
    pointD *pts;
    size_t n;
    Color color;
    char name[64];
    bool visible;
} series_t;

// stato selezione/tooltip (inserire in cima file, es. dopo typedefs)
typedef struct {
    bool active;          // true se un punto è selezionato
    size_t series_idx;    // indice della serie selezionata
    size_t pt_idx;        // indice del punto selezionato nella serie
    pointD pt;            // coordinate world del punto selezionato
} selection_t;

static selection_t g_selection = {0};

// area selection state
typedef struct {
    bool active;        // true mentre l'utente sta trascinando l'area
    bool finished;      // true se una area è stata selezionata (oltre al drag in corso)
    Vector2 start_px;   // corner iniziale in screen coords
    Vector2 end_px;     // corner corrente/finale in screen coords
} area_sel_t;

static area_sel_t g_area = {0};

typedef struct {
    bool has_any;
    // per semplicità memorizziamo un array dinamico di (series_idx, pt_idx)
    size_t *series_idx;
    size_t *pt_idx;
    size_t count;
} multi_selection_t;

static multi_selection_t g_multi = {0};

size_t nseries = 0;
series_t *series = NULL;

// utility functions: world <-> screen transformations
static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }
static inline double clamp_double(double v, double a, double b) { return v < a ? a : (v > b ? b : v); }

static Vector2 world_to_screen(const plotWidget *w, double wx, double wy)
{
    double sx = (wx - w->xmin) / (w->xmax - w->xmin) * w->viewport_px.width + w->viewport_px.x;
    double sy = (1.0 - (wy - w->ymin) / (w->ymax - w->ymin)) * w->viewport_px.height + w->viewport_px.y;
    return (Vector2){ (float)sx, (float)sy };
}

static void screen_to_world(const plotWidget *w, double sx, double sy, double *wx, double *wy)
{
    double nx = (sx - w->viewport_px.x) / w->viewport_px.width;
    double ny = (sy - w->viewport_px.y) / w->viewport_px.height;
    *wx = lerp(w->xmin, w->xmax, nx);
    *wy = lerp(w->ymax, w->ymin, ny); // note invert
}

// ritorna rect normalizzato (x<=x+w, y<=y+h)
static Rectangle norm_rect_from_points(Vector2 a, Vector2 b)
{
    float x = fminf(a.x, b.x);
    float y = fminf(a.y, b.y);
    float w = fabsf(a.x - b.x);
    float h = fabsf(a.y - b.y);
    return (Rectangle){ x, y, w, h };
}

// converte rect screen->world (ritorna rect world with xmin<=xmax, ymin<=ymax)
static void rect_screen_to_world_box(const plotWidget *w, Rectangle r, double *out_xmin, double *out_xmax, double *out_ymin, double *out_ymax)
{
    double x0, y0, x1, y1;
    screen_to_world(w, r.x, r.y, &x0, &y0);                     // top-left screen -> world
    screen_to_world(w, r.x + r.width, r.y + r.height, &x1, &y1); // bottom-right screen -> world
    // normalize to xmin/xmax, ymin/ymax
    *out_xmin = fmin(x0, x1);
    *out_xmax = fmax(x0, x1);
    *out_ymin = fmin(y0, y1);
    *out_ymax = fmax(y0, y1);
}

// ritorna distanza in pixel tra mouse screen pos e punto world
static double mouse_point_distance_px(const plotWidget *w, const pointD *p, Vector2 mouse)
{
    Vector2 sp = world_to_screen(w, p->x, p->y);
    double dx = sp.x - mouse.x;
    double dy = sp.y - mouse.y;
    return sqrt(dx*dx + dy*dy);
}

// cerca punto più vicino (entro threshold_px) e riempie out_* se trovato
static bool find_nearest_point(const plotWidget *w, const series_t *series, size_t nseries, Vector2 mouse, double threshold_px, size_t *out_series, size_t *out_ptidx, pointD *out_pt, double *out_dist_px)
{
    bool found = false;
    double bestd = 1e9;
    size_t best_si = 0, best_pi = 0;
    pointD best_pt = {0,0};

    for (size_t si = 0; si < nseries; ++si) {
        if (!series[si].visible) continue;
        const series_t *s = &series[si];
        for (size_t i = 0; i < s->n; ++i) {
            if (!isfinite(s->pts[i].x) || !isfinite(s->pts[i].y)) continue;
            double d = mouse_point_distance_px(w, &s->pts[i], mouse);
            if (d < bestd) {
                bestd = d;
                best_si = si; best_pi = i; best_pt = s->pts[i];
            }
        }
    }
    if (bestd <= threshold_px) {
        *out_series = best_si;
        *out_ptidx = best_pi;
        *out_pt = best_pt;
        *out_dist_px = bestd;
        return true;
    }
    return false;
}

// disegna tooltip (mouse + nearest values per serie)
static void draw_mouse_tooltip(const plotWidget *w, const series_t *series, size_t nseries, Vector2 mouse)
{
    // se mouse fuori viewport non disegnare
    if (mouse.x < w->viewport_px.x || mouse.x > w->viewport_px.x + w->viewport_px.width ||
        mouse.y < w->viewport_px.y || mouse.y > w->viewport_px.y + w->viewport_px.height)
        return;

    // coord world del mouse
    double mx, my;
    screen_to_world(w, mouse.x, mouse.y, &mx, &my);

    // costruisci testo: mouse coords
    char buf[256];
    int yoff = (int)(w->viewport_px.y + w->viewport_px.height + 6); // base per disegno sotto viewport
    snprintf(buf, sizeof(buf), "x=%.6g  y=%.6g", mx, my);
    DrawText(buf, (int)mouse.x + 12, (int)mouse.y - 6, 12, DARKGRAY);

    // per ogni serie visibile, trova punto con x più vicino (interpolazione opzionale)
    float ty = mouse.y + 16;
    for (size_t si = 0; si < nseries; ++si) {
        if (!series[si].visible) continue;
        // ricerca punto con x più vicino (soglia su distanza x in world)
        const series_t *s = &series[si];
        size_t best_i = SIZE_MAX;
        double best_dx = 1e9;
        for (size_t i = 0; i < s->n; ++i) {
            if (!isfinite(s->pts[i].x) || !isfinite(s->pts[i].y)) continue;
            double dx = fabs(s->pts[i].x - mx);
            if (dx < best_dx) { best_dx = dx; best_i = i; }
        }
        if (best_i != SIZE_MAX) {
            snprintf(buf, sizeof(buf), "%s: x=%.6g y=%.6g", s->name, s->pts[best_i].x, s->pts[best_i].y);
            DrawText(buf, (int)mouse.x + 12, (int)ty, 12, s->color);
            ty += 14;
        }
    }
}

// tick spacing utility: ottiene un "nice" step per ticks in base al range e numero minimo di ticks
static double nice_tick_step(double range, int target_ticks)
{
    if (range <= 0) return 0.0;
    double raw = range / target_ticks;
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

/* Format specifiers:
   'B' = base function: sin(x) + 0.2*sin(3*x)         (no args)
   'S' = sine float case: expects (float N, float step)  (two args)
   'C' = cubic example: x^3 / 50
*/
static double eval_func(char format, double x, ...)
{
    (void)x;
    if (format == 'B') {
        return sin(x) + 0.2 * sin(3.0 * x);
    } else if (format == 'C') {
        return (x*x*x) / 50.0;
    } else if (format == 'S') {
        va_list ap;
        va_start(ap, x);
        double a = va_arg(ap, double);
        double step = va_arg(ap, double);
        va_end(ap);
        return a * sin(step * x);
    }
    return 0.0;
}

// calculate distance point->segment in pixel
static double point_line_distance_px(const plotWidget *w, double x0, double y0, double x1, double y1, double xm, double ym)
{
    Vector2 p0 = world_to_screen(w, x0, y0);
    Vector2 p1 = world_to_screen(w, x1, y1);
    Vector2 pm = world_to_screen(w, xm, ym);
    double vx = p1.x - p0.x;
    double vy = p1.y - p0.y;
    double wx = pm.x - p0.x;
    double wy = pm.y - p0.y;
    double c1 = vx * wx + vy * wy;
    double c2 = vx * vx + vy * vy;
    double t = (c2 <= 1e-12) ? 0.0 : c1 / c2;
    t = clamp_double(t, 0.0, 1.0);
    double projx = p0.x + t * vx;
    double projy = p0.y + t * vy;
    double dx = pm.x - projx;
    double dy = pm.y - projy;
    return sqrt(dx*dx + dy*dy);
}

// recursion for subdivision and draw (for function plotting)
static void subdiv_draw(const plotWidget *w, char fmt, double x0, double x1, double y0, double y1, int depth, double extra1, double extra2)
{
    if (depth > SUBDIV_MAX_DEPTH) {
        Vector2 s0 = world_to_screen(w, x0, y0);
        Vector2 s1 = world_to_screen(w, x1, y1);
        DrawLineV(s0, s1, WHITE);
        return;
    }
    double xm = 0.5 * (x0 + x1);
    double ym = eval_func(fmt, xm, extra1, extra2);
    double dist = point_line_distance_px(w, x0, y0, x1, y1, xm, ym);
    if (dist > SUBDIV_THRESH_PX) {
        subdiv_draw(w, fmt, x0, xm, y0, ym, depth+1, extra1, extra2);
        subdiv_draw(w, fmt, xm, x1, ym, y1, depth+1, extra1, extra2);
    } else {
        Vector2 s0 = world_to_screen(w, x0, y0);
        Vector2 s1 = world_to_screen(w, x1, y1);
        DrawLineV(s0, s1, WHITE);
    }
}

// adaptive function drawing: partition visible interval in base steps depending on pixel
static void draw_adaptive_curve(const plotWidget *w, double xfunc_min, double xfunc_max)
{
    // Partition with number of base segments proportional to pixel width
    int base_segments = (int)fmax(32, w->viewport_px.width / 10.0);
    double dx = (xfunc_max - xfunc_min) / base_segments;
    for (int i = 0; i < base_segments; ++i) {
        double x0 = xfunc_min + i * dx;
        double x1 = x0 + dx;
        double y0 = eval_func('B', x0);
        double y1 = eval_func('B', x1);
        subdiv_draw(w, 'B', x0, x1, y0, y1, 0, 0.0, 0.0);
    }
}

// draw dynamic grid and ticks
static void draw_grid_and_ticks(const plotWidget *w)
{
    // border
    DrawRectangleLinesEx(w->viewport_px, 2, GRAY);

    double xrange = w->xmax - w->xmin;
    double yrange = w->ymax - w->ymin;
    if (xrange <= 0 || yrange <= 0)
        return;

    // compute nice tick steps
    int target_xticks = (int)fmax(4, w->viewport_px.width / 120.0);
    int target_yticks = (int)fmax(4, w->viewport_px.height / 80.0);
    double xstep = nice_tick_step(xrange, target_xticks);
    double ystep = nice_tick_step(yrange, target_yticks);

    // find first tick >= xmin
    double startx = ceil(w->xmin / xstep) * xstep;
    double starty = ceil(w->ymin / ystep) * ystep;

    // colors
    Color grid_col = (Color){200,200,200,80};
    Color x_axis_col = (Color){180,50,50,200};
    Color y_axis_col = (Color){50,180,50,200};
    Color tick_col = (Color){200,200,200,255};
    Color label_col = (Color){180,180,180,255};

    int sx = (int)fmax(0, floor(w->viewport_px.x));
    int sy = (int)fmax(0, floor(w->viewport_px.y));
    int sw = (int)fmax(0, ceil(w->viewport_px.width));
    int sh = (int)fmax(0, ceil(w->viewport_px.height));
    BeginScissorMode(sx, sy, sw, sh); {
        // vertical grid lines + ticks (x axis)
        for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
            Vector2 p1 = world_to_screen(w, tx, w->ymin);
            Vector2 p2 = world_to_screen(w, tx, w->ymax);
            DrawLineEx(p1, p2, 1.0f, grid_col);

            // tick on bottom
            Vector2 tl = (Vector2){ p1.x, p1.y };
            DrawLine((int)tl.x, (int)(tl.y), (int)tl.x, (int)(tl.y - 6), tick_col);
        }
        // horizontal grid lines + ticks (y axis)
        for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
            Vector2 p1 = world_to_screen(w, w->xmin, ty);
            Vector2 p2 = world_to_screen(w, w->xmax, ty);
            DrawLineEx(p1, p2, 1.0f, grid_col);

            // tick on left
            Vector2 tl = (Vector2){ p1.x, p1.y };
            DrawLine((int)(tl.x), (int)tl.y, (int)(tl.x + 6), (int)tl.y, tick_col);
        }

        // draw axis x=0 e y=0, if visible
        if (w->xmin <= 0 && w->xmax >= 0) {
            Vector2 a = world_to_screen(w, 0.0, w->ymin);
            Vector2 b = world_to_screen(w, 0.0, w->ymax);
            DrawLineEx(a, b, 1.0f, y_axis_col);
        }
        if (w->ymin <= 0 && w->ymax >= 0) {
            Vector2 a = world_to_screen(w, w->xmin, 0.0);
            Vector2 b = world_to_screen(w, w->xmax, 0.0);
            DrawLineEx(a, b, 1.0f, x_axis_col);
        }

        // world area info (lower right corner)
        char buf[128];
        snprintf(buf, sizeof(buf), "x:[%.3g, %.3g] y:[%.3g, %.3g]", w->xmin, w->xmax, w->ymin, w->ymax);
        int tw = MeasureText(buf, 10);
        DrawText(buf, (int)(w->viewport_px.x + w->viewport_px.width - tw - 6), (int)(w->viewport_px.y + w->viewport_px.height - 26), 10, DARKGRAY);

    } EndScissorMode();

    // horizontal ticks label
    for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
        Vector2 s0 = world_to_screen(w, tx, w->ymin);
        char buf[64];
        if (fabs(tx) < 1e-3 || fabs(tx) > 1e4) {
            snprintf(buf, sizeof(buf), "%.3g", tx);
        } else {
            // scegli precisione in base a xstep
            int prec = (int)ceil(-floor(log10(xstep)));
            if (prec < 0) prec = 0;
            char fmt[16]; snprintf(fmt, sizeof(fmt), "%%.%df", prec);
            snprintf(buf, sizeof(buf), fmt, tx);
        }
        Vector2 lblpos = (Vector2){ s0.x - MeasureText(buf, 10) / 2.0f, w->viewport_px.y + w->viewport_px.height + 2 };
        DrawText(buf, (int)lblpos.x, (int)lblpos.y, 10, label_col);
    }

    // vertical ticks label
    for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
        Vector2 s0 = world_to_screen(w, w->xmin, ty);
        char buf[64];
        if (fabs(ty) < 1e-3 || fabs(ty) > 1e4) {
            snprintf(buf, sizeof(buf), "%.3g", ty);
        } else {
            int prec = (int)ceil(-floor(log10(ystep)));
            if (prec < 0) prec = 0;
            char fmt[16]; snprintf(fmt, sizeof(fmt), "%%.%df", prec);
            snprintf(buf, sizeof(buf), fmt, ty);
        }
        Vector2 lblpos = (Vector2){ w->viewport_px.x - MeasureText(buf, 10) - 6, s0.y - 8 };
        DrawText(buf, (int)lblpos.x, (int)lblpos.y, 10, label_col);
    }
}

// draw curve data (series) with viewport clipping
static void draw_curve(const plotWidget *w, const pointD *pts, size_t n, Color color)
{
    if (!pts || n < 2) return;

    // set scissor/clipping to viewport
    int sx = (int)fmax(0, floor(w->viewport_px.x));
    int sy = (int)fmax(0, floor(w->viewport_px.y));
    int sw = (int)fmax(0, ceil(w->viewport_px.width));
    int sh = (int)fmax(0, ceil(w->viewport_px.height));
    BeginScissorMode(sx, sy, sw, sh); {

        // draw segments; skip NaNs
        for (size_t i = 0; i + 1 < n; ++i) {
            pointD a = pts[i];
            pointD b = pts[i+1];
            if (!isfinite(a.x) || !isfinite(a.y) || !isfinite(b.x) || !isfinite(b.y)) continue;
            // trivial visibility check: if both outside on one side, skip
            if ((a.x < w->xmin && b.x < w->xmin) || (a.x > w->xmax && b.x > w->xmax) ||
                (a.y < w->ymin && b.y < w->ymin) || (a.y > w->ymax && b.y > w->ymax)) continue;
            Vector2 sa = world_to_screen(w, a.x, a.y);
            Vector2 sb = world_to_screen(w, b.x, b.y);
            DrawLineV(sa, sb, color);
        }

    } EndScissorMode();
}

// wrapper to draw a named series (with legend handling elsewhere)
static void draw_series(const plotWidget *w, const series_t *s)
{
    if (!s->visible)
        return;
    draw_curve(w, s->pts, s->n, s->color);
}

// handle inputs: pan & zoom
static void handle_input(plotWidget *w)
{
    // Pan with middle mouse drag
    Vector2 mpos = GetMousePosition();
    if (IsMouseButtonPressed(MOUSE_MIDDLE_BUTTON)) {
        w->dragging = true;
        w->last_mouse = mpos;
    }
    if (IsMouseButtonReleased(MOUSE_MIDDLE_BUTTON)) {
        w->dragging = false;
    }
    if (w->dragging) {
        Vector2 d = Vector2Subtract(mpos, w->last_mouse);
        w->last_mouse = mpos;
        double dx = -d.x / w->viewport_px.width * (w->xmax - w->xmin);
        double dy = d.y / w->viewport_px.height * (w->ymax - w->ymin);
        w->xmin += dx; w->xmax += dx;
        w->ymin += dy; w->ymax += dy;
    }

    // Zoom with mouse wheel centered on mouse
    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f) {
        Vector2 mouse = GetMousePosition();
        if (mouse.x >= w->viewport_px.x && mouse.x <= w->viewport_px.x + w->viewport_px.width &&
            mouse.y >= w->viewport_px.y && mouse.y <= w->viewport_px.y + w->viewport_px.height) {
            double mx, my;
            screen_to_world(w, mouse.x, mouse.y, &mx, &my);
            double k = pow(1.15, -wheel);
            double nxmin = mx + (w->xmin - mx) * k;
            double nxmax = mx + (w->xmax - mx) * k;
            double nymin = my + (w->ymin - my) * k;
            double nymax = my + (w->ymax - my) * k;
            // avoid degenerate
            if (nxmax - nxmin > 1e-6 && nymax - nymin > 1e-6) {
                w->xmin = nxmin; w->xmax = nxmax;
                w->ymin = nymin; w->ymax = nymax;
            }
        }
    }

    // Area selection: Left mouse drag inside viewport
    mpos = GetMousePosition();
    bool inside = (mpos.x >= w->viewport_px.x && mpos.x <= w->viewport_px.x + w->viewport_px.width &&
                   mpos.y >= w->viewport_px.y && mpos.y <= w->viewport_px.y + w->viewport_px.height);

    if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && inside) {
        // start area drag: only if not clicking on legend (legend handler should run before)
        g_area.active = true;
        g_area.finished = false;
        g_area.start_px = mpos;
        g_area.end_px = mpos;
        // clear single selection while starting area selection
        g_selection.active = false;
    }

    if (g_area.active) {
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            g_area.end_px = mpos;
        }
        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON)) {
            g_area.active = false;
            g_area.finished = true;
            g_area.end_px = mpos;
            // After finishing, compute selected points
            // free previous
            if (g_multi.series_idx) { free(g_multi.series_idx); g_multi.series_idx = NULL; }
            if (g_multi.pt_idx) { free(g_multi.pt_idx); g_multi.pt_idx = NULL; }
            g_multi.count = 0;
            g_multi.has_any = false;

            // compute world bounding box of selection
            Rectangle srect = norm_rect_from_points(g_area.start_px, g_area.end_px);
            double wxmin, wxmax, wymin, wymax;
            rect_screen_to_world_box(w, srect, &wxmin, &wxmax, &wymin, &wymax);

            // collect points inside the world bbox
            // first pass count
            for (size_t si = 0; si < nseries; ++si) {
                if (!series[si].visible) continue;
                for (size_t pi = 0; pi < series[si].n; ++pi) {
                    pointD p = series[si].pts[pi];
                    if (!isfinite(p.x) || !isfinite(p.y)) continue;
                    if (p.x >= wxmin && p.x <= wxmax && p.y >= wymin && p.y <= wymax) {
                        g_multi.count++;
                    }
                }
            }
            if (g_multi.count > 0) {
                g_multi.series_idx = (size_t*)malloc(sizeof(size_t) * g_multi.count);
                g_multi.pt_idx = (size_t*)malloc(sizeof(size_t) * g_multi.count);
                size_t idx = 0;
                for (size_t si = 0; si < nseries; ++si) {
                    if (!series[si].visible) continue;
                    for (size_t pi = 0; pi < series[si].n; ++pi) {
                        pointD p = series[si].pts[pi];
                        if (!isfinite(p.x) || !isfinite(p.y)) continue;
                        if (p.x >= wxmin && p.x <= wxmax && p.y >= wymin && p.y <= wymax) {
                            g_multi.series_idx[idx] = si;
                            g_multi.pt_idx[idx] = pi;
                            idx++;
                        }
                    }
                }
                g_multi.has_any = true;
            } else {
                g_multi.has_any = false;
            }
        }
    }

    // Cancel area selection with Escape
    if (IsKeyPressed(KEY_Q)) {
        g_area.active = false;
        g_area.finished = false;
    }

    // Reset view with R
    if (IsKeyPressed(KEY_R)) {
        w->xmin = -10.0; w->xmax = 10.0; w->ymin = -3.0; w->ymax = 3.0;
    }


}

// draw widget's frame UI (frame, title)
static void draw_widget_frame(const plotWidget *w, const char *title)
{
    // border already drawn; draw title at top center of viewport
    int tw = MeasureText(title, 12);
    DrawText(title, (int)(w->viewport_px.x + (w->viewport_px.width - tw) / 2.0), (int)(w->viewport_px.y - 18), 12, LIGHTGRAY);
}

// example: generate points for noisy y = sin(x) + ...
static pointD* generate_example_curve(size_t *out_n)
{
    int n = 800;
    pointD *pts = (pointD*)malloc(sizeof(pointD) * n);
    double x0 = -20.0;
    double x1 = 20.0;
    for (int i = 0; i < n; ++i) {
        double t = (double)i / (n - 1);
        double x = lerp(x0, x1, t);
        double y = sin(x) + 0.2 * sin(3.0 * x) + ((rand() % 1000) / 1000.0 - 0.5) * 0.05;
        pts[i].x = x; pts[i].y = y;
    }
    *out_n = n;
    return pts;
}

// example: generate points for obstacle curve with gaps and different sampling
static pointD* generate_obstacle_curve(size_t *out_n)
{
    int n = 400;
    pointD *pts = (pointD*)malloc(sizeof(pointD) * n);
    double x0 = -10.0;
    double x1 = 30.0;
    for (int i = 0; i < n; ++i) {
        double t = (double)i / (n - 1);
        double x = lerp(x0, x1, t);
        double y;
        if (x > 5.0 && x < 8.0) {
            // gap (NaN)
            y = NAN;
        } else {
            y = 0.5 * cos(0.5 * x) + (rand()%1000)/1000.0*0.05;
        }
        pts[i].x = x; pts[i].y = y;
    }
    *out_n = n;
    return pts;
}

// helper: generate a linear series
static pointD* generate_linear_series(size_t *out_n, double x0, double x1, int samples)
{
    int n = samples;
    pointD *pts = (pointD*)malloc(sizeof(pointD) * n);
    for (int i = 0; i < n; ++i) {
        double t = (double)i / (n - 1);
        double x = lerp(x0, x1, t);
        double y = 0.05 * (x - x0);
        pts[i].x = x; pts[i].y = y;
    }
    *out_n = n;
    return pts;
}

// draw legend box (list of series names with colors). Returns area used.
static Rectangle draw_legend(const plotWidget *w, const series_t *series, size_t nseries)
{
    float pad = 8;
    float entry_h = 20;
    float box_w = 200;
    float box_h = pad*2 + nseries * entry_h;
    float x = w->viewport_px.x + w->viewport_px.width - box_w - 8;
    float y = w->viewport_px.y + 8;

    DrawRectangle((int)x, (int)y, (int)box_w, (int)box_h, (Color){250,250,250,16});
    DrawRectangle((int)x, (int)y, (int)box_w, (int)box_h, (Color){250,250,250,16});
    DrawRectangleLines((int)x, (int)y, (int)box_w, (int)box_h, GRAY);

    for (size_t i = 0; i < nseries; ++i) {
        float ey = y + pad + i * entry_h;
        // color swatch
        DrawRectangle((int)(x + 6), (int)(ey + 4), 12, 12, series[i].color);
        // name
        DrawText(series[i].name, (int)(x + 26), (int)(ey), 12, series[i].visible ? LIGHTGRAY : GRAY);
    }

    return (Rectangle){ x, y, box_w, box_h };
}

// check mouse clicks in legend to toggle series visibility
static void handle_legend_input(const plotWidget *w, series_t *series, size_t nseries, Rectangle legend_box)
{
    if (!IsMouseButtonPressed(MOUSE_LEFT_BUTTON))
        return;
    Vector2 m = GetMousePosition();
    if (m.x < legend_box.x || m.x > legend_box.x + legend_box.width ||
        m.y < legend_box.y || m.y > legend_box.y + legend_box.height)
        return;

    float pad = 8;
    float entry_h = 20;
    for (size_t i = 0; i < nseries; ++i) {
        float ey = legend_box.y + pad + i * entry_h;
        Rectangle entry_r = { legend_box.x, ey, legend_box.width, entry_h };
        if (CheckCollisionPointRec(m, entry_r)) {
            series[i].visible = !series[i].visible;
            break;
        }
    }
}

int main(void)
{
    srand((unsigned int)clock());

    const int screenW = 1000;
    const int screenH = 700;
    InitWindow(screenW, screenH, "Plot Widget - raylib");
    SetTargetFPS(60);

    // widget initialization
    plotWidget widget;
    widget.viewport_px = (Rectangle){ 60, 40, screenW - 120.0f, screenH - 100.0f };
    widget.xmin = -10.0; widget.xmax = 10.0;
    widget.ymin = -3.0; widget.ymax = 3.0;
    widget.dragging = false;

    // create series
    size_t n1, n2, n3;
    pointD *pts1 = generate_example_curve(&n1);
    pointD *pts2 = generate_obstacle_curve(&n2);
    pointD *pts3 = generate_linear_series(&n3, -30.0, 30.0, 300);

    nseries = 3;
    series = (series_t*)malloc(sizeof(series_t) * nseries);

    series[0].pts = pts1; 
    series[0].n = n1; 
    series[0].color = RED; 
    series[0].visible = true; 
    strncpy(series[0].name, "Noisy sin", sizeof(series[0].name)-1);
    series[1].pts = pts2; 
    series[1].n = n2; 
    series[1].color = BLUE; 
    series[1].visible = true; 
    strncpy(series[1].name, "Cos w/ gap", sizeof(series[1].name)-1);
    series[2].pts = pts3; 
    series[2].n = n3; 
    series[2].color = YELLOW; 
    series[2].visible = true; 
    strncpy(series[2].name, "Linear", sizeof(series[2].name)-1);

    // main loop
    while (!WindowShouldClose()) {
        BeginDrawing(); {
            ClearBackground(BLACK);

            // UI background per viewport
            DrawRectangleRec(widget.viewport_px, (Color){245,245,245,32});

            // Draw grid and ticks
            draw_grid_and_ticks(&widget);

            // Draw all series
            for (size_t i = 0; i < nseries; ++i) {
                draw_series(&widget, &series[i]);
            }

            // Draw adaptive analytic function clipped to viewport
            int sx = (int)fmax(0, floor(widget.viewport_px.x));
            int sy = (int)fmax(0, floor(widget.viewport_px.y));
            int sw = (int)fmax(0, ceil(widget.viewport_px.width));
            int sh = (int)fmax(0, ceil(widget.viewport_px.height));
            BeginScissorMode(sx, sy, sw, sh); {
                draw_adaptive_curve(&widget, widget.xmin, widget.xmax);
            } EndScissorMode();

            // Frame and title
            draw_widget_frame(&widget, "Example: multi-series plot (click legend to toggle)");

            // disegna area di drag in corso (wireframe) - dentro lo schermo normale (no scissor necessario ma meglio limitare al viewport)
            if (g_area.active) {
                Rectangle r = norm_rect_from_points(g_area.start_px, g_area.end_px);
                // draw semi-transparent fill
                DrawRectangleRec(r, (Color){ 200, 200, 255, 60 });
                DrawRectangleLines((int)r.x, (int)r.y, (int)r.width, (int)r.height, (Color){100,100,255,180});
            }

            // se area finished, evidenzia punti selezionati
            if (g_area.finished && g_multi.has_any) {
                // draw highlight circles for each selected point (clipped to viewport)
                BeginScissorMode((int)widget.viewport_px.x, (int)widget.viewport_px.y, (int)widget.viewport_px.width, (int)widget.viewport_px.height);
                for (size_t i = 0; i < g_multi.count; ++i) {
                    size_t si = g_multi.series_idx[i];
                    size_t pi = g_multi.pt_idx[i];
                    if (si >= nseries) continue;
                    pointD p = series[si].pts[pi];
                    Vector2 sp = world_to_screen(&widget, p.x, p.y);
                    DrawCircleV(sp, 4, WHITE);
                    DrawCircleV(sp, 2, series[si].color);
                }
                EndScissorMode();
                // optionally draw a small summary box near top-left of viewport
                char buf[256];
                int lines = 0;
                int max_lines = 6;
                int shown = (int)fmin((double)g_multi.count, (double)max_lines);
                snprintf(buf, sizeof(buf), "Selected: %zu pts", g_multi.count);
                DrawText(buf, (int)(widget.viewport_px.x + 6), (int)(widget.viewport_px.y + 6), 12, DARKGRAY);
                for (int i = 0; i < shown; ++i) {
                    size_t si = g_multi.series_idx[i];
                    size_t pi = g_multi.pt_idx[i];
                    pointD p = series[si].pts[pi];
                    snprintf(buf, sizeof(buf), "%s: x=%.6g y=%.6g", series[si].name, p.x, p.y);
                    DrawText(buf, (int)(widget.viewport_px.x + 6), (int)(widget.viewport_px.y + 24 + i * 14), 12, series[si].color);
                }
                if (g_multi.count > (size_t)shown) {
                    snprintf(buf, sizeof(buf), "... +%zu more", g_multi.count - shown);
                    DrawText(buf, (int)(widget.viewport_px.x + 6), (int)(widget.viewport_px.y + 24 + shown * 14), 12, DARKGRAY);
                }
            }

            // Legend (and get rectangle)
            Rectangle legend_box = draw_legend(&widget, series, nseries);
            // handle clicks on legend
            handle_legend_input(&widget, series, nseries, legend_box);

            // input handling (pan/zoom)
            handle_input(&widget);

            Vector2 mpos = GetMousePosition();
            const double pick_threshold_px = 8.0;
            size_t found_si=0, found_pi=0;
            pointD found_pt;
            double found_d;
            bool hit = find_nearest_point(&widget, series, nseries, mpos, pick_threshold_px, &found_si, &found_pi, &found_pt, &found_d);

            // tooltip (transient) si disegna poi con draw_mouse_tooltip
            // selection toggle on left click
            if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                // se click dentro legend l'handler legend lo gestisce prima; quindi qui consideriamo solo viewport
                if (mpos.x >= widget.viewport_px.x && mpos.x <= widget.viewport_px.x + widget.viewport_px.width &&
                    mpos.y >= widget.viewport_px.y && mpos.y <= widget.viewport_px.y + widget.viewport_px.height) {
                    if (hit) {
                        g_selection.active = true;
                        g_selection.series_idx = found_si;
                        g_selection.pt_idx = found_pi;
                        g_selection.pt = found_pt;
                    } else {
                        // click fuori da qualsivoglia punto annulla selezione
                        g_selection.active = false;
                    }
                }
            }
            // cancel selection with 'q'
            if (IsKeyPressed(KEY_Q))
                g_selection.active = false;

            BeginScissorMode(sx, sy, sw, sh); {
                // evidenzia punto selezionato permanente
                if (g_selection.active) {
                    // sicurezza: indice serie valido
                    if (g_selection.series_idx < nseries) {
                        series_t *s = &series[g_selection.series_idx];
                        Vector2 sp = world_to_screen(&widget, g_selection.pt.x, g_selection.pt.y);
                        // cerchio di evidenziazione (fuori scissor per visibilità) o dentro viewport
                        DrawCircleV(sp, 6, WHITE);
                        DrawCircleV(sp, 4, s->color);
                        // testo info vicino al punto (box)
                        char sbuf[128];
                        snprintf(sbuf, sizeof(sbuf), "%s\nx=%.6g\ny=%.6g", s->name, g_selection.pt.x, g_selection.pt.y);
                        int tw = MeasureTextEx(GetFontDefault(), sbuf, 12, 1).x;
                        int th = 48;
                        float bx = sp.x + 10;
                        float by = sp.y - 6;
                        DrawRectangle((int)bx, (int)by, tw + 8, th, (Color){255,255,255,220});
                        DrawText(sbuf, (int)bx + 4, (int)by + 4, 12, DARKGRAY);
                    }
                }

                // tooltip sotto il mouse (transient), disegnalo sempre per l'utente:
                draw_mouse_tooltip(&widget, series, nseries, mpos);
            } EndScissorMode();

            // instructions
            DrawText("Pan: Middle mouse drag    Zoom: Mouse wheel    Reset: R", 10, 10, 12, DARKGRAY);

        } EndDrawing();
    }

    // cleanup
    for (size_t i = 0; i < nseries; ++i) {
        free(series[i].pts);
    }
    free(series);

    CloseWindow();
    return 0;
}
