#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdbool.h>
#include <float.h>

#include <raylib.h>
#define RAYMATH_IMPLEMENTATION
#include <raymath.h>

#define STB_DS_IMPLEMENTATION
#include "stb_ds.h"

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
typedef struct
{
    bool active;          // true se un punto è selezionato
    size_t series_idx;    // indice della serie selezionata
    size_t pt_idx;        // indice del punto selezionato nella serie
    pointD pt;            // coordinate world del punto selezionato
} selection_t;

static selection_t g_selection = {};

// area selection state
typedef struct
{
    bool active;        // true mentre l'utente sta trascinando l'area
    bool finished;      // true se una area è stata selezionata (oltre al drag in corso)
    Vector2 start_px;   // corner iniziale in screen coords
    Vector2 end_px;     // corner corrente/finale in screen coords
} area_sel_t;

static area_sel_t g_area = {};

typedef struct
{
    bool has_any;
    // per semplicità memorizziamo un array dinamico di (series_idx, pt_idx)
    size_t *series_idx;
    size_t *pt_idx;
    size_t count;
} multi_selection_t;

static multi_selection_t g_multi = {};

//size_t nseries = 0;
//series_t *series = NULL;

typedef struct plot_widget_s
{
    Rectangle viewport_px;     // area in pixel sullo schermo dedicata al grafico
    double xmin, xmax;        // vista in coordinate mondo (x)
    double ymin, ymax;        // vista in coordinate mondo (y)
    bool show_grid;
    bool show_ticks;
    bool show_labels;
    bool show_legend;
    bool dragging;
    Vector2 last_mouse;
    char title[128];
    size_t nseries;
    series_t *series;
} plot_widget_t;

// utility functions: world <-> screen transformations
static inline double lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}
static inline double clamp_double(double v, double a, double b)
{
    return v < a ? a : (v > b ? b : v);
}

static inline bool mod_shift(void)
{
    return IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT);
}
static inline bool mod_ctrl(void)
{
    return IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL);
}

static Vector2 world_to_screen(const plot_widget_t *w, double wx, double wy)
{
    double sx = (wx - w->xmin) / (w->xmax - w->xmin) * w->viewport_px.width + w->viewport_px.x;
    double sy = (1.0 - (wy - w->ymin) / (w->ymax - w->ymin)) * w->viewport_px.height + w->viewport_px.y;
    return (Vector2){ (float)sx, (float)sy };
}

static void screen_to_world(const plot_widget_t *w, double sx, double sy, double *wx, double *wy)
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
static void rect_screen_to_world_box(const plot_widget_t *w, Rectangle r, double *out_xmin, double *out_xmax, double *out_ymin, double *out_ymax)
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
static double mouse_point_distance_px(const plot_widget_t *w, const pointD *p, Vector2 mouse)
{
    Vector2 sp = world_to_screen(w, p->x, p->y);
    double dx = sp.x - mouse.x;
    double dy = sp.y - mouse.y;
    return sqrt(dx*dx + dy*dy);
}

// cerca punto più vicino (entro threshold_px) e riempie out_* se trovato
static bool find_nearest_point(const plot_widget_t *w, const series_t *series, size_t nseries, Vector2 mouse, double threshold_px, size_t *out_series, size_t *out_ptidx, pointD *out_pt, double *out_dist_px)
{
//    bool found = false;
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
static void draw_mouse_tooltip(const plot_widget_t *w, const series_t *series, size_t nseries, Vector2 mouse)
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
//    int yoff = (int)(w->viewport_px.y + w->viewport_px.height + 6); // base per disegno sotto viewport
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
static double point_line_distance_px(const plot_widget_t *w, double x0, double y0, double x1, double y1, double xm, double ym)
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
static void subdiv_draw(const plot_widget_t *w, char fmt, double x0, double x1, double y0, double y1, int depth, double extra1, double extra2)
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
static void draw_adaptive_curve(const plot_widget_t *w, double xfunc_min, double xfunc_max)
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
static void draw_grid_and_ticks(const plot_widget_t *w)
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

    if (w->show_grid) {
       int sx = (int)fmax(0, floor(w->viewport_px.x));
        int sy = (int)fmax(0, floor(w->viewport_px.y));
        int sw = (int)fmax(0, ceil(w->viewport_px.width));
        int sh = (int)fmax(0, ceil(w->viewport_px.height));
        BeginScissorMode(sx, sy, sw, sh); {
            // vertical grid lines (x axis)
            for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
                Vector2 p1 = world_to_screen(w, tx, w->ymin);
                Vector2 p2 = world_to_screen(w, tx, w->ymax);
                DrawLineEx(p1, p2, 1.0f, grid_col);
            }
            // horizontal grid lines (y axis)
            for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
                Vector2 p1 = world_to_screen(w, w->xmin, ty);
                Vector2 p2 = world_to_screen(w, w->xmax, ty);
                DrawLineEx(p1, p2, 1.0f, grid_col);
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
    }

    if (w->show_ticks) {
        // vertical ticks (x axis)
        for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
            Vector2 p1 = world_to_screen(w, tx, w->ymin);
            Vector2 tl = (Vector2){ p1.x, p1.y };
            DrawLine((int)tl.x, (int)(tl.y), (int)tl.x, (int)tl.y, tick_col);
        }
        // horizontal ticks (y axis)
        for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
            Vector2 p1 = world_to_screen(w, w->xmin, ty);

            // tick on left
            Vector2 tl = (Vector2){ p1.x, p1.y };
            DrawLine((int)(tl.x), (int)tl.y, (int)(tl.x), (int)tl.y, tick_col);
        }
    }

    if (w->show_labels) {
        // vertical ticks labels (x axis)
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
            Vector2 lblpos = (Vector2){ s0.x - MeasureText(buf, 10) / 2.0f, w->viewport_px.y + w->viewport_px.height };
            DrawText(buf, (int)lblpos.x, (int)lblpos.y, 10, label_col);
        }
        // horizontal ticks labels (y axis)
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
            Vector2 lblpos = (Vector2){ w->viewport_px.x - MeasureText(buf, 10), s0.y - 8 };
            DrawText(buf, (int)lblpos.x, (int)lblpos.y, 10, label_col);
        }
    }
}
// draw curve data (series) with viewport clipping
static void draw_curve(const plot_widget_t *w, const pointD *pts, size_t n, Color color)
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
static void draw_series(const plot_widget_t *w, const series_t *s)
{
    if (!s->visible)
        return;
    draw_curve(w, s->pts, s->n, s->color);
}

void get_series_ranges(const series_t *s, int ns, double *xmin, double *xmax, double *ymin, double *ymax)
{
    *xmin = FLT_MAX; *xmax = FLT_MIN;
    *ymin = FLT_MAX; *ymax = FLT_MIN;
    for (int si = 0; si < ns; ++si) {
        pointD *pts = s[si].pts;
        size_t n = s[si].n;
        for (int pi = 0; pi < n; ++pi) {
            *xmin = fmin(*xmin, pts[pi].x);
            *xmax = fmax(*xmax, pts[pi].x);
            *ymin = fmin(*ymin, pts[pi].y);
            *ymax = fmax(*ymax, pts[pi].y);
        }
    }
}

// handle inputs: pan & zoom
static void handle_input(plot_widget_t *w)
{
    // Pan with middle mouse drag, with axis modifiers
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

        // world del delta (full)
        double dx_full = -d.x / w->viewport_px.width * (w->xmax - w->xmin);
        double dy_full =  d.y / w->viewport_px.height * (w->ymax - w->ymin);

        bool shift = mod_shift();
        bool ctrl  = mod_ctrl();

        if (shift && !ctrl) {
            // Shift -> pan only X
            w->xmin += dx_full; w->xmax += dx_full;
        } else if (ctrl && !shift) {
            // Ctrl -> pan only Y
            w->ymin += dy_full; w->ymax += dy_full;
        } else {
            // No modifier or both -> pan both axes
            w->xmin += dx_full; w->xmax += dx_full;
            w->ymin += dy_full; w->ymax += dy_full;
        }
    }

    // Zoom with mouse wheel centered on mouse, with axis modifiers
    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f) {
        Vector2 mouse = GetMousePosition();
        if (mouse.x >= w->viewport_px.x && mouse.x <= w->viewport_px.x + w->viewport_px.width &&
            mouse.y >= w->viewport_px.y && mouse.y <= w->viewport_px.y + w->viewport_px.height) {

            double mx, my;
            screen_to_world(w, mouse.x, mouse.y, &mx, &my);
            double k = pow(1.15, -wheel);

            bool shift = mod_shift();
            bool ctrl  = mod_ctrl();

            if (shift && !ctrl) {
                // Shift -> zoom only X (change xmin/xmax around mx)
                double nxmin = mx + (w->xmin - mx) * k;
                double nxmax = mx + (w->xmax - mx) * k;
                if (nxmax - nxmin > 1e-12) {
                    w->xmin = nxmin; w->xmax = nxmax;
                }
            } else if (ctrl && !shift) {
                // Ctrl -> zoom only Y (change ymin/ymax around my)
                double nymin = my + (w->ymin - my) * k;
                double nymax = my + (w->ymax - my) * k;
                if (nymax - nymin > 1e-12) {
                    w->ymin = nymin; w->ymax = nymax;
                }
            } else {
                // No modifier or both -> zoom both axes
                double nxmin = mx + (w->xmin - mx) * k;
                double nxmax = mx + (w->xmax - mx) * k;
                double nymin = my + (w->ymin - my) * k;
                double nymax = my + (w->ymax - my) * k;
                if (nxmax - nxmin > 1e-12 && nymax - nymin > 1e-12) {
                    w->xmin = nxmin; w->xmax = nxmax;
                    w->ymin = nymin; w->ymax = nymax;
                }
            }
        }
    }

    // Area selection: Left mouse drag inside viewport
    mpos = GetMousePosition();
    bool inside = (mpos.x >= w->viewport_px.x &&
                   mpos.x <= w->viewport_px.x + w->viewport_px.width &&
                   mpos.y >= w->viewport_px.y &&
                   mpos.y <= w->viewport_px.y + w->viewport_px.height);

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
            if (g_multi.series_idx) {
                free(g_multi.series_idx);
                g_multi.series_idx = NULL;
            }
            if (g_multi.pt_idx) {
                free(g_multi.pt_idx);
                g_multi.pt_idx = NULL;
            }
            g_multi.count = 0;
            g_multi.has_any = false;

            // compute world bounding box of selection
            Rectangle srect = norm_rect_from_points(g_area.start_px, g_area.end_px);
            double wxmin, wxmax, wymin, wymax;
            rect_screen_to_world_box(w, srect, &wxmin, &wxmax, &wymin, &wymax);

            // collect points inside the world bbox
            // first pass count
            for (size_t si = 0; si < w->nseries; ++si) {
                if (!w->series[si].visible) continue;
                for (size_t pi = 0; pi < w->series[si].n; ++pi) {
                    pointD p = w->series[si].pts[pi];
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
                for (size_t si = 0; si < w->nseries; ++si) {
                    if (!w->series[si].visible) continue;
                    for (size_t pi = 0; pi < w->series[si].n; ++pi) {
                        pointD p = w->series[si].pts[pi];
                        if (!isfinite(p.x) || !isfinite(p.y))
                            continue;
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
        w->xmin = -10.0; w->xmax = 10.0;
        w->ymin =  -3.0; w->ymax =  3.0;
    }

    // fill to series
    if (IsKeyPressed(KEY_F)) {
        double xmin, xmax;
        double ymin, ymax;
        get_series_ranges(w->series, w->nseries, &xmin, &xmax, &ymin, &ymax);
        w->xmin = xmin - (xmax - xmin) * 0.05;
        w->xmax = xmax + (xmax - xmin) * 0.05;
        w->ymin = ymin - (ymax - ymin) * 0.05;
        w->ymax = ymax + (ymax - ymin) * 0.05;
    }
}

// draw widget's frame UI (frame, title)
static void draw_widget_frame(const plot_widget_t *w)
{
    // border already drawn; draw title at top center of viewport
    int tw = MeasureText(w->title, 12);
    DrawText(w->title, (int)(w->viewport_px.x + (w->viewport_px.width - tw) / 2.0), (int)(w->viewport_px.y - 18), 12, LIGHTGRAY);
}

typedef enum
{
    IZK_MODEL_TONIC_SPIKING,
    IZK_MODEL_PHASIC_SPIKING,
    IZK_MODEL_TONIC_BURSTING,
    IZK_MODEL_PHASIC_BURSTING,
    IZK_MODEL_MIXED_MODE,
    IZK_MODEL_SPIKE_FREQUENCY_ADAPTATION,
    IZK_MODEL_CLASS_1,
    IZK_MODEL_CLASS_2,
    IZK_MODEL_SPIKE_LATENCY,
    IZK_MODEL_SUBTHRESHOLD_OSCILLATIONS,
    IZK_MODEL_RESONATOR,
    IZK_MODEL_INTEGRATOR,
    IZK_MODEL_REBOUND_SPIKE,
    IZK_MODEL_REBOUND_BURST,
    IZK_MODEL_THRESHOLD_VARIABILITY,
    IZK_MODEL_BISTABILITY,
    IZK_MODEL_DAP,
    IZK_MODEL_ACCOMODATION,
    IZK_MODEL_INHIBITION_INDUCED_SPIKING,
    IZK_MODEL_INHIBITION_INDUCED_BURSTING,
    IZK_MODEL_N
} izk_model_e;

typedef struct
{
    double a;
    double b;
    double c;
    double d;
    double I;
} izk_param;

izk_param izk_pars[] = {
        {  0.02,      0.2,     -65,     6.0,    14.0 }, // (A) tonic spiking
        {  0.02,      0.25,    -65,     6.0,     0.5 }, // (B) phasic spiking
        {  0.02,      0.2,     -50,     2.0,    15.0 }, // (C) tonic bursting
        {  0.02,      0.25,    -55,     0.05,    0.6 }, // (D) phasic bursting
        {  0.02,      0.2,     -55,     4.0,    10.0 }, // (E) mixed mode
        {  0.01,      0.2,     -65,     8.0,    30.0 }, // (F) spike frequency adaptation
        {  0.02,     -0.1,     -55,     6.0,     0.0 }, // (G) Class 1
        {  0.2,       0.26,    -65,     0.0,     0.0 }, // (H) Class 2
        {  0.02,      0.2,     -65,     6.0,     7.0 }, // (I) spike latency
        {  0.05,      0.26,    -60,     0.0,     0.0 }, // (J) subthreshold oscillations
        {  0.1,       0.26,    -60,    -1.0,     0.0 }, // (K) resonator
        {  0.02,     -0.1,     -55,     6.0,     0.0 }, // (L) integrator
        {  0.03,      0.25,    -60,     4.0,     0.0 }, // (M) rebound spike
        {  0.03,      0.25,    -52,     0.0,     0.0 }, // (N) rebound burst
        {  0.03,      0.25,    -60,     4.0,     0.0 }, // (O) threshold variability
        {  1.0,       1.5,     -60,     0.0,   -65.0 }, // (P) bistability
        {  1.0,       0.2,     -60,   -21.0,     0.0 }, // (Q) DAP
        {  0.02,      1.0,     -55,     4.0,     0.0 }, // (R) accomodation
        { -0.02,     -1.0,     -60,     8.0,    80.0 }, // (S) inhibition-induced spiking
        { -0.026,    -1.0,     -45,     0.0,    80.0 }  // (T) inhibition-induced bursting
};

typedef struct
{
    double V;
    double u;
    double tau;
    double t_start;
    double t_end;
    double T1;
    double T2;
    double T3;
    double T4;
} izk_init_param;

izk_init_param izk_inits[] = {
        { -70.0,   0.0,  0.25,  0.0,  100.0,   0.0,   0, 0, 0 }, // (A)
        { -64.0,   0.0,  0.25,  0.0,  200.0,  20.0,   0, 0, 0 }, // (B)
        { -70.0,   0.0,  0.25,  0.0,  220.0,  22.0,   0, 0, 0 }, // (C)
        { -64.0,   0.0,  0.2,   0.0,  200.0,  20.0,   0, 0, 0 }, // (D)
        { -70.0,   0.0,  0.25,  0.0,  160.0,   0.0,   0, 0, 0 }, // (E)
        { -70.0,   0.0,  0.25,  0.0,   85.0,   0.0,   0, 0, 0 }, // (F)
        { -60.0,   0.0,  0.25,  0.0,  300.0,  30.0,   0, 0, 0 }, // (G)
        { -64.0,   0.0,  0.25,  0.0,  300.0,  30.0,   0, 0, 0 }, // (H)
        { -70.0,   0.0,  0.2,   0.0,  100.0,   0.0,   0, 0, 0 }, // (I)
        { -62.0,   0.0,  0.25,  0.0,  200.0,   0.0,   0, 0, 0 }, // (J)
        { -62.0,   0.0,  0.25,  0.0,  400.0,   0.0,   0, 0, 0 }, // (K)
        { -60.0,   0.0,  0.25,  0.0,  100.0,   0.0,   0, 0, 0 }, // (L)
        { -64.0,   0.0,  0.2,   0.0,  200.0,  20.0,   0, 0, 0 }, // (M)
        { -64.0,   0.0,  0.2,   0.0,  200.0,  20.0,   0, 0, 0 }, // (N)
        { -64.0,   0.0,  0.25,  0.0,  100.0,   0.0,   0, 0, 0 }, // (O)
        { -61.0,   0.0,  0.25,  0.0,  300.0,   0.0, 216, 0, 0 }, // (P)
        { -70.0,   0.0,  0.1,   0.0,   50.0,  10.0,   0, 0, 0 }, // (Q)
        { -65.0, -16.0,  0.5,   0.0,  400.0,   0.0,   0, 0, 0 }, // (R)
        { -63.8,   0.0,  0.5,   0.0,  350.0,   0.0,   0, 0, 0 }, // (S)
        { -63.8,   0.0,  0.5,   0.0,  350.0,   0.0,   0, 0, 0 }, // (T)
};

static pointD* generate_izk_fig1(size_t *out_n, pointD *out_v[3], izk_model_e model)
{
}

// %%%%%%%%%%%%%%% (A) tonic spiking %%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_tonic_spike(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.2;
    const double c = -65.0;
    const double d = 6.0;

    double V = -70.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 10.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 14.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%%% (B) phasic spiking %%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_phasic_spike(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.25;
    const double c = -65.0;
    const double d = 6.0;

    double V = -64.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 20.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 0.5 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (C) tonic bursting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_tonic_bursting(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.2;
    const double c = -50.0;
    const double d = 2.0;

    double V = -70.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 220.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 22.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 15.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%% (D) phasic bursting %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_phasic_bursting(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.25;
    const double c = -55.0;
    const double d = 0.05;

    double V = -64.0;
    double u = b * V;

    const double tau = 0.2;
    const double t_start = 0.0;
    const double t_end = 200.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 20.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 0.6 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%% (E) mixed mode %%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_mixed_mode(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.2;
    const double c = -55.0;
    const double d = 4.0;

    double V = -70.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 160.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 10.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 10.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%% (F) spike freq. adapt %%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_spike_freq_adapt(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.01;
    const double b = 0.2;
    const double c = -65.0;
    const double d = 8.0;

    double V = -70.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 85.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 10.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? 30.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (G) Class 1 exc. %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_class1_exc(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = -0.1;
    const double c = -55.0;
    const double d = 6.0;

    double V = -60.0;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 300.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 30.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? (0.075*(t-T1)) : 0.0;

        V = V + tau * (0.04 * V * V + 4.1 * V + 108 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%%% (H) Class 2 exc. %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_class2_exc(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.2;
    const double b = 0.26;
    const double c = -65;
    const double d = 0;

    double V = -64;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 300.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 30.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = (t > T1) ? -0.5+(0.015*(t-T1)) : -0.5;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (I) spike latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_spike_latency(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 0.2;
    const double c = -65;
    const double d = 6;

    double V = -70;
    double u = b * V;

    const double tau = 0.2;
    const double t_start = 0.0;
    const double t_end = 100.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 10.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t > T1) & (t < T1 + 3.0)) ? 7.04 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (J) subthresh. osc. %%%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_subthr_osc(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.05;
    const double b = 0.26;
    const double c = -60;
    const double d = 0;

    double V = -62;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 200;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 10.0; /* soglia tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t>T1) & (t < T1+5)) ? 2.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%%% (K) resonator %%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_resonator(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.1;
    const double b = 0.26;
    const double c = -60;
    const double d = -1;

    double V = -62;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 400;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* soglie tempo per I */
    const double T1 = t_end / 10.0;
    const double T2 = T1+20;
    const double T3 = 0.7 * t_end;
    const double T4 = T3 + 40;

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t>T1) & (t < T1+4)) | ((t>T2) & (t < T2+4)) | ((t>T3) & (t < T3+4)) | ((t>T4) & (t < T4+4))
                ? 0.65 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%% (L) integrator %%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_integrator(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = -0.1;
    const double c = -55;
    const double d = 6;

    double V = -60;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 100;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* soglie tempo per I */
    const double T1 = t_end / 11.0;
    const double T2 = T1 + 5;
    const double T3 = 0.7 * t_end;
    const double T4 = T3 + 10;

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t>T1) & (t < T1+2)) | ((t>T2) & (t < T2+2)) | ((t>T3) & (t < T3+2)) | ((t>T4) & (t < T4+2))
                ? 9 : 0.0;

        V = V + tau * (0.04 * V * V + 4.1 * V + 108 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (M) rebound spike %%%%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_rebound_spike(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.03;
    const double b = 0.25;
    const double c = -60;
    const double d = 4;

    double V = -64;
    double u = b * V;

    const double tau = 0.2;
    const double t_start = 0.0;
    const double t_end = 200;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 20; /* soglie tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t>T1) & (t < T1+5)) ? -15 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (N) rebound burst %%%%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_rebound_burst(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.03;
    const double b = 0.25;
    const double c = -52;
    const double d = 0;

    double V = -64;
    double u = b * V;

    const double tau = 0.2;
    const double t_start = 0.0;
    const double t_end = 200;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 20; /* soglie tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = ((t>T1) & (t < T1+5)) ? -15 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%%%%% (O) thresh. variability %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_thresh_variability(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.03;
    const double b = 0.25;
    const double c = -60;
    const double d = 4;

    double V = -64;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 100;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I =
                ((t>10) & (t < 15)) | ((t>80) & (t < 85)) ?
                        1 :
                        (t>70) & (t < 75) ?
                                -6:
                                0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (P) bistability %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_bistability(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.1;
    const double b = 0.26;
    const double c = -60;
    const double d = 0;

    double V = -61;
    double u = b * V;

    const double tau = 0.25;
    const double t_start = 0.0;
    const double t_end = 300;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = t_end / 8; /* soglie tempo per I */
    const double T2 = 216;

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I =((t>T1) & (t < T1+5)) | ((t>T2) & (t < T2+5)) ? 1.24 : 0.24;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (Q) DAP %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_DAP(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 1.0;
    const double b = 0.2;
    const double c = -60;
    const double d = -21;

    double V = -70;
    double u = b * V;

    const double tau = 0.1;
    const double t_start = 0.0;
    const double t_end = 50;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    const double T1 = 10; /* soglie tempo per I */

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I =(abs(t-T1)<1.0) ? 20.0 : 0.0;

        V = V + tau * (0.04 * V * V + 5 * V + 140 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (R) accomodation %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_accomodation(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = 0.02;
    const double b = 1;
    const double c = -55.0;
    const double d = 4.0;

    double V = -65.0;
    double u = 16;

    const double tau = 0.5;
    const double t_start = 0.0;
    const double t_end = 400.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    double *II = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = 0;
        if (t < 200)
            I = t / 25.0;
        else if (t < 300.0)
            I = 0.0;
        else if (t < 312.5)
            I = (t - 300.0) / 12.5 * 4.0;
        else
            I = 0.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * (V + 65.0));

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;
        II[idx] = I;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (S) inhibition induced spiking %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_inh_induced_sp(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = -0.02;
    const double b = -1;
    const double c = -60.0;
    const double d = 8.0;

    double V = -63.8;
    double u = b * V;

    const double tau = 0.5;
    const double t_start = 0.0;
    const double t_end = 350.0;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = 0.0;
        if ((t < 50.0) | (t > 250.0))
            I = 80.0;
        else
            I = 75.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
}

// %%%%%%%%%%%%%% (T) inhibition induced bursting %%%%%%%%%%%%%%%%%%%%%%%%%%
static pointD* generate_izk_inh_induced_brst(size_t *out_n, pointD *out_v[3])
{
    /* Parametri */
    const double a = -0.026;
    const double b = -1;
    const double c = -45;
    const double d = -2;

    double V = -63.8;
    double u = b * V;

    const double tau = 0.5;
    const double t_start = 0.0;
    const double t_end = 350;
    const int steps = (int)round((t_end - t_start) / tau) + 1;

    out_v[0] = (pointD*)malloc(sizeof(pointD) * steps); // I
    out_v[1] = (pointD*)malloc(sizeof(pointD) * steps); // V
    out_v[2] = (pointD*)malloc(sizeof(pointD) * steps); // u

    /* Array dinamici per salvare risultati */
    double *tspan = malloc(sizeof(double) * steps);
    double *VV = malloc(sizeof(double) * steps);
    double *uu = malloc(sizeof(double) * steps);
    if (!tspan || !VV || !uu) {
        fprintf(stderr, "Memoria insufficiente\n");
        return NULL;
    }

    /* Simulazione forward Euler come nel MATLAB */
    int idx = 0;
    for (int i = 0; i < steps; ++i) {
        double t = t_start + i * tau;
        tspan[idx] = t;

        double I = 0.0;
        if ((t < 50.0) | (t > 250.0))
            I = 80.0;
        else
            I = 75.0;

        V = V + tau * (0.04 * V * V + 5.0 * V + 140.0 - u + I);
        u = u + tau * a * (b * V - u);

        if (V > 30.0) {
            VV[idx] = 30.0;
            V = c;
            u = u + d;
        } else {
            VV[idx] = V;
        }
        uu[idx] = u;

        out_v[0][i] = (pointD) {i, I};
        out_v[1][i] = (pointD) {i, V};
        out_v[2][i] = (pointD) {i, u};

        ++idx;
    }

    free(tspan);
    free(VV);
    free(uu);

    *out_n = steps;
    return out_v[1]; // V
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
static Rectangle draw_legend(const plot_widget_t *w, const series_t *series, size_t nseries)
{
    float pad = 8;
    float entry_h = 20;
    float box_w = 200;
    float box_h = pad*2 + nseries * entry_h;
    float x = w->viewport_px.x + w->viewport_px.width - box_w - 8;
    float y = w->viewport_px.y + 8;

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
static void handle_legend_input(const plot_widget_t *w, series_t *series, size_t nseries, Rectangle legend_box)
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

    const int screenW = 1150;
    const int screenH = 750;
    InitWindow(screenW, screenH, "Plot Widget - raylib");
    SetTargetFPS(60);

    // create widgets
    plot_widget_t *widgets = NULL;
    size_t widgets_n = 0;
    series_t s = {};

    // widget initialization

    float ww = 200; // screenW - 50; // width
    float wh = 150; // screenH - 50; // height
    int bwx = 25; // base x
    int bwy = 25; // base y
    int wsx = 25; // widget spacing x
    int wsy = 35; // widget spacing y
    int dwx = ww + wsx; // delta x
    int dwy = wh + wsy; // delta y
    int wx = bwx; // widget x
    int wy = bwy; // widget y

#if 0
    // %%%%%%%%%%%%%%% (A) tonic spiking %%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t tonic_spike =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(A) tonic spiking",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, tonic_spike);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t tonic_spike_n;
    pointD *tonic_spike_out_v[3] = {};
    generate_izk_tonic_bursting(&tonic_spike_n, tonic_spike_out_v);

    s = (series_t) { pts: tonic_spike_out_v[0], n: tonic_spike_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_spike_out_v[1], n: tonic_spike_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_spike_out_v[2], n: tonic_spike_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#else
    // %%%%%%%%%%%%%%% (A) tonic spiking %%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t tonic_spike =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(A) tonic spiking",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, tonic_spike);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t tonic_spike_n;
    pointD *tonic_spike_out_v[3] = {};
    generate_izk_tonic_spike(&tonic_spike_n, tonic_spike_out_v);

    s = (series_t) { pts: tonic_spike_out_v[0], n: tonic_spike_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_spike_out_v[1], n: tonic_spike_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_spike_out_v[2], n: tonic_spike_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%%% (B) phasic spiking %%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t phasic_spike =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(B) phasic spiking",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, phasic_spike);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t phasic_spike_n;
    pointD *phasic_spike_out_v[3] = {};
    generate_izk_phasic_spike(&phasic_spike_n, phasic_spike_out_v);

    s = (series_t) { pts: phasic_spike_out_v[0], n: phasic_spike_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: phasic_spike_out_v[1], n: phasic_spike_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: phasic_spike_out_v[2], n: phasic_spike_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%% (C) tonic bursting %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t tonic_burst =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(C) tonic bursting",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, tonic_burst);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t tonic_burst_n;
    pointD *tonic_burst_out_v[3] = {};
    generate_izk_tonic_bursting(&tonic_burst_n, tonic_burst_out_v);

    s = (series_t) { pts: tonic_burst_out_v[0], n: tonic_burst_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_burst_out_v[1], n: tonic_burst_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: tonic_burst_out_v[2], n: tonic_burst_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%% (D) phasic bursting %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t phasic_burst =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(D) phasic bursting",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, phasic_burst);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t phasic_burst_n;
    pointD *phasic_burst_out_v[3] = {};
    generate_izk_phasic_bursting(&phasic_burst_n, phasic_burst_out_v);

    s = (series_t) { pts: phasic_burst_out_v[0], n: phasic_burst_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: phasic_burst_out_v[1], n: phasic_burst_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: phasic_burst_out_v[2], n: phasic_burst_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%% (E) mixed mode %%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t mixed_mode =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(E) mixed mode",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, mixed_mode);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t mixed_mode_n;
    pointD *mixed_mode_out_v[3] = {};
    generate_izk_mixed_mode(&mixed_mode_n, mixed_mode_out_v);

    s = (series_t) { pts: mixed_mode_out_v[0], n: mixed_mode_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: mixed_mode_out_v[1], n: mixed_mode_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: mixed_mode_out_v[2], n: mixed_mode_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx = bwx;
    wy += dwy;

    // %%%%%%%%%%%%%%%% (F) spike freq. adapt %%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t spike_freq_adapt =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(F) spike freq. adapt",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, spike_freq_adapt);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t spike_freq_adapt_n;
    pointD *spike_freq_adapt_out_v[3] = {};
    generate_izk_spike_freq_adapt(&spike_freq_adapt_n, spike_freq_adapt_out_v);

    s = (series_t) { pts: spike_freq_adapt_out_v[0], n: spike_freq_adapt_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: spike_freq_adapt_out_v[1], n: spike_freq_adapt_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: spike_freq_adapt_out_v[2], n: spike_freq_adapt_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (G) Class 1 exc. %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t class1_exc =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(G) Class 1 exc.",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, class1_exc);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t class1_exc_n;
    pointD *class1_exc_out_v[3] = {};
    generate_izk_class1_exc(&class1_exc_n, class1_exc_out_v);

    s = (series_t) { pts: class1_exc_out_v[0], n: class1_exc_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: class1_exc_out_v[1], n: class1_exc_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: class1_exc_out_v[2], n: class1_exc_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%%% (H) Class 2 exc. %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t class2_exc =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(H) Class 2 exc.",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, class2_exc);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t class2_exc_n;
    pointD *class2_exc_out_v[3] = {};
    generate_izk_class2_exc(&class2_exc_n, class2_exc_out_v);

    s = (series_t) { pts: class2_exc_out_v[0], n: class2_exc_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: class2_exc_out_v[1], n: class2_exc_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: class2_exc_out_v[2], n: class2_exc_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (I) spike latency %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t spike_latency =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(I) spike latency",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, spike_latency);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t spike_latency_n;
    pointD *spike_latency_out_v[3] = {};
    generate_izk_spike_latency(&spike_latency_n, spike_latency_out_v);

    s = (series_t) { pts: spike_latency_out_v[0], n: spike_latency_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: spike_latency_out_v[1], n: spike_latency_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: spike_latency_out_v[2], n: spike_latency_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (J) subthresh. osc. %%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t subthr_osc =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(J) subthresh. osc.",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, subthr_osc);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t subthr_osc_n;
    pointD *subthr_osc_out_v[3] = {};
    generate_izk_subthr_osc(&subthr_osc_n, subthr_osc_out_v);

    s = (series_t) { pts: subthr_osc_out_v[0], n: subthr_osc_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: subthr_osc_out_v[1], n: subthr_osc_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: subthr_osc_out_v[2], n: subthr_osc_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx = bwx;
    wy += dwy;

    // %%%%%%%%%%%%%%%%%% (K) resonator %%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t resonator =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(K) resonator",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, resonator);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t resonator_n;
    pointD *resonator_out_v[3] = {};
    generate_izk_resonator(&resonator_n, resonator_out_v);

    s = (series_t) { pts: resonator_out_v[0], n: resonator_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: resonator_out_v[1], n: resonator_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: resonator_out_v[2], n: resonator_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%% (L) integrator %%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t integrator =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(L) integrator",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, integrator);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t integrator_n;
    pointD *integrator_out_v[3] = {};
    generate_izk_integrator(&integrator_n, integrator_out_v);

    s = (series_t) { pts: integrator_out_v[0], n: integrator_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: integrator_out_v[1], n: integrator_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: integrator_out_v[2], n: integrator_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (M) rebound spike %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t rebound_spike =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(M) rebound spike",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, rebound_spike);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t rebound_spike_n;
    pointD *rebound_spike_out_v[3] = {};
    generate_izk_rebound_spike(&rebound_spike_n, rebound_spike_out_v);

    s = (series_t) { pts: rebound_spike_out_v[0], n: rebound_spike_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: rebound_spike_out_v[1], n: rebound_spike_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: rebound_spike_out_v[2], n: rebound_spike_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (N) rebound burst %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t rebound_burst =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(N) rebound burst",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, rebound_burst);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t rebound_burst_n;
    pointD *rebound_burst_out_v[3] = {};
    generate_izk_rebound_burst(&rebound_burst_n, rebound_burst_out_v);

    s = (series_t) { pts: rebound_burst_out_v[0], n: rebound_burst_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: rebound_burst_out_v[1], n: rebound_burst_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: rebound_burst_out_v[2], n: rebound_burst_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%%%%% (O) thresh. variability %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t thresh_variability =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(O) thresh. variability",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, thresh_variability);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t thresh_variability_n;
    pointD *thresh_variability_out_v[3] = {};
    generate_izk_thresh_variability(&thresh_variability_n, thresh_variability_out_v);

    s = (series_t) { pts: thresh_variability_out_v[0], n: thresh_variability_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: thresh_variability_out_v[1], n: thresh_variability_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: thresh_variability_out_v[2], n: thresh_variability_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx = bwx;
    wy += dwy;

    // %%%%%%%%%%%%%% (P) bistability %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t bistability =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(P) bistability",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, bistability);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t bistability_n;
    pointD *bistability_out_v[3] = {};
    generate_izk_bistability(&bistability_n, bistability_out_v);

    s = (series_t) { pts: bistability_out_v[0], n: bistability_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: bistability_out_v[1], n: bistability_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: bistability_out_v[2], n: bistability_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%% (Q) DAP %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t DAP =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(Q) DAP",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, DAP);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t DAP_n;
    pointD *DAP_out_v[3] = {};
    generate_izk_DAP(&DAP_n, DAP_out_v);

    s = (series_t) { pts: DAP_out_v[0], n: DAP_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: DAP_out_v[1], n: DAP_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: DAP_out_v[2], n: DAP_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%% (R) accomodation %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t accomodation =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(R) accomodation",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, accomodation);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t accomodation_n;
    pointD *accomodation_out_v[3] = {};
    generate_izk_accomodation(&accomodation_n, accomodation_out_v);

    s = (series_t) { pts: accomodation_out_v[0], n: accomodation_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: accomodation_out_v[1], n: accomodation_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: accomodation_out_v[2], n: accomodation_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%% (S) inhibition induced spiking %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t inh_induced_sp =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(S) inh. induced sp.",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, inh_induced_sp);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t inh_induced_sp_n;
    pointD *inh_induced_sp_out_v[3] = {};
    generate_izk_inh_induced_sp(&inh_induced_sp_n, inh_induced_sp_out_v);

    s = (series_t) { pts: inh_induced_sp_out_v[0], n: inh_induced_sp_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: inh_induced_sp_out_v[1], n: inh_induced_sp_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: inh_induced_sp_out_v[2], n: inh_induced_sp_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;

    // %%%%%%%%%%%%%% (T) inhibition induced bursting %%%%%%%%%%%%%%%%%%%%%%%%%%
    plot_widget_t inh_induced_brst =
            {
                viewport_px: (Rectangle){ wx, wy, ww, wh },
                xmin:        -7.0,
                xmax:         7.0,
                ymin:        -1.0,
                ymax:         1.0,
                show_grid:    true,
                show_ticks:   true,
                show_labels:  true,
                show_legend:  false,
                dragging:     false,
                title:        "(T) inh. induced brst.",
                series:       NULL,
                nseries:      3
            };
    arrput(widgets, inh_induced_brst);
    widgets_n = arrlen(widgets) - 1;

    // create series
    size_t inh_induced_brst_n;
    pointD *inh_induced_brst_out_v[3] = {};
    generate_izk_inh_induced_brst(&inh_induced_brst_n, inh_induced_brst_out_v);

    s = (series_t) { pts: inh_induced_brst_out_v[0], n: inh_induced_brst_n, color: YELLOW, visible: true, name: "I" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: inh_induced_brst_out_v[1], n: inh_induced_brst_n, color: RED, visible: true, name: "V" };
    arrput(widgets[widgets_n].series, s);
    s = (series_t) { pts: inh_induced_brst_out_v[2], n: inh_induced_brst_n, color: BLUE, visible: true, name: "u" };
    arrput(widgets[widgets_n].series, s);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wx += dwx;
    wy += 0;
#endif

    // main loop
    while (!WindowShouldClose()) {
        BeginDrawing(); {
            ClearBackground(BLACK);

            for (int wi = 0; wi < arrlen(widgets); ++wi) {
                // UI background per viewport
                DrawRectangleRec(widgets[wi].viewport_px, (Color){245,245,245,32});

                // Draw grid and ticks
                draw_grid_and_ticks(&widgets[wi]);

                // Draw all series
                for (size_t i = 0; i < widgets[wi].nseries; ++i) {
                    draw_series(&widgets[wi], &widgets[wi].series[i]);
                }

                // Draw adaptive analytic function clipped to viewport
                int sx = (int)fmax(0, floor(widgets[wi].viewport_px.x));
                int sy = (int)fmax(0, floor(widgets[wi].viewport_px.y));
                int sw = (int)fmax(0, ceil(widgets[wi].viewport_px.width));
                int sh = (int)fmax(0, ceil(widgets[wi].viewport_px.height));
                //            BeginScissorMode(sx, sy, sw, sh); {
                //                draw_adaptive_curve(&widget, widget.xmin, widget.xmax);
                //            } EndScissorMode();

                // Frame and title
                draw_widget_frame(&widgets[wi]);

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
                    BeginScissorMode((int)widgets[wi].viewport_px.x, (int)widgets[wi].viewport_px.y, (int)widgets[wi].viewport_px.width, (int)widgets[wi].viewport_px.height);
                    for (size_t i = 0; i < g_multi.count; ++i) {
                        size_t si = g_multi.series_idx[i];
                        size_t pi = g_multi.pt_idx[i];
                        if (si >= widgets[wi].nseries) continue;
                        pointD p = widgets[wi].series[si].pts[pi];
                        Vector2 sp = world_to_screen(&widgets[wi], p.x, p.y);
                        DrawCircleLinesV(sp, 4, WHITE);
                        DrawCircleLinesV(sp, 2, widgets[wi].series[si].color);
                    }
                    EndScissorMode();
                    // optionally draw a small summary box near top-left of viewport
                    char buf[256];
                    //                int lines = 0;
                    int max_lines = 6;
                    int shown = (int)fmin((double)g_multi.count, (double)max_lines);
                    snprintf(buf, sizeof(buf), "Selected: %zu pts", g_multi.count);
                    DrawText(buf, (int)(widgets[wi].viewport_px.x + 6), (int)(widgets[wi].viewport_px.y + 6), 12, DARKGRAY);
                    for (int i = 0; i < shown; ++i) {
                        size_t si = g_multi.series_idx[i];
                        size_t pi = g_multi.pt_idx[i];
                        pointD p = widgets[wi].series[si].pts[pi];
                        snprintf(buf, sizeof(buf), "%s: x=%.6g y=%.6g", widgets[wi].series[si].name, p.x, p.y);
                        DrawText(buf, (int)(widgets[wi].viewport_px.x + 6), (int)(widgets[wi].viewport_px.y + 24 + i * 14), 12, widgets[wi].series[si].color);
                    }
                    if (g_multi.count > (size_t)shown) {
                        snprintf(buf, sizeof(buf), "... +%zu more", g_multi.count - shown);
                        DrawText(buf, (int)(widgets[wi].viewport_px.x + 6), (int)(widgets[wi].viewport_px.y + 24 + shown * 14), 12, DARKGRAY);
                    }
                }

                // Legend (and get rectangle)
                if (widgets[wi].show_legend) {
                    Rectangle legend_box = draw_legend(&widgets[wi], widgets[wi].series, widgets[wi].nseries);
                    // handle clicks on legend
                    handle_legend_input(&widgets[wi], widgets[wi].series, widgets[wi].nseries, legend_box);
                }

                // input handling (pan/zoom)
                handle_input(&widgets[wi]);

                Vector2 mpos = GetMousePosition();
                const double pick_threshold_px = 8.0;
                size_t found_si = 0, found_pi = 0;
                pointD found_pt;
                double found_d;
                bool hit = find_nearest_point(&widgets[wi], widgets[wi].series, widgets[wi].nseries, mpos, pick_threshold_px, &found_si, &found_pi, &found_pt, &found_d);

                // tooltip (transient) si disegna poi con draw_mouse_tooltip
                // selection toggle on left click
                if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) {
                    // se click dentro legend l'handler legend lo gestisce prima; quindi qui consideriamo solo viewport
                    if (mpos.x >= widgets[wi].viewport_px.x && mpos.x <= widgets[wi].viewport_px.x + widgets[wi].viewport_px.width &&
                            mpos.y >= widgets[wi].viewport_px.y && mpos.y <= widgets[wi].viewport_px.y + widgets[wi].viewport_px.height) {
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
                        if (g_selection.series_idx < widgets[wi].nseries) {
                            series_t *s = &widgets[wi].series[g_selection.series_idx];
                            Vector2 sp = world_to_screen(&widgets[wi], g_selection.pt.x, g_selection.pt.y);
                            // cerchio di evidenziazione (fuori scissor per visibilità) o dentro viewport
                            DrawCircleV(sp, 4, WHITE);
                            DrawCircleV(sp, 3, s->color);
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
                    draw_mouse_tooltip(&widgets[wi], widgets[wi].series, widgets[wi].nseries, mpos);
                } EndScissorMode();

                if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
                    char lab[64] = {};
                    if (mod_shift() && !mod_ctrl()) snprintf(lab, sizeof(lab), "Pan: X only");
                    else if (mod_ctrl() && !mod_shift()) snprintf(lab, sizeof(lab), "Pan: Y only");
                    else snprintf(lab, sizeof(lab), "Pan: XY");
                    DrawText(lab, (int)widgets[wi].viewport_px.x + 8, (int)widgets[wi].viewport_px.y + 8, 10, DARKGRAY);
                }

                // instructions
//                DrawText("Pan: Middle drag (Shift: X only, Ctrl: Y only)  Zoom: Wheel (Shift: X only, Ctrl: Y only)", 8, 8, 12, DARKGRAY);
            }
        } EndDrawing();
    }

    // cleanup
    for (int oi = 0; oi < 3; oi++) {
        //~ free(tonic_spike_out_v[oi]);
        //~ free(phasic_spike_out_v[oi]);
        //~ free(tonic_burst_out_v[oi]);
        //~ free(phasic_burst_out_v[oi]);
        //~ free(mixed_mode_out_v[oi]);
        //~ free(spike_freq_adapt_out_v[oi]);
        //~ free(class1_exc_out_v[oi]);
        //~ free(class2_exc_out_v[oi]);
        //~ free(spike_latency_out_v[oi]);
    }
    
    for (int wi = 0; wi < arrlen(widgets); ++wi) {
        arrfree(widgets[wi].series);
    }
    arrfree(widgets);

    CloseWindow();
    return 0;
}
