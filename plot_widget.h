#ifndef PLOT_WIDGET_H
#define PLOT_WIDGET_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <raylib.h>
#include <raymath.h>

typedef struct {
    double x, y;
} pointD;

typedef struct {
    pointD *pts;
    size_t n;
    Color color;
    char name[64];
    bool visible;
} series_t;

typedef struct {
    bool active;
    size_t widget_idx;
    size_t series_idx;
    size_t pt_idx;
    pointD pt;
} selection_t;

typedef struct {
    bool active;
    bool finished;
    size_t widget_idx;
    Vector2 start_px;
    Vector2 end_px;
} area_sel_t;

typedef struct {
    bool has_any;
    size_t widget_idx;
    size_t *series_idx;
    size_t *pt_idx;
    size_t count;
} multi_selection_t;

typedef struct plot_widget_s {
    Rectangle viewport_px;
    double xmin, xmax;
    double ymin, ymax;
    bool show_grid;
    bool show_ticks;
    bool show_labels;
    bool show_legend;
    bool dragging;
    Vector2 last_mouse;
    char title[128];
    size_t nseries;
    series_t *series;

    /* per-widget interactive state */
    selection_t selection;
    area_sel_t area;
    multi_selection_t multi;

    /* fullscreen toggle support: saved viewport before expanding to full window */
    Rectangle saved_viewport;
} plot_widget_t;

/* utility inlines */
static inline double lerp(double a, double b, double t) { return a + (b - a) * t; }
static inline double clamp_double(double v, double a, double b) { return v < a ? a : (v > b ? b : v); }
static inline bool mod_shift(void) { return IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT); }
static inline bool mod_ctrl(void)  { return IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_RIGHT_CONTROL); }

static inline Rectangle norm_rect_from_points(Vector2 a, Vector2 b)
{
    float x = fminf(a.x, b.x), y = fminf(a.y, b.y);
    return (Rectangle){ x, y, fabsf(a.x - b.x), fabsf(a.y - b.y) };
}

/* coordinate transforms */
Vector2 world_to_screen(const plot_widget_t *w, double wx, double wy);
void   screen_to_world(const plot_widget_t *w, double sx, double sy, double *wx, double *wy);
void   rect_screen_to_world_box(const plot_widget_t *w, Rectangle r, double *oxmin, double *oxmax, double *oymin, double *oymax);

/* range / nearest-point */
void  get_series_ranges(const series_t *s, int ns, double *xmin, double *xmax, double *ymin, double *ymax);
bool  find_nearest_point(const plot_widget_t *w, const series_t *series, size_t nseries, Vector2 mouse, double thresh_px, size_t *out_si, size_t *out_pi, pointD *out_pt, double *out_dist);

/* drawing */
void draw_grid_and_ticks(const plot_widget_t *w);
void draw_curve(const plot_widget_t *w, const pointD *pts, size_t n, Color color);
void draw_series(const plot_widget_t *w, const series_t *s);
void draw_widget_frame(const plot_widget_t *w);
void draw_mouse_tooltip(const plot_widget_t *w, const series_t *series, size_t nseries, Vector2 mouse);
Rectangle draw_legend(const plot_widget_t *w, const series_t *series, size_t nseries);

/* input */
void handle_input(plot_widget_t *w, int w_idx);
void handle_legend_input(plot_widget_t *w, series_t *series, size_t nseries, Rectangle legend_box);

/* helpers */
double nice_tick_step(double range, int target_ticks);

#endif /* PLOT_WIDGET_H */
