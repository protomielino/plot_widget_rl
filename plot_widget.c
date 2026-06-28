#define RAYMATH_IMPLEMENTATION
#include "plot_widget.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ------------------------------------------------------------------ */
/*  coordinate transforms                                             */
/* ------------------------------------------------------------------ */

Vector2 world_to_screen(const plot_widget_t *w, double wx, double wy)
{
    double sx = (wx - w->xmin) / (w->xmax - w->xmin) * w->viewport_px.width + w->viewport_px.x;
    double sy = (1.0 - (wy - w->ymin) / (w->ymax - w->ymin)) * w->viewport_px.height + w->viewport_px.y;
    return (Vector2){ (float)sx, (float)sy };
}

void screen_to_world(const plot_widget_t *w, double sx, double sy,
                     double *wx, double *wy)
{
    double nx = (sx - w->viewport_px.x) / w->viewport_px.width;
    double ny = (sy - w->viewport_px.y) / w->viewport_px.height;
    *wx = lerp(w->xmin, w->xmax, nx);
    *wy = lerp(w->ymax, w->ymin, ny);
}

void rect_screen_to_world_box(const plot_widget_t *w, Rectangle r,
                              double *oxmin, double *oxmax,
                              double *oymin, double *oymax)
{
    double x0, y0, x1, y1;
    screen_to_world(w, r.x, r.y, &x0, &y0);
    screen_to_world(w, r.x + r.width, r.y + r.height, &x1, &y1);
    *oxmin = fmin(x0, x1);
    *oxmax = fmax(x0, x1);
    *oymin = fmin(y0, y1);
    *oymax = fmax(y0, y1);
}

/* ------------------------------------------------------------------ */
/*  internal helpers                                                  */
/* ------------------------------------------------------------------ */

static double mouse_point_distance_px(const plot_widget_t *w,
                                      const pointD *p, Vector2 mouse)
{
    Vector2 sp = world_to_screen(w, p->x, p->y);
    double dx = sp.x - mouse.x, dy = sp.y - mouse.y;
    return sqrt(dx * dx + dy * dy);
}



/* ------------------------------------------------------------------ */
/*  public helpers                                                     */
/* ------------------------------------------------------------------ */

double nice_tick_step(double range, int target_ticks)
{
    if (range <= 0) return 0.0;
    double raw = range / target_ticks;
    double mag = pow(10.0, floor(log10(raw)));
    double norm = raw / mag;
    double nice;
    if (norm < 1.5)      nice = 1.0;
    else if (norm < 3.0) nice = 2.0;
    else if (norm < 7.0) nice = 5.0;
    else                 nice = 10.0;
    return nice * mag;
}

void get_series_ranges(const series_t *s, int ns,
                       double *xmin, double *xmax,
                       double *ymin, double *ymax)
{
    *xmin = DBL_MAX; *xmax = -DBL_MAX;
    *ymin = DBL_MAX; *ymax = -DBL_MAX;
    for (int si = 0; si < ns; ++si) {
        const pointD *pts = s[si].pts;
        size_t n = s[si].n;
        for (size_t pi = 0; pi < n; ++pi) {
            if (!isfinite(pts[pi].x) || !isfinite(pts[pi].y)) continue;
            *xmin = fmin(*xmin, pts[pi].x);
            *xmax = fmax(*xmax, pts[pi].x);
            *ymin = fmin(*ymin, pts[pi].y);
            *ymax = fmax(*ymax, pts[pi].y);
        }
    }
}

bool find_nearest_point(const plot_widget_t *w, const series_t *series,
                        size_t nseries, Vector2 mouse, double thresh_px,
                        size_t *out_si, size_t *out_pi,
                        pointD *out_pt, double *out_dist)
{
    double bestd = 1e9;
    size_t best_si = 0, best_pi = 0;
    pointD best_pt = {0, 0};

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
    if (bestd <= thresh_px) {
        *out_si = best_si; *out_pi = best_pi;
        *out_pt = best_pt; *out_dist = bestd;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  drawing                                                            */
/* ------------------------------------------------------------------ */

void draw_grid_and_ticks(const plot_widget_t *w)
{
    /* border */
    DrawRectangleLinesEx(w->viewport_px, 2, GRAY);

    double xrange = w->xmax - w->xmin;
    double yrange = w->ymax - w->ymin;
    if (xrange <= 0 || yrange <= 0) return;

    /* nice tick steps */
    int target_xticks = (int)fmax(4, w->viewport_px.width / 120.0);
    int target_yticks = (int)fmax(4, w->viewport_px.height / 80.0);
    double xstep = nice_tick_step(xrange, target_xticks);
    double ystep = nice_tick_step(yrange, target_yticks);

    double startx = ceil(w->xmin / xstep - 1e-12) * xstep;
    double starty = ceil(w->ymin / ystep - 1e-12) * ystep;

    Color grid_col  = (Color){200, 200, 200, 80};
    Color x_axis_col = (Color){180, 50, 50, 200};
    Color y_axis_col = (Color){50, 180, 50, 200};
    Color tick_col  = (Color){200, 200, 200, 255};
    Color label_col = (Color){180, 180, 180, 255};

    if (w->show_grid) {
        int sx = (int)fmax(0, floor(w->viewport_px.x));
        int sy = (int)fmax(0, floor(w->viewport_px.y));
        int sw = (int)fmax(0, ceil(w->viewport_px.width));
        int sh = (int)fmax(0, ceil(w->viewport_px.height));
        BeginScissorMode(sx, sy, sw, sh); {
            /* vertical grid lines (x) */
            for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
                Vector2 p1 = world_to_screen(w, tx, w->ymin);
                Vector2 p2 = world_to_screen(w, tx, w->ymax);
                DrawLineEx(p1, p2, 1.0f, grid_col);
            }
            /* horizontal grid lines (y) */
            for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
                Vector2 p1 = world_to_screen(w, w->xmin, ty);
                Vector2 p2 = world_to_screen(w, w->xmax, ty);
                DrawLineEx(p1, p2, 1.0f, grid_col);
            }
            /* axes x=0 and y=0 */
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
            /* world range label (lower right) */
            char buf[128];
            snprintf(buf, sizeof(buf), "x:[%.3g, %.3g] y:[%.3g, %.3g]",
                     w->xmin, w->xmax, w->ymin, w->ymax);
            int tw = MeasureText(buf, 10);
            DrawText(buf,
                     (int)(w->viewport_px.x + w->viewport_px.width - tw - 6),
                     (int)(w->viewport_px.y + w->viewport_px.height - 26),
                     10, DARKGRAY);
        } EndScissorMode();
    }

    if (w->show_ticks) {
        int tick_len = 6;
        for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
            Vector2 p1 = world_to_screen(w, tx, w->ymin);
            DrawLine((int)p1.x, (int)p1.y, (int)p1.x, (int)p1.y + tick_len, tick_col);
        }
        for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
            Vector2 p1 = world_to_screen(w, w->xmin, ty);
            DrawLine((int)p1.x, (int)p1.y, (int)p1.x - tick_len, (int)p1.y, tick_col);
        }
    }

    if (w->show_labels) {
        /* x-axis labels */
        for (double tx = startx; tx <= w->xmax + 1e-12; tx += xstep) {
            Vector2 s0 = world_to_screen(w, tx, w->ymin);
            char buf[64];
            snprintf(buf, sizeof(buf), "%.4g", tx);
            float lx = s0.x - MeasureText(buf, 10) / 2.0f;
            DrawText(buf, (int)lx,
                     (int)(w->viewport_px.y + w->viewport_px.height) + 3,
                     10, label_col);
        }
        /* y-axis labels */
        for (double ty = starty; ty <= w->ymax + 1e-12; ty += ystep) {
            Vector2 s0 = world_to_screen(w, w->xmin, ty);
            char buf[64];
            snprintf(buf, sizeof(buf), "%.4g", ty);
            float lx = w->viewport_px.x - MeasureText(buf, 10) - 3;
            DrawText(buf, (int)lx, (int)s0.y - 8, 10, label_col);
        }
    }
}

void draw_curve(const plot_widget_t *w, const pointD *pts,
                size_t n, Color color)
{
    if (!pts || n < 2) return;

    int sx = (int)fmax(0, floor(w->viewport_px.x));
    int sy = (int)fmax(0, floor(w->viewport_px.y));
    int sw = (int)fmax(0, ceil(w->viewport_px.width));
    int sh = (int)fmax(0, ceil(w->viewport_px.height));
    BeginScissorMode(sx, sy, sw, sh); {
        for (size_t i = 0; i + 1 < n; ++i) {
            pointD a = pts[i], b = pts[i + 1];
            if (!isfinite(a.x) || !isfinite(a.y) ||
                !isfinite(b.x) || !isfinite(b.y))
                continue;
            if ((a.x < w->xmin && b.x < w->xmin) ||
                (a.x > w->xmax && b.x > w->xmax) ||
                (a.y < w->ymin && b.y < w->ymin) ||
                (a.y > w->ymax && b.y > w->ymax))
                continue;
            Vector2 sa = world_to_screen(w, a.x, a.y);
            Vector2 sb = world_to_screen(w, b.x, b.y);
            DrawLineV(sa, sb, color);
        }
    } EndScissorMode();
}

void draw_series(const plot_widget_t *w, const series_t *s)
{
    if (!s->visible) return;
    draw_curve(w, s->pts, s->n, s->color);
}

void draw_widget_frame(const plot_widget_t *w)
{
    Rectangle tb = title_bar_rect(w);
    DrawRectangleRec(tb, (Color){245, 245, 245, 16});
    DrawLineEx(
        (Vector2){ tb.x, tb.y + tb.height },
        (Vector2){ tb.x + tb.width, tb.y + tb.height },
        1.0f, (Color){245, 245, 245, 40});
    int tw = MeasureText(w->title, 12);
    DrawText(w->title,
             (int)(w->viewport_px.x + (w->viewport_px.width - tw) / 2.0),
             (int)(w->viewport_px.y - 16), 12, LIGHTGRAY);
}

void draw_resize_handle(const plot_widget_t *w)
{
    Rectangle hr = resize_handle_rect(w);
    float hx = hr.x, hy = hr.y, hs = hr.width;
    Vector2 p1 = { hx + hs, hy };
    Vector2 p2 = { hx + hs, hy + hs };
    Vector2 p3 = { hx,      hy + hs };
    DrawTriangle(p1, p2, p3, (Color){200, 200, 200, 90});
    DrawLineEx((Vector2){ hx + hs - 5, hy + 1 }, (Vector2){ hx + hs - 1, hy + 5 }, 1.0f, (Color){180, 180, 180, 120});
    DrawLineEx((Vector2){ hx + hs - 8, hy + 1 }, (Vector2){ hx + hs - 1, hy + 8 }, 1.0f, (Color){180, 180, 180, 120});
}

void draw_mouse_tooltip(const plot_widget_t *w, const series_t *series,
                        size_t nseries, Vector2 mouse)
{
    if (mouse.x < w->viewport_px.x ||
        mouse.x > w->viewport_px.x + w->viewport_px.width ||
        mouse.y < w->viewport_px.y ||
        mouse.y > w->viewport_px.y + w->viewport_px.height)
        return;

    double mx, my;
    screen_to_world(w, mouse.x, mouse.y, &mx, &my);

    char buf[256];
    snprintf(buf, sizeof(buf), "x=%.6g  y=%.6g", mx, my);
    DrawText(buf, (int)mouse.x + 12, (int)mouse.y - 6, 12, DARKGRAY);

    float ty = mouse.y + 16;
    for (size_t si = 0; si < nseries; ++si) {
        if (!series[si].visible) continue;
        const series_t *s = &series[si];
        size_t best_i = SIZE_MAX;
        double best_dx = 1e9;
        for (size_t i = 0; i < s->n; ++i) {
            if (!isfinite(s->pts[i].x) || !isfinite(s->pts[i].y)) continue;
            double dx = fabs(s->pts[i].x - mx);
            if (dx < best_dx) { best_dx = dx; best_i = i; }
        }
        if (best_i != SIZE_MAX) {
            snprintf(buf, sizeof(buf), "%s: x=%.6g y=%.6g",
                     s->name, s->pts[best_i].x, s->pts[best_i].y);
            DrawText(buf, (int)mouse.x + 12, (int)ty, 12, s->color);
            ty += 14;
        }
    }
}

Rectangle draw_legend(const plot_widget_t *w, const series_t *series,
                      size_t nseries)
{
    float pad = 8, entry_h = 20;
    float box_w = 200;
    float box_h = pad * 2 + nseries * entry_h;
    float x = w->viewport_px.x + w->viewport_px.width - box_w - 8;
    float y = w->viewport_px.y + 8;

    DrawRectangle((int)x, (int)y, (int)box_w, (int)box_h,
                  (Color){250, 250, 250, 16});
    DrawRectangleLines((int)x, (int)y, (int)box_w, (int)box_h, GRAY);

    for (size_t i = 0; i < nseries; ++i) {
        float ey = y + pad + i * entry_h;
        DrawRectangle((int)(x + 6), (int)(ey + 4), 12, 12, series[i].color);
        DrawText(series[i].name, (int)(x + 26), (int)(ey), 12,
                 series[i].visible ? LIGHTGRAY : GRAY);
    }
    return (Rectangle){ x, y, box_w, box_h };
}

/* ------------------------------------------------------------------ */
/*  input handling                                                     */
/* ------------------------------------------------------------------ */

void handle_input(plot_widget_t *w, int w_idx, bool suppress_left_click)
{
    Vector2 mpos = GetMousePosition();
    bool inside = CheckCollisionPointRec(mpos, w->viewport_px);

    /* -------- pan (middle mouse) -------- */
    if (IsMouseButtonPressed(MOUSE_MIDDLE_BUTTON) && inside) {
        w->dragging = true;
        w->last_mouse = mpos;
    }
    if (IsMouseButtonReleased(MOUSE_MIDDLE_BUTTON))
        w->dragging = false;

    if (w->dragging) {
        Vector2 d = Vector2Subtract(mpos, w->last_mouse);
        w->last_mouse = mpos;

        double dx_full = -d.x / w->viewport_px.width  * (w->xmax - w->xmin);
        double dy_full =  d.y / w->viewport_px.height * (w->ymax - w->ymin);

        bool shift = mod_shift(), ctrl = mod_ctrl();
        if (shift && !ctrl) {
            w->xmin += dx_full; w->xmax += dx_full;
        } else if (ctrl && !shift) {
            w->ymin += dy_full; w->ymax += dy_full;
        } else {
            w->xmin += dx_full; w->xmax += dx_full;
            w->ymin += dy_full; w->ymax += dy_full;
        }
    }

    /* -------- zoom (mouse wheel) -------- */
    float wheel = GetMouseWheelMove();
    if (wheel != 0.0f && inside) {
        double mx, my;
        screen_to_world(w, mpos.x, mpos.y, &mx, &my);
        double k = pow(1.15, -wheel);

        bool shift = mod_shift(), ctrl = mod_ctrl();
        if (shift && !ctrl) {
            double nxmin = mx + (w->xmin - mx) * k;
            double nxmax = mx + (w->xmax - mx) * k;
            if (nxmax - nxmin > 1e-12) { w->xmin = nxmin; w->xmax = nxmax; }
        } else if (ctrl && !shift) {
            double nymin = my + (w->ymin - my) * k;
            double nymax = my + (w->ymax - my) * k;
            if (nymax - nymin > 1e-12) { w->ymin = nymin; w->ymax = nymax; }
        } else {
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

    /* -------- resize (right-mouse drag on bottom-right handle) -------- */
    if (IsMouseButtonPressed(MOUSE_BUTTON_RIGHT)) {
        if (CheckCollisionPointRec(mpos, resize_handle_rect(w))) {
            w->resizing = true;
            w->last_mouse = mpos;
        }
    }
    if (IsMouseButtonReleased(MOUSE_BUTTON_RIGHT))
        w->resizing = false;
    if (w->resizing) {
        Vector2 d = Vector2Subtract(mpos, w->last_mouse);
        w->last_mouse = mpos;
        float nw = w->viewport_px.width  + d.x;
        float nh = w->viewport_px.height + d.y;
        if (nw >= 80) w->viewport_px.width  = nw;
        if (nh >= 60) w->viewport_px.height = nh;
    }

    /* -------- area selection (left-mouse drag) -------- */
    if (!suppress_left_click && IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && inside) {
        w->area.active = true;
        w->area.finished = false;
        w->area.widget_idx = w_idx;
        w->area.start_px = mpos;
        w->area.end_px = mpos;
        w->multi.has_any = false;
        w->multi.count = 0;
        w->multi.widget_idx = w_idx;
        w->selection.active = false;
    }

    if (w->area.active) {
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
            w->area.end_px = mpos;

        if (IsMouseButtonReleased(MOUSE_LEFT_BUTTON) &&
            w->area.widget_idx == (size_t)w_idx) {
            w->area.active = false;
            if (inside) {
                w->area.finished = true;
                w->area.end_px = mpos;
                w->multi.widget_idx = w_idx;

                free(w->multi.series_idx); w->multi.series_idx = NULL;
                free(w->multi.pt_idx);     w->multi.pt_idx = NULL;
                w->multi.count = 0;
                w->multi.has_any = false;

                Rectangle srect = norm_rect_from_points(w->area.start_px,
                                                        w->area.end_px);
                double wxmin, wxmax, wymin, wymax;
                rect_screen_to_world_box(w, srect,
                                         &wxmin, &wxmax, &wymin, &wymax);

                /* count points inside bounding box */
                for (size_t si = 0; si < w->nseries; ++si) {
                    if (!w->series[si].visible) continue;
                    for (size_t pi = 0; pi < w->series[si].n; ++pi) {
                        pointD p = w->series[si].pts[pi];
                        if (!isfinite(p.x) || !isfinite(p.y)) continue;
                        if (p.x >= wxmin && p.x <= wxmax &&
                            p.y >= wymin && p.y <= wymax)
                            w->multi.count++;
                    }
                }

                if (w->multi.count > 0) {
                    w->multi.series_idx = malloc(sizeof(size_t) * w->multi.count);
                    w->multi.pt_idx     = malloc(sizeof(size_t) * w->multi.count);
                    size_t idx = 0;
                    for (size_t si = 0; si < w->nseries; ++si) {
                        if (!w->series[si].visible) continue;
                        for (size_t pi = 0; pi < w->series[si].n; ++pi) {
                            pointD p = w->series[si].pts[pi];
                            if (!isfinite(p.x) || !isfinite(p.y)) continue;
                            if (p.x >= wxmin && p.x <= wxmax &&
                                p.y >= wymin && p.y <= wymax) {
                                w->multi.series_idx[idx] = si;
                                w->multi.pt_idx[idx] = pi;
                                idx++;
                            }
                        }
                    }
                    w->multi.has_any = true;
                } else {
                    w->multi.has_any = false;
                }
            } else {
                w->area.finished = false;
            }
        }
    }

    /* Cancel area selection with Q */
    if (IsKeyPressed(KEY_Q)) {
        w->area.active = false;
        w->area.finished = false;
    }

    /* Reset view with R (on hovered widget only – caller handles w_idx) */
    if (IsKeyPressed(KEY_R) && inside) {
        w->xmin = -10.0; w->xmax = 10.0;
        w->ymin =  -3.0; w->ymax =  3.0;
    }
}

void handle_legend_input(plot_widget_t *w, series_t *series,
                         size_t nseries, Rectangle legend_box)
{
    (void)w;
    if (!IsMouseButtonPressed(MOUSE_LEFT_BUTTON)) return;
    Vector2 m = GetMousePosition();
    if (m.x < legend_box.x || m.x > legend_box.x + legend_box.width ||
        m.y < legend_box.y || m.y > legend_box.y + legend_box.height)
        return;

    float pad = 8, entry_h = 20;
    for (size_t i = 0; i < nseries; ++i) {
        float ey = legend_box.y + pad + i * entry_h;
        Rectangle entry_r = { legend_box.x, ey, legend_box.width, entry_h };
        if (CheckCollisionPointRec(m, entry_r)) {
            series[i].visible = !series[i].visible;
            break;
        }
    }
}
