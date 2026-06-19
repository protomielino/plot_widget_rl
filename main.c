#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "plot_widget.h"
#include "izhikevich.h"

/* ------------------------------------------------------------------ */
/*  simple dynamic array for widgets                                   */
/* ------------------------------------------------------------------ */
typedef struct {
    plot_widget_t *items;
    size_t count;
    size_t capacity;
} widget_arr_t;

static void arr_append(widget_arr_t *a, plot_widget_t item)
{
    if (a->count >= a->capacity) {
        a->capacity = a->capacity ? a->capacity * 2 : 32;
        a->items = realloc(a->items, a->capacity * sizeof(plot_widget_t));
    }
    a->items[a->count++] = item;
}

/* ------------------------------------------------------------------ */
/*  auto-fit viewport to data                                          */
/* ------------------------------------------------------------------ */
static void auto_fit_view(plot_widget_t *w)
{
    double xmin, xmax, ymin, ymax;
    get_series_ranges(w->series, (int)w->nseries, &xmin, &xmax, &ymin, &ymax);
    double xpad = (xmax - xmin) * 0.05;
    double ypad = (ymax - ymin) * 0.05;
    if (xpad < 1e-12) xpad = 1.0;
    if (ypad < 1e-12) ypad = 1.0;
    w->xmin = xmin - xpad;
    w->xmax = xmax + xpad;
    w->ymin = ymin - ypad;
    w->ymax = ymax + ypad;
}

/* ------------------------------------------------------------------ */
/*  allocate 3 series (I, V, u) in a widget                            */
/* ------------------------------------------------------------------ */
static void assign_3_series(plot_widget_t *w, size_t n,
                            pointD *i_arr, pointD *v_arr, pointD *u_arr)
{
    w->series = malloc(3 * sizeof(series_t));
    w->nseries = 3;
    w->series[0] = (series_t){ .pts = i_arr, .n = n, .color = YELLOW, .visible = true, .name = "I" };
    w->series[1] = (series_t){ .pts = v_arr, .n = n, .color = RED,    .visible = true, .name = "V" };
    w->series[2] = (series_t){ .pts = u_arr, .n = n, .color = BLUE,   .visible = true, .name = "u" };
}

/* ------------------------------------------------------------------ */
/*  helper: create one widget + generate data + auto-fit               */
/* ------------------------------------------------------------------ */
static void add_widget(widget_arr_t *da, Rectangle vp, const char *title,
                       pointD *(*gen)(size_t *, pointD **))
{
    plot_widget_t w = {
        .viewport_px  = vp,
        .xmin = -10.0, .xmax = 10.0,
        .ymin = -3.0,  .ymax = 3.0,
        .show_grid   = true,
        .show_ticks  = true,
        .show_labels = true,
        .show_legend = false,
        .dragging    = false,
        .title       = "",
        .series      = NULL,
        .nseries     = 0,
        .selection   = {0},
        .area        = {0},
        .multi       = {0},
    };
    snprintf(w.title, sizeof(w.title), "%s", title);

    size_t n;
    pointD *out[3] = { NULL, NULL, NULL };
    gen(&n, out);

    assign_3_series(&w, n, out[0], out[1], out[2]);
    auto_fit_view(&w);

    arr_append(da, w);
}

int main(void)
{
    srand((unsigned int)clock());

    const int screenW = 1150;
    const int screenH = 750;
    InitWindow(screenW, screenH, "Plot Widget - raylib");
    SetTargetFPS(60);

    /* layout constants */
    const float ww = 200, wh = 150;
    const int bwx = 25, bwy = 25;
    const int wsx = 25, wsy = 35;
    const int dwx = ww + wsx, dwy = wh + wsy;

    widget_arr_t widgets = {0};

    /* ---------- row 1 ---------- */
    add_widget(&widgets, (Rectangle){ bwx + 0*dwx, bwy + 0*dwy, ww, wh }, "(A) tonic spiking",              generate_izk_tonic_spike);
    add_widget(&widgets, (Rectangle){ bwx + 1*dwx, bwy + 0*dwy, ww, wh }, "(B) phasic spiking",             generate_izk_phasic_spike);
    add_widget(&widgets, (Rectangle){ bwx + 2*dwx, bwy + 0*dwy, ww, wh }, "(C) tonic bursting",            generate_izk_tonic_bursting);
    add_widget(&widgets, (Rectangle){ bwx + 3*dwx, bwy + 0*dwy, ww, wh }, "(D) phasic bursting",           generate_izk_phasic_bursting);
    add_widget(&widgets, (Rectangle){ bwx + 4*dwx, bwy + 0*dwy, ww, wh }, "(E) mixed mode",                generate_izk_mixed_mode);

    /* ---------- row 2 ---------- */
    add_widget(&widgets, (Rectangle){ bwx + 0*dwx, bwy + 1*dwy, ww, wh }, "(F) spike freq. adapt",         generate_izk_spike_freq_adapt);
    add_widget(&widgets, (Rectangle){ bwx + 1*dwx, bwy + 1*dwy, ww, wh }, "(G) Class 1 exc.",              generate_izk_class1_exc);
    add_widget(&widgets, (Rectangle){ bwx + 2*dwx, bwy + 1*dwy, ww, wh }, "(H) Class 2 exc.",              generate_izk_class2_exc);
    add_widget(&widgets, (Rectangle){ bwx + 3*dwx, bwy + 1*dwy, ww, wh }, "(I) spike latency",             generate_izk_spike_latency);
    add_widget(&widgets, (Rectangle){ bwx + 4*dwx, bwy + 1*dwy, ww, wh }, "(J) subthresh. osc.",           generate_izk_subthr_osc);

    /* ---------- row 3 ---------- */
    add_widget(&widgets, (Rectangle){ bwx + 0*dwx, bwy + 2*dwy, ww, wh }, "(K) resonator",                 generate_izk_resonator);
    add_widget(&widgets, (Rectangle){ bwx + 1*dwx, bwy + 2*dwy, ww, wh }, "(L) integrator",                generate_izk_integrator);
    add_widget(&widgets, (Rectangle){ bwx + 2*dwx, bwy + 2*dwy, ww, wh }, "(M) rebound spike",             generate_izk_rebound_spike);
    add_widget(&widgets, (Rectangle){ bwx + 3*dwx, bwy + 2*dwy, ww, wh }, "(N) rebound burst",             generate_izk_rebound_burst);
    add_widget(&widgets, (Rectangle){ bwx + 4*dwx, bwy + 2*dwy, ww, wh }, "(O) thresh. variability",       generate_izk_thresh_variability);

    /* ---------- row 4 ---------- */
    add_widget(&widgets, (Rectangle){ bwx + 0*dwx, bwy + 3*dwy, ww, wh }, "(P) bistability",               generate_izk_bistability);
    add_widget(&widgets, (Rectangle){ bwx + 1*dwx, bwy + 3*dwy, ww, wh }, "(Q) DAP",                       generate_izk_DAP);
    add_widget(&widgets, (Rectangle){ bwx + 2*dwx, bwy + 3*dwy, ww, wh }, "(R) accomodation",              generate_izk_accomodation);
    add_widget(&widgets, (Rectangle){ bwx + 3*dwx, bwy + 3*dwy, ww, wh }, "(S) inh. induced sp.",          generate_izk_inh_induced_sp);
    add_widget(&widgets, (Rectangle){ bwx + 4*dwx, bwy + 3*dwy, ww, wh }, "(T) inh. induced brst.",        generate_izk_inh_induced_brst);

    /* ---------------------------------------------------------------- */
    /*  main loop                                                       */
    /* ---------------------------------------------------------------- */
    int fullscreen_idx = -1;

    while (!WindowShouldClose()) {
        Vector2 gpos = GetMousePosition();
        int active_wi = -1;
        for (int wi = 0; wi < (int)widgets.count; ++wi) {
            if (CheckCollisionPointRec(gpos, widgets.items[wi].viewport_px)) {
                active_wi = wi;
                break;
            }
        }

        BeginDrawing(); {
            ClearBackground(BLACK);

            for (int wi = 0; wi < (int)widgets.count; ++wi) {
                if (fullscreen_idx >= 0 && wi != fullscreen_idx) continue;
                plot_widget_t *w = &widgets.items[wi];

                /* background */
                DrawRectangleRec(w->viewport_px, (Color){245, 245, 245, 32});

                /* grid / ticks / labels */
                draw_grid_and_ticks(w);

                /* series */
                for (size_t i = 0; i < w->nseries; ++i)
                    draw_series(w, &w->series[i]);

                /* frame / title */
                draw_widget_frame(w);

                /* highlight active widget */
                if (wi == active_wi)
                    DrawRectangleLinesEx(w->viewport_px, 2,
                                         (Color){255, 255, 100, 160});

                /* area drag rectangle (in progress) */
                if (w->area.active && w->area.widget_idx == (size_t)wi) {
                    Rectangle r = norm_rect_from_points(w->area.start_px,
                                                        w->area.end_px);
                    DrawRectangleRec(r, (Color){200, 200, 255, 60});
                    DrawRectangleLines((int)r.x, (int)r.y,
                                       (int)r.width, (int)r.height,
                                       (Color){100, 100, 255, 180});
                }

                /* highlight multi-selected points */
                if (w->area.finished && w->multi.has_any &&
                    w->multi.widget_idx == (size_t)wi) {
                    int sx = (int)w->viewport_px.x;
                    int sy = (int)w->viewport_px.y;
                    int sw = (int)w->viewport_px.width;
                    int sh = (int)w->viewport_px.height;
                    BeginScissorMode(sx, sy, sw, sh); {
                        for (size_t i = 0; i < w->multi.count; ++i) {
                            size_t si = w->multi.series_idx[i];
                            size_t pi = w->multi.pt_idx[i];
                            if (si >= w->nseries) continue;
                            pointD p = w->series[si].pts[pi];
                            Vector2 sp = world_to_screen(w, p.x, p.y);
                            DrawCircleLinesV(sp, 4, WHITE);
                            DrawCircleLinesV(sp, 2, w->series[si].color);
                        }
                    } EndScissorMode();

                    /* summary box */
                    char buf[256];
                    int max_lines = 6;
                    int shown = (int)fmin((double)w->multi.count,
                                          (double)max_lines);
                    snprintf(buf, sizeof(buf), "Selected: %zu pts",
                             w->multi.count);
                    DrawText(buf, (int)w->viewport_px.x + 6,
                             (int)w->viewport_px.y + 6, 12, DARKGRAY);
                    for (int i = 0; i < shown; ++i) {
                        size_t si = w->multi.series_idx[i];
                        size_t pi = w->multi.pt_idx[i];
                        pointD p = w->series[si].pts[pi];
                        snprintf(buf, sizeof(buf), "%s: x=%.6g y=%.6g",
                                 w->series[si].name, p.x, p.y);
                        DrawText(buf, (int)w->viewport_px.x + 6,
                                 (int)w->viewport_px.y + 24 + i * 14,
                                 12, w->series[si].color);
                    }
                    if (w->multi.count > (size_t)shown) {
                        snprintf(buf, sizeof(buf), "... +%zu more",
                                 w->multi.count - shown);
                        DrawText(buf, (int)w->viewport_px.x + 6,
                                 (int)w->viewport_px.y + 24 + shown * 14,
                                 12, DARKGRAY);
                    }
                }

                /* legend */
                if (w->show_legend) {
                    Rectangle lb = draw_legend(w, w->series, w->nseries);
                    handle_legend_input(w, w->series, w->nseries, lb);
                }

                /* pan/zoom/selection input */
                handle_input(w, wi);

                /* single-point selection tooltip */
                Vector2 mpos = GetMousePosition();
                const double pick_thresh = 8.0;
                size_t found_si, found_pi;
                pointD found_pt;
                double found_d;
                bool hit = find_nearest_point(w, w->series, w->nseries,
                                              mpos, pick_thresh,
                                              &found_si, &found_pi,
                                              &found_pt, &found_d);

                if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) &&
                    wi == active_wi) {
                    if (CheckCollisionPointRec(mpos, w->viewport_px)) {
                        if (hit) {
                            w->selection.active = true;
                            w->selection.widget_idx = wi;
                            w->selection.series_idx = found_si;
                            w->selection.pt_idx = found_pi;
                            w->selection.pt = found_pt;
                        } else {
                            w->selection.active = false;
                        }
                    }
                }

                /* cancel selection with Q */
                if (IsKeyPressed(KEY_Q))
                    w->selection.active = false;

                /* highlight single selection point + tooltip */
                int csx = (int)fmax(0, floor(w->viewport_px.x));
                int csy = (int)fmax(0, floor(w->viewport_px.y));
                int csw = (int)fmax(0, ceil(w->viewport_px.width));
                int csh = (int)fmax(0, ceil(w->viewport_px.height));
                BeginScissorMode(csx, csy, csw, csh); {
                    if (w->selection.active &&
                        w->selection.widget_idx == (size_t)wi) {
                        if (w->selection.series_idx < w->nseries) {
                            series_t *s = &w->series[w->selection.series_idx];
                            Vector2 sp = world_to_screen(w,
                                w->selection.pt.x, w->selection.pt.y);
                            DrawCircleV(sp, 4, WHITE);
                            DrawCircleV(sp, 3, s->color);
                            char sbuf[128];
                            snprintf(sbuf, sizeof(sbuf), "%s\nx=%.6g\ny=%.6g",
                                     s->name,
                                     w->selection.pt.x, w->selection.pt.y);
                            int tw = MeasureTextEx(GetFontDefault(),
                                                   sbuf, 12, 1).x;
                            float bx = sp.x + 10, by = sp.y - 6;
                            DrawRectangle((int)bx, (int)by, tw + 8, 48,
                                          (Color){255, 255, 255, 220});
                            DrawText(sbuf, (int)bx + 4, (int)by + 4,
                                     12, DARKGRAY);
                        }
                    }
                    /* mouse tooltip */
                    draw_mouse_tooltip(w, w->series, w->nseries, mpos);
                } EndScissorMode();

                /* pan hint text while dragging */
                if (IsMouseButtonDown(MOUSE_MIDDLE_BUTTON)) {
                    char lab[64] = {};
                    if (mod_shift() && !mod_ctrl())
                        snprintf(lab, sizeof(lab), "Pan: X only");
                    else if (mod_ctrl() && !mod_shift())
                        snprintf(lab, sizeof(lab), "Pan: Y only");
                    else
                        snprintf(lab, sizeof(lab), "Pan: XY");
                    DrawText(lab, (int)w->viewport_px.x + 8,
                             (int)w->viewport_px.y + 8, 10, DARKGRAY);
                }
            }

            /* f / F : auto-fit view */
            if (IsKeyPressed(KEY_F)) {
                int start = 0, end = (int)widgets.count;
                if (!mod_shift() && active_wi >= 0) {
                    start = active_wi;
                    end   = active_wi + 1;
                }
                for (int fi = start; fi < end; ++fi)
                    auto_fit_view(&widgets.items[fi]);
            }

            /* TAB : toggle fullscreen on active widget */
            if (IsKeyPressed(KEY_TAB)) {
                if (fullscreen_idx >= 0) {
                    for (int fi = 0; fi < (int)widgets.count; ++fi)
                        widgets.items[fi].viewport_px = widgets.items[fi].saved_viewport;
                    fullscreen_idx = -1;
                } else if (active_wi >= 0) {
                    for (int fi = 0; fi < (int)widgets.count; ++fi)
                        widgets.items[fi].saved_viewport = widgets.items[fi].viewport_px;
                    widgets.items[active_wi].viewport_px = (Rectangle){
                        55, 25, screenW - 65, screenH - 42
                    };
                    fullscreen_idx = active_wi;
                }
            }

        } EndDrawing();
    }

    /* ---------------------------------------------------------------- */
    /*  cleanup                                                         */
    /* ---------------------------------------------------------------- */
    for (size_t wi = 0; wi < widgets.count; ++wi) {
        plot_widget_t *w = &widgets.items[wi];
        for (size_t si = 0; si < w->nseries; ++si)
            free(w->series[si].pts);
        free(w->series);
        free(w->multi.series_idx);
        free(w->multi.pt_idx);
    }
    free(widgets.items);

    CloseWindow();
    return 0;
}
