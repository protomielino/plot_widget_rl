# plot_widget_rl

Interactive 2D plot widget system built with [raylib](https://www.raylib.com/).

## Features

- **Multiple plot widgets** in a single window, each with independent viewport, title, and data series
- **Move widgets** by click+drag on the title bar
- **Resize widgets** by right-click+drag on the bottom-right corner handle
- **Redistribute layout** with G key — auto-grids all widgets filling the window
- **Pan & zoom** via middle-mouse drag and scroll wheel, with axis-lock modifiers (Shift = X only, Ctrl = Y only)
- **Adaptive grid and ticks** — tick spacing auto-adjusts; compact `%g` labels
- **Cursor tooltip** — hover for world coordinates and nearest point values per series
- **Point selection** — left-click on a data point to highlight it; Q to deselect
- **Area selection** — left-click drag to select multiple points in a rectangle
- **Legend** with per-series visibility toggle
- **NaN-safe rendering** — gaps in data (NaN/Inf) are skipped
- **Scissor-clipped** viewports for correct boundaries

## Controls

| Input | Action |
|---|---|
| Left-click drag on title bar | Move widget |
| Right-click drag on corner handle | Resize widget |
| Middle-mouse drag | Pan (Shift = X only, Ctrl = Y only) |
| Scroll wheel | Zoom centered on cursor (Shift = X only, Ctrl = Y only) |
| Left-click on point | Select / highlight point |
| Left-click drag on plot | Area selection (multi-point) |
| G | Redistribute widgets in a grid |
| Q | Clear selection / cancel area drag |
| R | Reset pan and zoom on hovered widget |
| Tab | Toggle fullscreen on active widget |
| F | Auto-fit view to data on hovered widget |
| Shift+F | Auto-fit all widgets |

## Widget data

Four widgets show simple mathematical functions (256 points each):

- **Sine & Cosine** — sin(x), cos(x), sin(x)·cos(x)
- **Polynomial** — x²−4, x³/10, x
- **Damped Osc.** — e⁻ˣ/⁵ sin(x), e⁻ˣ/¹⁰ cos(x), e⁻ˣ/⁷ sin(2x)
- **Trigonometric** — sin(2x), sin(3x), cos(x)·sin(2x)

## Dependencies

- [raylib](https://www.raylib.com/) (>= 4.0)
- C99 compiler (gcc, clang)

## Build & Run

```sh
make
./plot_widget_rl
```
