# plot_widget_rl

Interactive 2D plot widget system built with [raylib](https://www.raylib.com/). Provides real-time, zoomable and pannable graphing designed for scientific computing and neural simulation visualization.

## Features

- **Multiple plot widgets** in a single window, each with independent viewport, title, and series
- **Per-widget input** — only the widget under the mouse cursor (the *active* widget) responds to pan, zoom, click, and key shortcuts; the active widget is highlighted with a yellow border
- **Pan & zoom** via middle-mouse drag and scroll wheel, with axis-lock modifiers (Shift = X only, Ctrl = Y only)
- **Adaptive grid and ticks** — tick spacing automatically adjusts using "nice" step values (1, 2, 5 × 10^n); labels use compact `%g` formatting
- **Cursor tooltip** — hover to see world coordinates and nearest point values for every visible series
- **Point selection** — left-click on any data point to highlight it and display its coordinates; press Q to deselect
- **Area selection** — left-click drag to select multiple points within a rectangle; highlighted in the viewport
- **Legend** with per-series visibility toggle (click the series name)
- **NaN-safe rendering** — gaps in data (NaN/Inf) are skipped gracefully
- **Clipping** via raylib scissor mode for correct viewport boundaries

## Controls

| Input | Action |
|---|---|
| Middle-mouse drag | Pan (Shift = X only, Ctrl = Y only) |
| Scroll wheel | Zoom centered on cursor (Shift = X only, Ctrl = Y only) |
| Left-click on point | Select / highlight point |
| Left-click drag | Area selection (multi-point) |
| Q | Clear selection / cancel area drag |
| R | Reset pan and zoom on the hovered widget |
| Tab | Toggle fullscreen on the active widget |
| f | Zoom-fill: auto-fit view to data on the hovered widget |
| F (Shift+f) | Zoom-fill all widgets |

## Izhikevich Neuron Models

The project simulates all 20 firing patterns from Izhikevich's 2003 paper
*"Simple model of spiking neurons"* (IEEE Transactions on Neural Networks).
Each pattern is displayed in its own widget with three series: injected
current **I**, membrane potential **V**, and recovery variable **u**.

Models included: tonic spiking, phasic spiking, tonic bursting, phasic
bursting, mixed mode, spike frequency adaptation, Class 1/2 excitable,
spike latency, subthreshold oscillations, resonator, integrator, rebound
spike/burst, threshold variability, bistability, DAP, accommodation,
inhibition-induced spiking/bursting.

## Dependencies

- [raylib](https://www.raylib.com/) (>= 4.0)
- C99 compiler (gcc, clang)

## Build

```sh
make
```

## Source files

- `main.c` — program entry point, widget setup, main loop
- `plot_widget.h` / `plot_widget.c` — reusable plot widget types and rendering
- `izhikevich.h` / `izhikevich.c` — Izhikevich neuron model generators
```

## Run

```sh
./plot_widget_rl
```
