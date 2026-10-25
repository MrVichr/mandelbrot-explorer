# AGENTS.md - Mandelbrot Explorer

Primary orientation document for any AI agent or developer picking up work
across sessions.

---

## What this project is

A Mandelbrot / Julia / Laguerre fractal explorer: a Qt 6 / QML desktop app
(`Mandelbrot.pro`, `main.cpp`, `main.qml`) whose computation core can also be
built **headlessly** against `fake_qt/` stub headers, so the math can be
exercised and regression-tested from a plain C++ binary (`mandel_runner`).

---

## Building

### GUI (real Qt)

`Mandelbrot.pro` with `QT += quick`; application type is `QGuiApplication`.
Build out-of-tree (e.g. `build-gui/`) to keep qmake output out of the source
directory.

`main.qml` uses **QtQuick.Controls 2** menus (not `Qt.labs.platform`): native
platform menus on Wayland caused `setGrabPopup ... does not match the current
topmost grabbing popup` errors and erratic popups.

```bash
cd build-gui
qmake6 ../Mandelbrot.pro
make -j$(nproc)
./Mandelbrot
```

### Headless (`mandel_runner`)

The hand-written headless build lives in **`Makefile-headless`** (not
`Makefile`, so an in-tree qmake run cannot overwrite it).

```bash
make -f Makefile-headless                   # build ./mandel_runner
make -f Makefile-headless check-text-man    # also -jul / -lag
make -f Makefile-headless check-images-man  # bless-* targets record references
```

`CXX` defaults to `g++-11`; override with `CXX=g++` where that is not
installed.  `check-*` targets exit 0 with `PASS` or `SKIP` (no blessed copy
yet) and non-zero on `FAIL`.

Pieces of the headless build:

| File | Role |
|---|---|
| `fake_qt/` | Stub Qt headers. `QThread` is real (`std::thread` + task queue, `postTask`, `setCurrentThread` for the main thread) — **do not revert to stubs**. `QObject` has real signal/slot dispatch with thread affinity (same thread → direct call, cross-thread → `postTask`). `QReadWriteLock` = `std::shared_mutex`, `QImage` = `std::vector<QRgb>`, `QElapsedTimer` = `std::chrono`. |
| `runner/moc_signals.cpp` | Signal bodies for `MandelEvaluatorThread`, `MandelModel`, `JuliaModel`, `LaguerreModel`. |
| `runner/main_runner.cpp` | Headless driver. Options: `--fractal man\|jul\|lag`, `--mode text\|images`, `--cx`, `--cy`, `--scale`, `--width`, `--height`, `--title`, `--bless`. Text mode writes per-pixel JSON; images mode writes one PPM per paint style. |
| `runner/check_images.py` | Regression wrapper: runs `mandel_runner`, compares `current/` against `blessed/`, reports `PASS` / `FAIL` / `SKIP` (no blessed copy yet). See its docstring for usage. |

`blessed/` and `current/` are git-ignored, local only.

### Critical usage note

**`setImageSize(W, H)` must be called BEFORE `setView_double(...)`.**
`MandelModel::transformStore` loops `for (int newy=(new_height-1)&0xfffffff; ...)`,
i.e. ~268M no-op iterations when the height is still 0.

---

## Key files

| File | Role |
|---|---|
| `MandelModel.cpp/hpp`, `JuliaModel.cpp/hpp`, `LaguerreModel.cpp/hpp` | Top-level models: image store, worker threads, view position, epoch counter, paint styles |
| `MandelEvaluator.cpp/hpp` | Per-thread worker: `MandelEvaluatorThread` (QThread subclass) + `MandelEvaluator<BASE>` (`evaluateMandel`, `evaluateJulia`, `mandel_until_bailout`, `julia_until_bailout`, `MandelLoopEvaluator`, `Bulb`, `LaguerreStep`) |
| `MandelMath.hpp` / `MandelMath_var.cpp` | Number types: `number<BASE>` for `double`, `__float128`, `dd_real`, `dq_real`, `real642`, and `number_any` (variant wrapper); `complex<BASE>`; `Scratchpad` |
| `double_double.cpp/hpp`, `multiprec.cpp/hpp` | Extended-precision arithmetic |
| `MandelImageCombiner.cpp/hpp`, `ShareableImageWrapper.cpp/hpp`, `atomic_update_rect.cpp/hpp` | GUI-side image composition |

---

## Threading architecture

```
MandelModel (main thread)
  ├─ startNewEpoch()
  │    └─ emit triggerComputeThreaded(epoch)      [Queued →]
  │         MandelEvaluatorThread::doMandelThreaded(epoch)   (worker OS thread)
  │           MandelEvaluator<T>::doMandel_threaded()
  │           emit thread.doneMandelThreaded(&thread)        [Queued →]
  └─ doneWorkInThread()   (main thread) → giveWorkThreaded() / next epoch
```

`MandelEvaluatorThread` is `moveToThread`'d onto its own thread; the main
thread needs its own event loop (or a `QThread` wrapper in the headless build)
to drain the queued replies.

---

## Number / Scratchpad design (current)

- `number<BASE>` holds only `BASE store`; `number<number_any>` is a
  `std::variant`.  Constructors take `NumberType` (`number(NumberType)`,
  `complex(NumberType)`, and likewise `MandelPoint`, `JuliaPoint`,
  `LaguerrePoint`, `ComputeMandelParams`, `ComputeJuliaParams`,
  `ComputeLaguerreParams`).
- Temporaries come from a per-type `static thread_local Scratchpad *tl_spad`,
  registered with `number<BASE>::set_spad(ntype, &tmp)` at worker-thread
  startup (`MandelEvaluator.cpp`) and in `Scratchpad` constructors.
- `get_spad()` is per-instance: inlined `tl_spad` for concrete types; for
  `number_any` it dispatches on the variant index and follows the concrete
  spad's `inner` back-pointer.
- **GCC TLS pitfall**: `set_spad`/`get_spad` specialisations must stay
  declared in `MandelMath.hpp` and defined out-of-line in `MandelMath_var.cpp`;
  inline copies produce a TLS wrapper with a null init pointer → SIGSEGV.
- **Known limitation**: two models with the same concrete precision on the
  main thread share `number<number_any>::tl_spad` (the last constructed wins).

---

## To do

1. **Refactor how math is written.**  Arithmetic currently goes through
   explicit method calls and `*_tmp()` scratch pointers (e.g. `getMag_tmp()`
   returns a pointer into the scratchpad).  Consider operator overloading
   (`+ - * /`, compound assignment) and returning by value (NRVO) to make
   formulas readable and remove pointer-lifetime questions.  Leftovers from
   the scratchpad cleanup that fit into this work:
   - `Scratchpad::ntype` is redundant with the variant index (still read via
     `->ntype` in the models, `MandelMath_var.cpp` and `MandelEvaluator.cpp`).
   - `MandelEvaluator.hpp` still keeps a `number<number_any>::Scratchpad *spad`.
2. **Speed up initial iterations with a bilinear approximation.**  Determine
   how many initial iterations can be skipped by approximating the orbit with
   a bilinear map (series-approximation style), with an error bound that
   decides when to switch to exact iteration.
3. **Float + GPU.**  If (2) works, move the bulk of the per-pixel iteration
   to `float` (perturbation around a high-precision reference) and run it on
   the GPU.
4. **Julia interior distance estimate is wrong.**  The Wikipedia formulas
   do not give correct results; local patches did not fix it.  Code:
   `MandelEvaluator<BASE>::evaluateJulia()` (`interior.inte`, `phi`, around
   `juliaData.store->interior.*` in `MandelEvaluator.cpp`), and the hover
   readout / paint styles `Phi`, `Phi1`, `Phi2`, `F1`, `FirstUnder1` in
   `JuliaModel.cpp`.
