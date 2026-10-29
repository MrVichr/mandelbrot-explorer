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
directory.  VS Code: `.vscode/c_cpp_properties.json` has a "Qt6 GUI" and a
"Headless (fake_qt)" IntelliSense configuration (gnu++20); `launch.json` has
gdb configs for the GUI and the runner, with build tasks in `tasks.json`.

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
make -f Makefile-headless check-all         # every fractal x mode x precision
make -f Makefile-headless bless-all         # record all references
make -f Makefile-headless bench             # timing of a larger Mandelbrot view (double, dd, ddheap)
```

`CXX` defaults to `g++-11`; override with `CXX=g++` where that is not
installed.  `check-*` targets exit 0 with `PASS` or `SKIP` (no blessed copy
yet) and non-zero on `FAIL`.  Single-case targets take `PREC=double|float128|dd|qd|real642|ddheap`.
`check-all` prints one line per case, full log in `current/check-all.log`;
Mandelbrot float128/real642 cases take minutes (soft-float).

Both builds use `-ffp-contract=off` so the compiler cannot fuse `a*b+c` into
FMA; the double-double error-free transforms (`two_prod`, `two_sum`) and
bit-exact regression checks depend on that.  Explicit `std::fma` is unaffected:
both builds also pass `-mfma`, which defines `FP_FAST_FMA` and makes
`two_prod`/`two_sqr` use `std::fma` (one instruction), and `mul_double`/`mul`/`sqr`
use the FMA algorithms of Joldes, Muller, Popescu 2017 ("Tight and rigorous
error bounds for basic building blocks of double-word arithmetic":
DWTimesFP3, DWTimesDW3).  Without `-mfma` they fall back to Dekker's split and
the QD formulas; dd/qd results then differ in the last bits (1 Laguerre pixel
in the check views), so references are build-specific.  An `-mfma` binary
needs an FMA CPU (Haswell / Piledriver or newer).

Pieces of the headless build:

| File | Role |
|---|---|
| `fake_qt/` | Stub Qt headers. `QThread` is real (`std::thread` + task queue, `postTask`, `setCurrentThread` for the main thread) — **do not revert to stubs**. `QObject` has real signal/slot dispatch with thread affinity (same thread → direct call, cross-thread → `postTask`); `~QObject` removes all connections where it is sender or receiver (as real Qt does — without it, switching precision was a use-after-free). `QReadWriteLock` = `std::shared_mutex`, `QImage` = `std::vector<QRgb>`, `QElapsedTimer` = `std::chrono`. |
| `runner/moc_signals.cpp` | Signal bodies for `MandelEvaluatorThread`, `MandelModel`, `JuliaModel`, `LaguerreModel`. |
| `runner/main_runner.cpp` | Headless driver. Options: `--fractal man\|jul\|lag`, `--mode text\|images`, `--cx`, `--cy`, `--scale`, `--width`, `--height`, `--title`, `--bless`, `--precision double\|float128\|dd\|qd\|real642\|ddheap`, `--c-re`, `--c-im`, `--period` (Julia/Laguerre), `--point`. `--cx/--cy/--c-re/--c-im` are decimal strings parsed in the selected precision (Makefile-headless passes Julia/Laguerre c as exact binary64 expansions so the references stay valid). `--point` evaluates only (`--cx`, `--cy`) as a 1x1 image at scale 2^-300 (so the centre is not rounded to the pixel grid), skips `paintOrbit`, and prints the pixel's store and point values; launch config "Headless: single point" debugs it. Text mode writes per-pixel JSON; images mode writes one PPM per paint style. Output names carry a precision tag (e.g. `Man_dd_...`) except for double. Prints `timing:` to stderr. The precision-change slot is queued, so the runner calls `processEvents()` right after `setselectedPrecision`. |
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
| `MandelMath.hpp` / `MandelMath_var.cpp` | Number types: `number<BASE>` for `double`, `__float128`, `dd_real`, `dq_real`, `real642`, `dd_heap`, and `number_any` (variant wrapper); `complex<BASE>`; `number_traits<BASE>` |
| `dd_heap.hpp` | `dd_heap`: dd_real kept on the heap with a per-thread free list; the template for future heap-backed types |
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

## Number design (current)

- Value semantics, no scratchpad.  `number<BASE>` is one generic template
  holding only `BASE store`; every operation forwards to
  `number_traits<BASE>` (specialised for `double`, `__float128`, `dd_real`,
  `dq_real`, `real642`, `dd_heap`), all inline in `MandelMath.hpp`.  The five
  POD-backed `number<BASE>` and `complex<BASE>` are trivially copyable (static_asserted),
  so temporaries are free.
- `number<number_any>` is a full class specialisation over a `std::variant`,
  with out-of-line members in `MandelMath_var.cpp`.  It is used by the models
  and by the main-thread `MandelEvaluator<number_any>` in each Orbit, so it
  must support the same operator set as the concrete types.  Its move leaves
  the source empty (swap); `complex` deliberately has no implicit move so a
  moved-from `complex<number_any>` keeps its value.
- Operators come from the CRTP base `number_operators<N>` (hidden friends):
  `+= -= *=` (N or double), binary `+ - *`, unary `-`, `==`, and free
  `sqr`, `ldexp`, `recip`, `sqrt`.  `complex` has the same plus by-value
  `getMag()`, `getMag1()`, `getDist1()`, `mulreT()`, `ccw()`, `dist2()`,
  `recip_prepared(mag)`.  Old method names (`add`, `mul`, ...) remain as
  in-place wrappers.
- Hot double-double primitives (`two_sum`, `two_prod`, `add`, `mul`, `sqr`,
  `recip`, `sqrt`, compares, `ldexp2`) are inline in `double_double.hpp`.
- Rule for any future heap-backed BASE (e.g. MPFR): follow `dd_heap`
  (`dd_heap.hpp`).  It pools its storage internally (thread_local free list;
  assignment copies the value so long-lived objects keep their buffer; moved-from
  objects hold nullptr and may only be assigned or destroyed) and declares a POD
  storage form via `number_storage<BASE>` (dd_real for dd_heap), used by
  `readFrom`/`writeTo`/`convert_block` and the models' `points` arrays
  (`delete[] (dd_real *)points`).  Results of `ddheap` are bit-identical to `dd`;
  the Mandelbrot check view is ~17% slower.
- Adding a precision touches: `NumberType`, `NumberTypeFromBase`, the
  `number_any` variant, `number_traits`, every `Assigner` and `convert_block`
  in `MandelMath_var.cpp`, `MandelEvaluatorThread::owner`, the models'
  precision enums / `threads` variants / point-array switches, `main.qml`,
  `PRECISION_NAMES` in the runner, `check_images.py` and `PRECISIONS` in
  Makefile-headless.  `MandelMath::visit` handles up to 8 alternatives.
- Constructors still take `NumberType` (`number(NumberType)`,
  `complex(NumberType)`, `MandelPoint`, `JuliaPoint`, `LaguerrePoint`,
  `Compute*Params`).
- Decimal I/O: `toString()` returns the shortest decimal that parses back to
  the same value (double: `std::to_chars`; float128/dd/qd/ddheap: libquadmath
  via `float128_toString`; real642 prints `val2`); `fromString(QString)` parses
  at full precision (dd from a float128 parse, so 113 bits) and returns false on
  garbage.  Both are independent of `LC_NUMERIC` (Qt sets it from the
  environment).  Both builds link `-lquadmath`.
- GUI view presets in `main.qml` are strings (`viewRe`, `viewIm`, `cRe`, `cIm`;
  `viewZoom` stays a number), applied via `MandelModel::setView_str` /
  `makeViewInfo`.  `getViewJson()` on each model gives the current view in the
  same format; the popup menu's "Copy view" puts it on the clipboard
  (`MandelImageCombiner::copyToClipboard`).
- float128 note: the old `(uint16_t*)&s` exponent hack in `lshift` violated
  strict aliasing and was miscompiled once inlined (`lshift_float128` now
  uses `memcpy`).  Older float128 references were produced by that
  miscompile; current references were blessed from the original sources
  built with `-fno-strict-aliasing`.

---

## To do

1. **Finish the math-style migration.**  The hot loops (`MandelLoopEvaluator`,
   `*_until_bailout`, evaluateMandel/Julia steps, `LaguerreStep`) use
   operators; the rest of the evaluator and the models still use the in-place
   method style and can be converted gradually.  Numeric follow-ups (change
   results, so re-bless): `operator/` where a reciprocal is not reused.
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
