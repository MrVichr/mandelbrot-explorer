# AGENTS.md — Mandelbrot Explorer (headless C++)

This file is the primary orientation document for any AI agent or developer
picking up work across sessions.

---

## What this project is

A Mandelbrot / Julia / Laguerre fractal explorer, originally a Qt/QML desktop
app (`Mandelbrot.pro`).  The goal is to **run its computation core headlessly**,
without Qt installed, so we can exercise and test the math directly from a
plain C++ binary.

---

## Current state (as of 2026-03)

### What is DONE — everything needed to run headlessly

| Area | Status |
|---|---|
| `fake_qt/` headers | All Qt types compile cleanly (`make` builds without errors). |
| `QThread` | Backed by `std::thread` + mutex/condition-var task queue. `start()` spawns OS thread. `exec()` drains task queue. `postTask(fn)` enqueues work. `wait()` joins. `currentThread()` works. `setCurrentThread(t)` lets main thread register itself without spawning an OS thread. |
| `QReadWriteLock` | `std::shared_mutex`-backed. |
| `QImage` | Real `std::vector<QRgb>` pixel buffer. `isNull()`, `width()`, `height()`, `pixel()`, `setPixel()`, `bits()`, `scanLine()` all functional. |
| `QElapsedTimer` | `std::chrono`-backed. |
| Signal/slot dispatch | `QObject::connect` registers `FakeConnEntry`. `emitSignal` dispatches: same thread → direct call; cross-thread → `targetThread->postTask(lambda)`. `fakeSubtuple<SArgs...>` handles slots with fewer params than signal. |
| `runner/moc_signals.cpp` | Signal body definitions for 4 classes: `MandelEvaluatorThread` (3 signals), `MandelModel` (4), `JuliaModel` (7), `LaguerreModel` (3). |
| `Makefile` | `make` → `./mandel_runner`. All 11 project files + runner files. Regression targets: `bless-text-{man,jul,lag}` / `check-text-{man,jul,lag}` write/compare per-pixel JSON; `bless-images-{man,jul,lag}` / `check-images-{man,jul,lag}` write/compare per-style PPMs.  `check-*` targets exit 0 with `SKIP` summary when no blessed artefacts exist yet. |
| `runner/main_runner.cpp` | Headless test driver. Creates fake main thread, `MandelModelInspector` (exposes protected members), calls `setImageSize` then `setView_double`, pumps event loop, prints per-pixel iter counts + rstate + ARGB colours. Snapshot regression: `--bless` writes `runner/mandel_runner.snap`; normal run compares and exits 1 on mismatch. |
| End-to-end computation | **Working.** Snapshot test passes (`./mandel_runner` exits 0). Symmetric Mandelbrot at 40×20 with correct iteration counts, rstate categories, and ARGB colours. |

### Critical usage note

**`setImageSize(W, H)` must be called BEFORE `setView_double(...)`.**

`MandelModel::transformStore` has:
```cpp
for (int newy=(new_height-1)&0xfffffff; newy>=0; newy--)  // ~268M iters if new_height==0!
```
Calling `setView_double` with `imageHeight==0` causes ~268 million no-op iterations (runtime: ~1 second). Always size the image first.

---

## Scratchpad refactoring — current state and remaining work

The original design passed a `Scratchpad *spad` through every constructor so
that `number<BASE>` objects could find their scratch temporaries.  This caused
widespread constructor-parameter boilerplate and left `number_any` variant
members monostate (uninitialised) whenever a `Scratchpad *` was not available.

### What has been done

| Change | Where |
|---|---|
| `tl_spad` is `static thread_local Scratchpad*` per concrete type | `MandelMath_var.cpp` (one definition per type, lines ~4659) |
| `set_spad` for concrete types is **out-of-line** in `MandelMath_var.cpp` | Explicit specialisation declarations in `MandelMath.hpp` lines ~354 suppress inline copies in other TUs; this fixes a GCC TLS-wrapper SIGSEGV on first worker-thread access |
| `set_spad(&tmp)` called at the top of all three `do*ThreadedSplit` functions | `MandelEvaluator.cpp` — ensures each worker OS thread registers its scratchpad before any `number_any` construction |
| `complex(Scratchpad *spad)` initialises `re(spad), im(spad)` | `MandelMath.hpp` line ~535 — ensures `number<number_any>(spad)` runs for both components, emplacing the variant |
| `MandelPoint`, `JuliaPoint`, `LaguerrePoint` each have a `(Store*, Scratchpad*)` constructor overload | `MandelEvaluator.hpp/.cpp` — all `complex`/`number` member fields initialised with `spad`, so they are never monostate |
| `wtiPoint` in all three `PrecisionRecord` constructors uses the spad overload | `MandelModel.cpp`, `JuliaModel.cpp`, `LaguerreModel.cpp` — passes `&orbit.evaluator.tmp` which is already constructed at that point |

### Testing improvements — done

1. **CLI parameters for position/zoom/size** — `runner/main_runner.cpp` accepts
   `--cx`, `--cy`, `--scale`, `--width`, `--height`, `--title`, `--fractal`,
   `--mode`, and `--bless`.

2. **One PPM per paintStyle** — `--mode images` cycles through all paint-style
   values, calls `writeToImage`, and writes PPMs (+ PNGs via Pillow) into
   `current/` or `blessed/`.

3. **Python regression wrapper** — `runner/check_images.py` drives
   `mandel_runner`, compares output against `blessed/`, and reports
   `PASS` / `FAIL` / `SKIP` (SKIP when no blessed artefact exists yet).

4. **Makefile regression targets** — `bless-text-*` / `check-text-*` for JSON
   output; `bless-images-*` / `check-images-*` for PPM image output.  All
   `check-*` targets exit 0 with a `SKIP` summary when no blessed copy exists.

### What still needs work (scratchpad refactoring)

1. **`number<BASE>` still carries `Scratchpad *tmp`** — the 8-byte member
   pointer exists on every `number<double>` instance.  Now that `tl_spad` is
   `thread_local`, this pointer is redundant for all operations that already
   call `set_spad`.  Removing it would halve the size of `number<double>` and
   eliminate the remaining constructor plumbing.  The main obstacle is that
   `Scratchpad` itself contains `number<BASE>` members (`tmp1`–`tmp4`), so the
   two types are mutually recursive; careful ordering is needed.

2. **`Scratchpad` should be static**
   - set_spad should only be called near MandelEvaluator::moveToThread.

3. **`Scratchpad` constructor parameter still threads through many call sites**
   — `PrecisionRecord` members (`position`, `params`, `lagu_c`, etc.) all
   accept `Scratchpad *` in their constructors.  Once `tl_spad` is fully
   trusted these could use `tl_spad` directly and drop the parameter.

4. **`Scratchpad::ntype` is redundant** — with `NUMBER_IS_VARIANT=1` the
   variant's `.index()` already encodes the type.  The field can be removed
   once callers that read `spad->ntype` are updated to use
   `number<number_any>::get_spad()->ntype` or the variant index instead.

5. **Delete the `NUMBER_IS_VARIANT=0` branch** — two parallel implementations
   exist (variant vs. virtual dispatch).  The virtual dispatch branch is dead
   code; delete everything inside the `#else` block down to `#endif`.

6. **Return by value from `_tmp()` methods** — `getMag_tmp()` etc. return
   `const number<BASE>*` into a scratch pool.  NRVO makes returning by value
   safe and eliminates pointer-lifetime questions at every call site.

---

## Key files

| File | Role |
|---|---|
| `MandelModel.cpp/hpp` | Top-level model: owns image store, worker threads, view position, epoch counter |
| `MandelEvaluator.cpp/hpp` | Per-thread worker; `MandelEvaluatorThread` (QThread subclass) + `MandelEvaluator<BASE>` template |
| `MandelMath.hpp` / `MandelMath_var.cpp` | Number types: `number<double>`, `number<dd_real>`, `number_any` (variant wrapper) |
| `double_double.cpp/hpp` | `dd_real` / `dq_real` extended-precision arithmetic |
| `fake_qt/QThread` | Real `std::thread`-backed QThread — **DO NOT revert to stubs** |
| `fake_qt/QObject` | Real signal/slot implementation with thread-affinity dispatch |
| `fake_qt/QReadWriteLock` | `std::shared_mutex`-backed |
| `fake_qt/QImage` | Real pixel-buffer backing (`std::vector<QRgb>`) |
| `runner/moc_signals.cpp` | Signal body definitions (compile as part of `mandel_runner`) |
| `runner/main_runner.cpp` | Headless test driver — **primary way to exercise the engine** |
| `runner/check_images.py` | Python regression wrapper: drives `mandel_runner`, compares against `blessed/`, reports PASS/FAIL/SKIP |
| `Makefile` | `make` builds `./mandel_runner`; `check-*` / `bless-*` targets run regression checks |

---

## Threading architecture (important context)

```
MandelModel (main thread)
  │
  ├─ MandelModel::startNewEpoch()
  │    └─ emit triggerComputeThreaded(epoch)          [QueuedConnection →]
  │                                                      MandelEvaluatorThread::doMandelThreaded(epoch)
  │                                                         runs on worker thread OS thread
  │                                                         MandelEvaluator<T>::doMandel_threaded()
  │                                                         emit thread.doneMandelThreaded(&thread) [QueuedConnection →]
  │                                                                                MandelModel::doneWorkInThread()
  │                                                                                (runs on main thread)
  └─ MandelModel::doneWorkInThread()
       giveWorkThreaded() / triggers next epoch
```

`MandelEvaluatorThread` has `moveToThread(&thread)` called on it, so it is
"owned" by its own OS thread.  Signals into it are `QueuedConnection` →
posted to its task queue.  `MandelModel` lives on the main thread; signals
back to it are also `QueuedConnection` → posted to the main thread's queue.

The main thread must therefore have its own `QThread` wrapper (or equivalent
`processEvents()` capability) to drain those posted slots.

---

## How to run the regression checks

```bash
cd /workspaces/LNX/MAND/mandelbrot-explorer
make check-text-man     # compare Mandelbrot JSON against blessed/
make check-images-man   # compare Mandelbrot PPMs against blessed/
make check-text-jul     # Julia text check
make check-images-jul   # Julia image check
make check-text-lag     # Laguerre text check
make check-images-lag   # Laguerre image check
```

Targets exit 0 (`PASS` or `SKIP`) and non-zero on `FAIL`.
Run the corresponding `bless-*` target first to record reference output.
