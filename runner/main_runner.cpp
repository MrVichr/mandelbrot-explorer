// runner/main_runner.cpp — headless test driver for MandelModel / JuliaModel / LaguerreModel.
//
// Modes (--mode):
//   text    — Compute the fractal, then write a JSON file to the output
//             directory with per-point data and colour values (one per paint style).
//             Also prints a VT100 render of paintStyleCls to stdout.
//
//   images  — Compute the fractal, then write one PPM per paint style to the
//             output directory.
//
// Output directories (set by --bless):
//   default: ./current/
//   --bless: ./blessed/
//
// CLI options:
//   --fractal man|jul|lag  Fractal type             (default man)
//   --cx <num>             Centre real part         (default -0.5)
//   --cy <num>             Centre imaginary part    (default 0.0)
//   --scale <float>        Scale                    (default 0.125)
//   --width  <int>         Image width in pixels    (default 40)
//   --height <int>         Image height in pixels   (default 20)
//   --title <str>          Label                    (default "Wide view")
//   --c-re <num>           Julia/Laguerre parameter c, real part
//   --c-im <num>           Julia/Laguerre parameter c, imaginary part
//   --period <int>         Julia/Laguerre period
//   --precision double|float128|dd|qd|real642|ddheap (default double)
//   --mode text|images                              (default text)
//   --bless                Write to ./blessed/ instead of ./current/
//   --point                Evaluate only the point (--cx, --cy) and print the result;
//                          ignores --scale/--width/--height/--mode, writes no files
//
// <num> values are decimal strings parsed in the selected precision, so digits
// beyond double are kept.
//
// Generated filenames are prefixed: Man_..., Jul_..., Lag_...; non-double
// precisions add a tag, e.g. Man_dd_...
// Setup and compute wall times are printed to stderr.
//
// Build: make  (from the mandelbrot-explorer directory)
// NOTE: setImageSize must be called BEFORE setView_double / setView (see AGENTS.md).

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cctype>
#include <chrono>
#include <string>
#include <thread>
#include <vector>
#include <sys/stat.h>

#include "fake_qt/QThread"

#include "MandelModel.hpp"
#include "JuliaModel.hpp"
#include "LaguerreModel.hpp"
#include "ShareableImageWrapper.hpp"

// ── Expose protected members ─────────────────────────────────────────────────
struct MandelModelInspector : public MandelModel {
    const MandelPointStore *store()  const { return pointStore; }
    int width()  const { return imageWidth; }
    int height() const { return imageHeight; }
    QString center() const { return precisionRecord->position.center.toString(); }
    MandelPoint<MandelMath::number_any> &point(int index) {
        precisionRecord->wtiPoint.readFrom(precisionRecord->points, index*MandelPoint<MandelMath::number_any>::LEN);
        return precisionRecord->wtiPoint;
    }
};

struct JuliaModelInspector : public JuliaModel {
    const JuliaPointStore *store()  const { return pointStore; }
    int width()  const { return imageWidth; }
    int height() const { return imageHeight; }
    QString center() const { return precisionRecord->position.center.toString(); }
    QString c() const { return precisionRecord->params.c.toString(); }
    JuliaPoint<MandelMath::number_any> &point(int index) {
        precisionRecord->wtiPoint.readFrom(precisionRecord->points, index*JuliaPoint<MandelMath::number_any>::LEN);
        return precisionRecord->wtiPoint;
    }
    void setView_str(const char *cx, const char *cy, double scale) {
        MandelMath::complex<MandelMath::number_any> center(precisionRecord->orbit.evaluator.ntype);
        center.fromString(cx, cy);
        setView(center, scale);
    }
};

struct LaguerreModelInspector : public LaguerreModel {
    const LaguerrePointStore *store()  const { return pointStore; }
    int width()  const { return imageWidth; }
    int height() const { return imageHeight; }
    QString center() const { return precisionRecord->position.center.toString(); }
    QString c() const { return precisionRecord->params.c.toString(); }
    LaguerrePoint<MandelMath::number_any> &point(int index) {
        precisionRecord->wtiPoint.readFrom(precisionRecord->points, index*LaguerrePoint<MandelMath::number_any>::LEN);
        return precisionRecord->wtiPoint;
    }
    void setView_str(const char *cx, const char *cy, double scale) {
        MandelMath::complex<MandelMath::number_any> center(precisionRecord->orbit.evaluator.ntype);
        center.fromString(cx, cy);
        setView(center, scale);
    }
};

// ── Main-thread QThread ───────────────────────────────────────────────────────
struct MainQThread : public QThread {
    MainQThread() {
        QThread::setCurrentThread(this);
        fakeMainThread() = this;
    }
    void start() {}
};

// ── Paint-style tables ────────────────────────────────────────────────────────
struct PaintStyleEntry {
    int         value;   // stored as int to be fractal-agnostic
    const char *name;
};

// Mandelbrot — 8 styles
static const PaintStyleEntry MAN_PAINT_STYLES[] = {
    { MandelModel::paintStyleKind,       "Kind"      },
    { MandelModel::paintStyleCls,        "Cls"       },
    { MandelModel::paintStyleExter,      "Exter"     },
    { MandelModel::paintStyleExterAngle, "ExterAngle"},
    { MandelModel::paintStyleInter,      "Inter"     },
    { MandelModel::paintStyleNear,       "Near"      },
    { MandelModel::paintStyleFZ,         "FZ"        },
    { MandelModel::paintStyleFC,         "FC"        },
};
static constexpr int N_MAN_STYLES = 8;
static constexpr int MAN_CLS_IDX  = 1;

// Julia — 13 styles
static const PaintStyleEntry JUL_PAINT_STYLES[] = {
    { JuliaModel::paintStyleKind,         "Kind"       },
    { JuliaModel::paintStyleCls,          "Cls"        },
    { JuliaModel::paintStyleExter,        "Exter"      },
    { JuliaModel::paintStyleExterAngle,   "ExterAngle" },
    { JuliaModel::paintStyleInter,        "Inter"      },
    { JuliaModel::paintStyleNear0,        "Near0"      },
    { JuliaModel::paintStyleNearZ,        "NearZ"      },
    { JuliaModel::paintStylePhi,          "Phi"        },
    { JuliaModel::paintStylePhi1,         "Phi1"       },
    { JuliaModel::paintStylePhi2,         "Phi2"       },
    { JuliaModel::paintStyleF1,           "F1"         },
    { JuliaModel::paintStyleFirstUnder1,  "FirstUnder1"},
    { JuliaModel::paintStyleNearMR,       "NearMR"     },
};
static constexpr int N_JUL_STYLES = 13;
static constexpr int JUL_CLS_IDX  = 1;

// Laguerre — 4 styles
static const PaintStyleEntry LAG_PAINT_STYLES[] = {
    { LaguerreModel::paintStyleCls,    "Cls"    },
    { LaguerreModel::paintStyleNthFz,  "NthFz"  },
    { LaguerreModel::paintStyleNthFz1, "NthFz1" },
    { LaguerreModel::paintStyleRoots,  "Roots"  },
};
static constexpr int N_LAG_STYLES = 4;
static constexpr int LAG_CLS_IDX  = 0;

// ── Run options shared by the compute functions ──────────────────────────────

// Same order as the models' precision enums.
static const char *const PRECISION_NAMES[] = { "double", "float128", "dd", "qd", "real642", "ddheap" };
static constexpr int N_PRECISIONS = 6;

static int    g_precision = 0;
static const char *g_cx = "-0.5";
static const char *g_cy = "0";
static const char *g_c_re = nullptr;
static const char *g_c_im = nullptr;
static int    g_period = 0;
static bool   g_point = false;

// ── Helpers ──────────────────────────────────────────────────────────────────

static std::string file_prefix(const char *fractal_prefix)
{
    std::string s(fractal_prefix);
    if (g_precision != 0)
        s += std::string(PRECISION_NAMES[g_precision]) + "_";
    return s;
}

template <typename MODEL>
static void select_precision(MODEL &model, MainQThread &mainThread)
{
    if (g_precision != 0)
    {
        model.setselectedPrecision(static_cast<typename MODEL::precision>(g_precision));
        // The change slot is queued; run it now, before the old workers get any work.
        mainThread.processEvents();
    }
}

using Clock = std::chrono::steady_clock;

static double ms_since(Clock::time_point t0)
{
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

static void report_timing(double setup_ms, double compute_ms)
{
    std::fprintf(stderr, "timing: setup %.1f ms, compute %.1f ms\n", setup_ms, compute_ms);
}

static std::string sanitize_title(const char *title)
{
    std::string s(title);
    for (char &c : s)
        if (!std::isalnum(static_cast<unsigned char>(c)))
            c = '_';
    return s;
}

static const char *rstate_name(MandelPointStore::ResultState rs)
{
    switch (rs) {
        case MandelPointStore::ResultState::stUnknown:  return "rU";
        case MandelPointStore::ResultState::stOutside:  return "rO";
        case MandelPointStore::ResultState::stOutAngle: return "rA";
        case MandelPointStore::ResultState::stBoundary: return "rB";
        case MandelPointStore::ResultState::stMisiur:   return "rS";
        case MandelPointStore::ResultState::stDiverge:  return "rD";
        case MandelPointStore::ResultState::stPeriod2:  return "r2";
        case MandelPointStore::ResultState::stPeriod3:  return "r3";
        case MandelPointStore::ResultState::stMaxIter:  return "rM";
    }
    return "r?";
}

static const char *rstate_name(JuliaPointStore::ResultState rs)
{
    switch (rs) {
        case JuliaPointStore::ResultState::stUnknown:  return "rU";
        case JuliaPointStore::ResultState::stOutside:  return "rO";
        case JuliaPointStore::ResultState::stOutAngle: return "rA";
        case JuliaPointStore::ResultState::stBoundary: return "rB";
        case JuliaPointStore::ResultState::stMisiur:   return "rS";
        case JuliaPointStore::ResultState::stDiverge:  return "rD";
        case JuliaPointStore::ResultState::stPeriod2:  return "r2";
        case JuliaPointStore::ResultState::stPeriodBoundary: return "rI";
        case JuliaPointStore::ResultState::stMaxIter:  return "rM";
    }
    return "r?";
}

static const char *state_name(LaguerrePointStore::State s)
{
    switch (s) {
        case LaguerrePointStore::State::stUnknown:  return "sU";
        case LaguerrePointStore::State::stWorking:  return "sW";
        case LaguerrePointStore::State::stResolved: return "sR";
        case LaguerrePointStore::State::stFail:     return "sF";
    }
    return "s?";
}

// ── PPM writer ─────────────────────────────────────────────────────────────────
static bool write_ppm(const char *path, const QImage &img)
{
    FILE *f = std::fopen(path, "wb");
    if (!f) { std::perror(path); return false; }
    int W = img.width(), H = img.height();
    std::fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            unsigned argb = static_cast<unsigned>(img.pixel(x, y));
            unsigned char rgb[3] = {
                static_cast<unsigned char>((argb >> 16) & 0xff),
                static_cast<unsigned char>((argb >>  8) & 0xff),
                static_cast<unsigned char>( argb        & 0xff),
            };
            std::fwrite(rgb, 1, 3, f);
        }
    }
    std::fclose(f);
    return true;
}

// ── Wait for model threads ────────────────────────────────────────────────────
template <typename MODEL>
static bool wait_for_completion(MODEL &model, MainQThread &mainThread)
{
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(3000);
    while (model.getThreadsWorking() > 0) {
        mainThread.processEvents();
        if (std::chrono::steady_clock::now() > deadline) {
            std::fprintf(stderr, "Timeout: %d threads still working\n",
                         model.getThreadsWorking());
            return false;
        }
        std::this_thread::yield();
    }
    return true;
}

// ── Per-fractal compute ───────────────────────────────────────────────────────
// The view centre comes from g_cx / g_cy; in --point mode paintOrbit is skipped so
// only the worker thread evaluates the point.
static bool compute_mandelbrot(MandelModelInspector &model, MainQThread &mainThread,
                                int W, int H, double scale)
{
    auto t0 = Clock::now();
    select_precision(model, mainThread);
    // IMPORTANT: setImageSize before setView_str (see AGENTS.md).
    model.setImageSize(W, H);
    model.setView_str(g_cx, g_cy, scale);
    QImage image;
    ShareableImageWrapper wrapper(&image);
    if (!g_point)
        model.paintOrbit(wrapper, W/2, H/2); // paintOrbit must work
    // Use ExterAngle during evaluation so exterior-angle data is populated.
    model.setselectedPaintStyle(MandelModel::paintStyleExterAngle);
    double setup_ms = ms_since(t0);
    t0 = Clock::now();
    model.startRunning();
    bool ok = wait_for_completion(model, mainThread);
    report_timing(setup_ms, ms_since(t0));
    return ok;
}

// Julia/Laguerre parameters go through MandelModel::makeViewInfo, as in the GUI.
static ShareableViewInfo make_view_info(MandelModel &mandel, MainQThread &mainThread,
                                        int W, int H, double scale, const char *caption)
{
    select_precision(mandel, mainThread);
    mandel.setImageSize(W, H);
    mandel.setView_str(g_cx, g_cy, scale);
    QVariantMap params;
    params["viewRe"] = g_cx;
    params["viewIm"] = g_cy;
    params["viewZoom"] = scale;
    params["cRe"] = g_c_re;
    params["cIm"] = g_c_im;
    params["period"] = g_period;
    params["caption"] = caption;
    return mandel.makeViewInfo(params);
}

static bool compute_julia(JuliaModelInspector &model, MainQThread &mainThread,
                           int W, int H, double scale)
{
    auto t0 = Clock::now();
    select_precision(model, mainThread);
    model.setImageSize(W, H);
    model.setView_str(g_cx, g_cy, scale);
    QImage image;
    ShareableImageWrapper wrapper(&image);
    if (!g_point)
        model.paintOrbit(wrapper, W/2, H/2); // paintOrbit must work
    MandelModel mandel;
    ShareableViewInfo viewInfo=make_view_info(mandel, mainThread, W, H, scale, "Main view");
    model.setParams(viewInfo);
    model.setselectedPaintStyle(JuliaModel::paintStyleExterAngle);
    double setup_ms = ms_since(t0);
    t0 = Clock::now();
    model.startRunning();
    bool ok = wait_for_completion(model, mainThread);
    report_timing(setup_ms, ms_since(t0));
    return ok;
}

static bool compute_laguerre(LaguerreModelInspector &model, MainQThread &mainThread,
                              int W, int H, double scale)
{
    auto t0 = Clock::now();
    select_precision(model, mainThread);
    model.setImageSize(W, H);
    model.setView_str(g_cx, g_cy, scale);
    QImage image;
    ShareableImageWrapper wrapper(&image);
    MandelModel mandel;
    ShareableViewInfo viewInfo=make_view_info(mandel, mainThread, W, H, scale, "Evaluate");
    model.setParams(viewInfo);
    if (!g_point)
        model.paintOrbit(wrapper, W/2, H/2); // paintOrbit must work
    double setup_ms = ms_since(t0);
    t0 = Clock::now();
    model.startRunning();
    bool ok = wait_for_completion(model, mainThread);
    report_timing(setup_ms, ms_since(t0));
    return ok;
}

// ── Point mode ────────────────────────────────────────────────────────────────
static const char *qs(const QString &s)
{
    static std::string buf[8];
    static int next = 0;
    std::string &b = buf[next++ % 8];
    b = s.toStdString();
    return b.c_str();
}

static void print_point_man(MandelModelInspector &model)
{
    const MandelPointStore &p = model.store()[0];
    MandelPoint<MandelMath::number_any> &n = model.point(0);
    std::printf("c          = %s\n", qs(model.center()));
    std::printf("rstate=%s period=%d iter=%d newton_iter=%d near0iter_1=%d nearziter_0=%d surehand=%d\n",
                rstate_name(p.rstate), p.period, p.iter, p.newton_iter, p.near0iter_1, p.nearziter_0, p.surehand);
    std::printf("exterior_hits=%.17g exterior_avoids=%.17g interior_hits=%.17g\n",
                p.exterior_hits, p.exterior_avoids, p.interior.hits);
    std::printf("f          = %s\n", qs(n.f.toString()));
    std::printf("fc_c       = %s\n", qs(n.fc_c.toString()));
    std::printf("fz_r       = %s\n", qs(n.fz_r.toString()));
    std::printf("root       = %s\n", qs(n.root.toString()));
    std::printf("extangle   = %s\n", qs(n.extangle.toString()));
}

static void print_point_jul(JuliaModelInspector &model)
{
    const JuliaPointStore &p = model.store()[0];
    JuliaPoint<MandelMath::number_any> &n = model.point(0);
    std::printf("c          = %s\n", qs(model.c()));
    std::printf("z          = %s\n", qs(model.center()));
    std::printf("rstate=%s iter=%d near0iter_1=%d nearziter_0=%d nearmriter=%d\n",
                rstate_name(p.rstate), p.iter, p.near0iter_1, p.nearziter_0, p.nearmriter);
    std::printf("exterior none=%.17g/%.17g blend=%.17g/%.17g (mixed/empty)\n",
                p.exterior.none_mixed, p.exterior.none_empty, p.exterior.blend_mixed, p.exterior.blend_empty);
    std::printf("interior mixed=%.17g full=%.17g phi=%.17g%+.17gi first_under_1=%d\n",
                p.interior.mixed, p.interior.full, p.interior.phi_re, p.interior.phi_im, p.interior.first_under_1);
    std::printf("f          = %s\n", qs(n.f.toString()));
    std::printf("fz_z       = %s\n", qs(n.fz_z.toString()));
    std::printf("extangle   = %s\n", qs(n.extangle.toString()));
}

static void print_point_lag(LaguerreModelInspector &model)
{
    const LaguerrePointStore &p = model.store()[0];
    LaguerrePoint<MandelMath::number_any> &n = model.point(0);
    std::printf("c          = %s\n", qs(model.c()));
    std::printf("z          = %s\n", qs(model.center()));
    std::printf("state=%s iter=%d firstM=%.17g firstStep=%.17g%+.17gi naiveChoice=%d\n",
                state_name(p.state.load()), p.iter, p.firstM, p.firstStep_re, p.firstStep_im, (int)p.naiveChoice);
    std::printf("r          = %s\n", qs(n.r.toString()));
    std::printf("fz_r       = %s\n", qs(n.fz_r.toString()));
    std::printf("nth_fz     = %s\n", qs(n.nth_fz.toString()));
}

template <typename MODEL, typename COMPUTE, typename PRINT>
static void run_point(MODEL &model, MainQThread &mainThread, COMPUTE compute, PRINT print)
{
    // Fine enough that Position::setView does not round the point to the pixel grid.
    const double scale = std::ldexp(1.0, -300);
    std::printf("precision  = %s\n", PRECISION_NAMES[g_precision]);
    if (!compute(model, mainThread, 1, 1, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");
    print(model);
}

// ── Generic VT100 render helper ───────────────────────────────────────────────
static void print_vt100(const QImage &img, const char *style_name)
{
    int W = img.width(), H = img.height();
    std::printf("VT100 render (paintStyle=%s):\n", style_name);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            unsigned argb = static_cast<unsigned>(img.pixel(x, y));
            std::printf("\033[48;2;%d;%d;%dm  \033[0m",
                        (argb >> 16) & 0xff,
                        (argb >>  8) & 0xff,
                         argb        & 0xff);
        }
        std::putchar('\n');
    }
}

// ── Generic images mode ───────────────────────────────────────────────────────
// Works for any model whose paintStyle enum values are stored in PaintStyleEntry.
// prefix: "Man_", "Jul_", or "Lag_"
template<typename MODEL>
static void run_images_generic(MODEL &model,
                                const char *title,
                                double cx, double cy, double scale,
                                int W, int H,
                                const char *out_dir,
                                const char *prefix,
                                const PaintStyleEntry *styles, int n_styles,
                                int cls_idx)
{
    std::printf("\n=== %s  %dx%d  centre=(%g, %g)  scale=%.5g ===\n",
                title, W, H, cx, cy, scale);

    QImage img(W, H, QImage::Format::Format_ARGB32);
    ShareableImageWrapper wrapper(&img);
    std::string stem = sanitize_title(title);

    for (int si = 0; si < n_styles; ++si) {
        model.setselectedPaintStyle(
            static_cast<typename MODEL::paintStyle>(styles[si].value));
        int rc = model.writeToImage(wrapper);
        if (rc != 0) {
            std::fprintf(stderr, "  writeToImage(ps=%s) returned %d\n",
                         styles[si].name, rc);
            continue;
        }

        char ppm_path[1024];
        std::snprintf(ppm_path, sizeof(ppm_path),
                      "%s/%s%s_%s.ppm",
                      out_dir, prefix, stem.c_str(), styles[si].name);
        if (write_ppm(ppm_path, img))
            std::printf("  Wrote %s\n", ppm_path);

        if (si == cls_idx)
            print_vt100(img, styles[si].name);
    }
}

// ── Text mode — Mandelbrot ────────────────────────────────────────────────────
static void run_text_mode_man(MandelModelInspector &model, MainQThread &mainThread,
                               const char *title,
                               double cx, double cy, double scale,
                               int W, int H,
                               const char *out_dir)
{
    if (!compute_mandelbrot(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");

    const MandelPointStore *ps = model.store();
    QImage img(W, H, QImage::Format::Format_ARGB32);
    ShareableImageWrapper wrapper(&img);

    std::vector<std::vector<unsigned>> colors(N_MAN_STYLES,
                                              std::vector<unsigned>(W * H, 0));
    for (int si = 0; si < N_MAN_STYLES; ++si) {
        model.setselectedPaintStyle(
            static_cast<MandelModel::paintStyle>(MAN_PAINT_STYLES[si].value));
        int rc = model.writeToImage(wrapper);
        if (rc != 0) {
            std::fprintf(stderr, "  writeToImage(ps=%s) returned %d\n",
                         MAN_PAINT_STYLES[si].name, rc);
        } else {
            for (int y = 0; y < H; ++y)
                for (int x = 0; x < W; ++x)
                    colors[si][y * W + x] = static_cast<unsigned>(img.pixel(x, y));
        }
    }

    std::printf("\n=== %s  %dx%d  centre=(%.15g, %.15g)  scale=%.5g ===\n",
                title, W, H, cx, cy, scale);
    // Use Cls colour buffer for VT100 render.
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            unsigned argb = colors[MAN_CLS_IDX][y * W + x];
            std::printf("\033[48;2;%d;%d;%dm  \033[0m",
                        (argb >> 16) & 0xff, (argb >> 8) & 0xff, argb & 0xff);
        }
        std::putchar('\n');
    }

    char json_path[1024];
    std::snprintf(json_path, sizeof(json_path),
                  "%s/%s%s_text.json", out_dir, file_prefix("Man_").c_str(), sanitize_title(title).c_str());

    FILE *f = std::fopen(json_path, "w");
    if (!f) { std::perror(json_path); return; }

    std::fprintf(f, "{\n");
    std::fprintf(f, "  \"fractal\": \"man\",\n");
    std::fprintf(f, "  \"title\": \"%s\",\n", title);
    std::fprintf(f, "  \"cx\": %.17g,\n", cx);
    std::fprintf(f, "  \"cy\": %.17g,\n", cy);
    std::fprintf(f, "  \"scale\": %.17g,\n", scale);
    std::fprintf(f, "  \"width\": %d,\n", W);
    std::fprintf(f, "  \"height\": %d,\n", H);
    std::fprintf(f, "  \"paint_styles\": [");
    for (int si = 0; si < N_MAN_STYLES; ++si)
        std::fprintf(f, "%s\"%s\"", si ? ", " : "", MAN_PAINT_STYLES[si].name);
    std::fprintf(f, "],\n");

    std::fprintf(f, "  \"grid\": [\n");
    for (int y = 0; y < H; ++y) {
        std::fprintf(f, "    [");
        for (int x = 0; x < W; ++x) {
            const MandelPointStore &p = ps[y * W + x];
            std::fprintf(f,
                "\n      {\"y\":%d,\"x\":%d,"
                "\"rstate\":\"%s\","
                "\"period\":%d,"
                "\"iter\":%d,"
                "\"exterior_avoids\":%.6g,"
                "\"interior_hits\":%.6g,"
                "\"colors\":[",
                y, x,
                rstate_name(p.rstate),
                p.period, p.iter,
                p.exterior_avoids, p.interior.hits);
            for (int si = 0; si < N_MAN_STYLES; ++si)
                std::fprintf(f, "%s\"%08X\"", si ? "," : "",
                             colors[si][y * W + x]);
            std::fprintf(f, "]}%s", (x < W - 1) ? "," : "");
        }
        std::fprintf(f, "\n    ]%s\n", (y < H - 1) ? "," : "");
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    std::printf("  Wrote %s\n", json_path);
}

// ── Images mode — Mandelbrot ──────────────────────────────────────────────────
static void run_images_mode_man(MandelModelInspector &model, MainQThread &mainThread,
                                 const char *title,
                                 double cx, double cy, double scale,
                                 int W, int H,
                                 const char *out_dir)
{
    if (!compute_mandelbrot(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");
    run_images_generic(model, title, cx, cy, scale, W, H, out_dir,
                       file_prefix("Man_").c_str(), MAN_PAINT_STYLES, N_MAN_STYLES, MAN_CLS_IDX);
}

// ── Text mode — Julia ─────────────────────────────────────────────────────────
static void run_text_mode_jul(JuliaModelInspector &model, MainQThread &mainThread,
                               const char *title,
                               double cx, double cy, double scale,
                               int W, int H,
                               const char *out_dir)
{
    if (!compute_julia(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");

    const JuliaPointStore *ps = model.store();
    QImage img(W, H, QImage::Format::Format_ARGB32);
    ShareableImageWrapper wrapper(&img);

    std::vector<std::vector<unsigned>> colors(N_JUL_STYLES,
                                              std::vector<unsigned>(W * H, 0));
    for (int si = 0; si < N_JUL_STYLES; ++si) {
        model.setselectedPaintStyle(
            static_cast<JuliaModel::paintStyle>(JUL_PAINT_STYLES[si].value));
        int rc = model.writeToImage(wrapper);
        if (rc != 0) {
            std::fprintf(stderr, "  writeToImage(ps=%s) returned %d\n",
                         JUL_PAINT_STYLES[si].name, rc);
        } else {
            for (int y = 0; y < H; ++y)
                for (int x = 0; x < W; ++x)
                    colors[si][y * W + x] = static_cast<unsigned>(img.pixel(x, y));
        }
    }

    std::printf("\n=== %s  %dx%d  centre=(%.15g, %.15g)  scale=%.5g ===\n",
                title, W, H, cx, cy, scale);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            unsigned argb = colors[JUL_CLS_IDX][y * W + x];
            std::printf("\033[48;2;%d;%d;%dm  \033[0m",
                        (argb >> 16) & 0xff, (argb >> 8) & 0xff, argb & 0xff);
        }
        std::putchar('\n');
    }

    char json_path[1024];
    std::snprintf(json_path, sizeof(json_path),
                  "%s/%s%s_text.json", out_dir, file_prefix("Jul_").c_str(), sanitize_title(title).c_str());

    FILE *f = std::fopen(json_path, "w");
    if (!f) { std::perror(json_path); return; }

    std::fprintf(f, "{\n");
    std::fprintf(f, "  \"fractal\": \"jul\",\n");
    std::fprintf(f, "  \"title\": \"%s\",\n", title);
    std::fprintf(f, "  \"cx\": %.17g,\n", cx);
    std::fprintf(f, "  \"cy\": %.17g,\n", cy);
    std::fprintf(f, "  \"scale\": %.17g,\n", scale);
    std::fprintf(f, "  \"width\": %d,\n", W);
    std::fprintf(f, "  \"height\": %d,\n", H);
    std::fprintf(f, "  \"paint_styles\": [");
    for (int si = 0; si < N_JUL_STYLES; ++si)
        std::fprintf(f, "%s\"%s\"", si ? ", " : "", JUL_PAINT_STYLES[si].name);
    std::fprintf(f, "],\n");

    std::fprintf(f, "  \"grid\": [\n");
    for (int y = 0; y < H; ++y) {
        std::fprintf(f, "    [");
        for (int x = 0; x < W; ++x) {
            const JuliaPointStore &p = ps[y * W + x];
            std::fprintf(f,
                "\n      {\"y\":%d,\"x\":%d,"
                "\"rstate\":\"%s\","
                "\"iter\":%d,"
                "\"colors\":[",
                y, x,
                rstate_name(p.rstate),
                p.iter);
            for (int si = 0; si < N_JUL_STYLES; ++si)
                std::fprintf(f, "%s\"%08X\"", si ? "," : "",
                             colors[si][y * W + x]);
            std::fprintf(f, "]}%s", (x < W - 1) ? "," : "");
        }
        std::fprintf(f, "\n    ]%s\n", (y < H - 1) ? "," : "");
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    std::printf("  Wrote %s\n", json_path);
}

// ── Images mode — Julia ───────────────────────────────────────────────────────
static void run_images_mode_jul(JuliaModelInspector &model, MainQThread &mainThread,
                                 const char *title,
                                 double cx, double cy, double scale,
                                 int W, int H,
                                 const char *out_dir)
{
    if (!compute_julia(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");
    run_images_generic(model, title, cx, cy, scale, W, H, out_dir,
                       file_prefix("Jul_").c_str(), JUL_PAINT_STYLES, N_JUL_STYLES, JUL_CLS_IDX);
}

// ── Text mode — Laguerre ──────────────────────────────────────────────────────
static void run_text_mode_lag(LaguerreModelInspector &model, MainQThread &mainThread,
                               const char *title,
                               double cx, double cy, double scale,
                               int W, int H,
                               const char *out_dir)
{
    if (!compute_laguerre(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");

    const LaguerrePointStore *ps = model.store();
    QImage img(W, H, QImage::Format::Format_ARGB32);
    ShareableImageWrapper wrapper(&img);

    std::vector<std::vector<unsigned>> colors(N_LAG_STYLES,
                                              std::vector<unsigned>(W * H, 0));
    for (int si = 0; si < N_LAG_STYLES; ++si) {
        model.setselectedPaintStyle(
            static_cast<LaguerreModel::paintStyle>(LAG_PAINT_STYLES[si].value));
        int rc = model.writeToImage(wrapper);
        if (rc != 0) {
            std::fprintf(stderr, "  writeToImage(ps=%s) returned %d\n",
                         LAG_PAINT_STYLES[si].name, rc);
        } else {
            for (int y = 0; y < H; ++y)
                for (int x = 0; x < W; ++x)
                    colors[si][y * W + x] = static_cast<unsigned>(img.pixel(x, y));
        }
    }

    std::printf("\n=== %s  %dx%d  centre=(%.15g, %.15g)  scale=%.5g ===\n",
                title, W, H, cx, cy, scale);
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            unsigned argb = colors[LAG_CLS_IDX][y * W + x];
            std::printf("\033[48;2;%d;%d;%dm  \033[0m",
                        (argb >> 16) & 0xff, (argb >> 8) & 0xff, argb & 0xff);
        }
        std::putchar('\n');
    }

    char json_path[1024];
    std::snprintf(json_path, sizeof(json_path),
                  "%s/%s%s_text.json", out_dir, file_prefix("Lag_").c_str(), sanitize_title(title).c_str());

    FILE *f = std::fopen(json_path, "w");
    if (!f) { std::perror(json_path); return; }

    std::fprintf(f, "{\n");
    std::fprintf(f, "  \"fractal\": \"lag\",\n");
    std::fprintf(f, "  \"title\": \"%s\",\n", title);
    std::fprintf(f, "  \"cx\": %.17g,\n", cx);
    std::fprintf(f, "  \"cy\": %.17g,\n", cy);
    std::fprintf(f, "  \"scale\": %.17g,\n", scale);
    std::fprintf(f, "  \"width\": %d,\n", W);
    std::fprintf(f, "  \"height\": %d,\n", H);
    std::fprintf(f, "  \"paint_styles\": [");
    for (int si = 0; si < N_LAG_STYLES; ++si)
        std::fprintf(f, "%s\"%s\"", si ? ", " : "", LAG_PAINT_STYLES[si].name);
    std::fprintf(f, "],\n");

    std::fprintf(f, "  \"grid\": [\n");
    for (int y = 0; y < H; ++y) {
        std::fprintf(f, "    [");
        for (int x = 0; x < W; ++x) {
            const LaguerrePointStore &p = ps[y * W + x];
            std::fprintf(f,
                "\n      {\"y\":%d,\"x\":%d,"
                "\"state\":\"%s\","
                "\"iter\":%d,"
                "\"firstM\":%.6g,"
                "\"colors\":[",
                y, x,
                state_name(p.state.load()),
                p.iter,
                p.firstM);
            for (int si = 0; si < N_LAG_STYLES; ++si)
                std::fprintf(f, "%s\"%08X\"", si ? "," : "",
                             colors[si][y * W + x]);
            std::fprintf(f, "]}%s", (x < W - 1) ? "," : "");
        }
        std::fprintf(f, "\n    ]%s\n", (y < H - 1) ? "," : "");
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    std::printf("  Wrote %s\n", json_path);
}

// ── Images mode — Laguerre ────────────────────────────────────────────────────
static void run_images_mode_lag(LaguerreModelInspector &model, MainQThread &mainThread,
                                 const char *title,
                                 double cx, double cy, double scale,
                                 int W, int H,
                                 const char *out_dir)
{
    if (!compute_laguerre(model, mainThread, W, H, scale))
        std::fprintf(stderr, "Warning: computation timed out.\n");
    run_images_generic(model, title, cx, cy, scale, W, H, out_dir,
                       file_prefix("Lag_").c_str(), LAG_PAINT_STYLES, N_LAG_STYLES, LAG_CLS_IDX);
}

// ── Argument parsing ──────────────────────────────────────────────────────────

enum class Mode    { Text, Images };
enum class Fractal { Man, Jul, Lag };

struct Args {
    const char *cx      = "-0.5";
    const char *cy      = "0";
    double      scale   =  0.125;
    int         width   =  40;
    int         height  =  20;
    const char *title   = "Wide view";
    Mode        mode    = Mode::Text;
    Fractal     fractal = Fractal::Man;
    bool        bless   = false;
    bool        point   = false;
    const char *c_re    = nullptr;   // nullptr/0 = per-fractal default
    const char *c_im    = nullptr;
    int         period  = 0;
    int         precision = 0;
};

static void print_usage(const char *prog)
{
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "  --fractal man|jul|lag  Fractal type       (default man)\n"
        "  --cx <num>       Centre real part         (default -0.5)\n"
        "  --cy <num>       Centre imag part         (default 0.0)\n"
        "  --scale <float>  Scale                    (default 0.125)\n"
        "  --width <int>    Width in pixels          (default 40)\n"
        "  --height <int>   Height in pixels         (default 20)\n"
        "  --title <str>    Header label             (default 'Wide view')\n"
        "  --c-re <num>     Julia/Laguerre c, real   (default 0.3 / -0.7067658320364081)\n"
        "  --c-im <num>     Julia/Laguerre c, imag   (default 0.5 / 0.3540325897413918)\n"
        "  --period <int>   Julia/Laguerre period    (default 4 / 1632)\n"
        "  --precision double|float128|dd|qd|real642|ddheap (default double)\n"
        "  --mode text|images                        (default text)\n"
        "  --bless          Write to ./blessed/ instead of ./current/\n"
        "  --point          Only evaluate the point (--cx, --cy) and print the result\n"
        "  <num> is a decimal string, parsed in the selected precision\n",
        prog);
}

static const char *number_arg(const char *prog, const char *opt, const char *val)
{
    MandelMath::number<MandelMath::real642> check; // parses as both double and __float128
    if (!check.fromString(val)) {
        std::fprintf(stderr, "Not a number for %s: %s\n", opt, val);
        print_usage(prog);
        std::exit(1);
    }
    return val;
}

static Args parse_args(int argc, char *argv[])
{
    Args a;
    for (int i = 1; i < argc; ++i) {
        if      (std::strcmp(argv[i], "--bless") == 0) { a.bless = true; }
        else if (std::strcmp(argv[i], "--point") == 0) { a.point = true; }
        else if (i + 1 < argc) {
            if      (std::strcmp(argv[i], "--cx")      == 0) { a.cx     = number_arg(argv[0], argv[i], argv[i+1]); ++i; }
            else if (std::strcmp(argv[i], "--cy")      == 0) { a.cy     = number_arg(argv[0], argv[i], argv[i+1]); ++i; }
            else if (std::strcmp(argv[i], "--scale")   == 0) { a.scale  = std::atof(argv[++i]); }
            else if (std::strcmp(argv[i], "--width")   == 0) { a.width  = std::atoi(argv[++i]); }
            else if (std::strcmp(argv[i], "--height")  == 0) { a.height = std::atoi(argv[++i]); }
            else if (std::strcmp(argv[i], "--title")   == 0) { a.title  = argv[++i]; }
            else if (std::strcmp(argv[i], "--c-re")    == 0) { a.c_re   = number_arg(argv[0], argv[i], argv[i+1]); ++i; }
            else if (std::strcmp(argv[i], "--c-im")    == 0) { a.c_im   = number_arg(argv[0], argv[i], argv[i+1]); ++i; }
            else if (std::strcmp(argv[i], "--period")  == 0) { a.period = std::atoi(argv[++i]); }
            else if (std::strcmp(argv[i], "--precision") == 0) {
                ++i;
                a.precision = -1;
                for (int p = 0; p < N_PRECISIONS; ++p)
                    if (std::strcmp(argv[i], PRECISION_NAMES[p]) == 0)
                        a.precision = p;
                if (a.precision < 0) {
                    std::fprintf(stderr, "Unknown precision: %s\n", argv[i]);
                    print_usage(argv[0]);
                    std::exit(1);
                }
            }
            else if (std::strcmp(argv[i], "--mode")    == 0) {
                ++i;
                if      (std::strcmp(argv[i], "text")   == 0) { a.mode = Mode::Text; }
                else if (std::strcmp(argv[i], "images") == 0) { a.mode = Mode::Images; }
                else {
                    std::fprintf(stderr, "Unknown mode: %s  (expected text or images)\n", argv[i]);
                    print_usage(argv[0]);
                    std::exit(1);
                }
            }
            else if (std::strcmp(argv[i], "--fractal") == 0) {
                ++i;
                if      (std::strncmp(argv[i], "man", 3) == 0) { a.fractal = Fractal::Man; }
                else if (std::strncmp(argv[i], "jul", 3) == 0) { a.fractal = Fractal::Jul; }
                else if (std::strncmp(argv[i], "lag", 3) == 0) { a.fractal = Fractal::Lag; }
                else {
                    std::fprintf(stderr, "Unknown fractal: %s  (expected man, jul, or lag)\n", argv[i]);
                    print_usage(argv[0]);
                    std::exit(1);
                }
            }
            else {
                std::fprintf(stderr, "Unknown option: %s\n", argv[i]);
                print_usage(argv[0]);
                std::exit(1);
            }
        } else {
            std::fprintf(stderr, "Missing value for: %s\n", argv[i]);
            print_usage(argv[0]);
            std::exit(1);
        }
    }
    return a;
}

// ── main ──────────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    Args a = parse_args(argc, argv);

    const bool lag = a.fractal == Fractal::Lag;
    g_precision = a.precision;
    g_cx     = a.cx;
    g_cy     = a.cy;
    g_c_re   = a.c_re ? a.c_re : lag ? "-0.7067658320364081" : "0.3";
    g_c_im   = a.c_im ? a.c_im : lag ?  "0.3540325897413918" : "0.5";
    g_period = a.period > 0 ? a.period : lag ? 1632 : 4;
    g_point  = a.point;
    // double copies only for the printed headers and JSON
    const double cx = std::atof(a.cx);
    const double cy = std::atof(a.cy);

    MainQThread mainThread;

    if (g_point) {
        if (a.fractal == Fractal::Man) {
            MandelModelInspector model;
            run_point(model, mainThread, compute_mandelbrot, print_point_man);
        } else if (a.fractal == Fractal::Jul) {
            JuliaModelInspector model;
            run_point(model, mainThread, compute_julia, print_point_jul);
        } else {
            LaguerreModelInspector model;
            run_point(model, mainThread, compute_laguerre, print_point_lag);
        }
        return 0;
    }

    const char *out_dir = a.bless ? "blessed" : "current";
    mkdir(out_dir, 0755);   // ignore error if already exists

    if (a.fractal == Fractal::Man) {
        MandelModelInspector model;
        if (a.mode == Mode::Text)
            run_text_mode_man(model, mainThread,
                              a.title, cx, cy, a.scale, a.width, a.height,
                              out_dir);
        else
            run_images_mode_man(model, mainThread,
                                a.title, cx, cy, a.scale, a.width, a.height,
                                out_dir);
    } else if (a.fractal == Fractal::Jul) {
        JuliaModelInspector model;
        if (a.mode == Mode::Text)
            run_text_mode_jul(model, mainThread,
                              a.title, cx, cy, a.scale, a.width, a.height,
                              out_dir);
        else
            run_images_mode_jul(model, mainThread,
                                a.title, cx, cy, a.scale, a.width, a.height,
                                out_dir);
    } else {  // Fractal::Lag
        LaguerreModelInspector model;
        if (a.mode == Mode::Text)
            run_text_mode_lag(model, mainThread,
                              a.title, cx, cy, a.scale, a.width, a.height,
                              out_dir);
        else
            run_images_mode_lag(model, mainThread,
                                a.title, cx, cy, a.scale, a.width, a.height,
                                out_dir);
    }

    std::printf("\nDone.\n");
    return 0;
}

