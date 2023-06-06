#ifndef JULIAMODEL_H
#define JULIAMODEL_H

#include <QObject>
#include <QImage>
#include <QReadWriteLock>

#include "atomic_update_rect.hpp"
#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class JuliaModel: public QObject
{
  Q_OBJECT
protected:
  /*enum IndexIntoWorker
  {
    iiw_reImToPixel=10
  };*/
public:
  JuliaModel();
  ~JuliaModel();
  Q_INVOKABLE void startRunning();
  void transformStore(void *old_points, JuliaPointStore *old_store, int old_width, int old_height, MandelMath::complex<MandelMath::number_any> const &old_c,
                      void *new_points, JuliaPointStore *new_store, int new_width, int new_height, MandelMath::complex<MandelMath::number_any> const &new_c,
                      int inlog, int new_step_log);
  void recomputeRoot(int max_effort);
  Q_INVOKABLE void setParams(ShareableViewInfo viewInfo);
  void setView(MandelMath::complex<MandelMath::number_any> const &c, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  Q_INVOKABLE void pause(bool pause);
  void startNewEpoch();
  //void giveWorkAll();
  Q_INVOKABLE int writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, MandelMath::complex<MandelMath::number_any> const &z, MandelMath::number<MandelMath::number_any> *tmp);
  Q_INVOKABLE void paintOrbit(ShareableImageWrapper image, int x, int y);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();
  Q_INVOKABLE QString getTextXY();
  Q_INVOKABLE QString getTextInfoGen();
  Q_INVOKABLE QString getTextInfoSpec();
  struct Params
  {
    static constexpr int LEN=4;
    int period; //TODO: 0 doesn't seem to be the best choice for "no interior", MAXINT would better
    MandelMath::complex<MandelMath::number_any> c;
    MandelMath::complex<MandelMath::number_any> root;
    Params(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Params *source);
  };


  enum paintStyle
  {
    paintStyleKind=0,
    paintStyleCls=1,
    paintStyleExter=2,
    paintStyleExterAngle=3,
    paintStyleInter=4,
    paintStyleNear0=5,
    paintStyleNearZ=6,
    paintStylePhi=7,
    paintStylePhi1=8,
    paintStylePhi2=9,
    paintStyleF1=10,
    paintStyleFirstUnder1=11,
    paintStyleNearMR=12,
  };
  Q_ENUM(paintStyle);
  paintStyle _selectedPaintStyle;
  Q_PROPERTY(paintStyle selectedPaintStyle READ getselectedPaintStyle WRITE setselectedPaintStyle NOTIFY selectedPaintStyleChanged)
  paintStyle getselectedPaintStyle() { return _selectedPaintStyle; }
  void setselectedPaintStyle(paintStyle ps) { _selectedPaintStyle=ps; invalidateMainImage(); }
  int _threadsWorking;
  Q_PROPERTY(int threadsWorking READ getThreadsWorking CONSTANT)
  int getThreadsWorking() { return _threadsWorking; }
  Q_PROPERTY(int threadsMax READ getThreadCount CONSTANT)
  int getThreadCount() { return precisionRecord->threadCount; }

  enum precision
  {
    EprecisionDouble=0,
#if !ONLY_DOUBLE_WORKER
    EprecisionFloat128=1,
    EprecisionDDouble=2,
    EprecisionQDouble=3,
    EprecisionReal642=4
#endif
  };
  Q_ENUM(precision);
  precision _selectedPrecision;
  Q_PROPERTY(precision selectedPrecision READ getselectedPrecision WRITE setselectedPrecision NOTIFY selectedPrecisionChange)
  precision getselectedPrecision() { return _selectedPrecision; }
  void setselectedPrecision(precision ps) { _selectedPrecision=ps; emit selectedPrecisionChange(); }

  int _extAngleZoom;
  Q_PROPERTY(int extAngleZoom READ getextAngleZoom WRITE setextAngleZoom NOTIFY extAngleZoomChange)
  int getextAngleZoom() { return _extAngleZoom; }
  void setextAngleZoom(int zoom) { _extAngleZoom=zoom; invalidateMainImage(); }

  enum ext_de_patch_algo
  { //the Internet says that enum names need to start with a capital letter to be accessible from QML :-?
    EedepaNone=0,
    EedepaAlways=1,
    EedepaConditional=2,
    EedepaBlend=3,
    EedepaVerifyH=4,
    EedepaVerifyA=5,
  };
  Q_ENUM(ext_de_patch_algo);
  ext_de_patch_algo _selectedEdePatchAlgo;
  Q_PROPERTY(ext_de_patch_algo selectedEdePatchAlgo READ getselectedEdePatchAlgo WRITE setselectedEdePatchAlgo NOTIFY selectedEdePatchAlgoChange)
  ext_de_patch_algo getselectedEdePatchAlgo() { return _selectedEdePatchAlgo; }
  void setselectedEdePatchAlgo(ext_de_patch_algo ps) { _selectedEdePatchAlgo=ps; invalidateMainImage(); emit selectedEdePatchAlgoChange(); }

  enum exterior_coloring
  {
    EextcolScaledWaves=0,
    EextcolRainbow=1,
  };
  Q_ENUM(exterior_coloring);
  exterior_coloring _selectedExteriorColoring;
  Q_PROPERTY(exterior_coloring selectedExteriorColoring READ getselectedExteriorColoring WRITE setselectedExteriorColoring NOTIFY selectedExteriorColoringChange)
  exterior_coloring getselectedExteriorColoring() { return _selectedExteriorColoring; }
  void setselectedExteriorColoring(exterior_coloring ps) { _selectedExteriorColoring=ps; invalidateMainImage(); emit selectedExteriorColoringChange(); }

  bool _orbit_frozen;
  Q_PROPERTY(bool orbit_frozen READ getorbit_frozen WRITE setorbit_frozen NOTIFY orbit_frozenChange)
  bool getorbit_frozen() { return _orbit_frozen; }
  void setorbit_frozen(bool of) { _orbit_frozen=of; emit orbit_frozenChange(); }

  QVector<int> periodToIndexCache;
  int periodToIndex(int period);
protected:
  QReadWriteLock threading_mutex;
  template <typename BASE>
  int giveWorkThreaded(MandelEvaluator<BASE> *me);
  template <typename BASE>
  int doneWorkThreaded(MandelEvaluator<BASE> *me, bool giveWork);

  AtomicUpdateRect image_dirty;
  void invalidateMainImage()
  {
    image_dirty.invalidate();
  }
public slots:
  void doneWorkInThread(MandelEvaluatorThread *me);
  void selectedPrecisionChanged();
  void selectedEdePatchAlgoChanged();
  void selectedExteriorColoringChanged();
signals:
  void selectedPaintStyleChanged();
  void selectedPrecisionChange();
  void selectedEdePatchAlgoChange();
  void selectedExteriorColoringChange();
  void orbit_frozenChange();
  void extAngleZoomChange();
  void triggerJuliaThreaded(int epoch, int juliaPeriod);
protected:
  //MandelMath::worker_multi::Allocator<MandelMath::worker_multi> *storeAllocator;
  //MandelMath::worker_multi *storeWorker; //pointStore
  JuliaPointStore *pointStore;
  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  int nextGivenPointIndex;
  int effortBonus;
  constexpr static int MAX_EFFORT=22;
  //constexpr static int MAX_EFFORT=13;
  QElapsedTimer timerWriteToImage;

  struct Position
  {
    static constexpr int LEN=2;
    MandelMath::complex<MandelMath::number_any> center;
    int step_log;
    double step_size; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Position *source);
    ~Position();
    //void assign(Position *src);
    void setView(MandelMath::complex<MandelMath::number_any> const &c, double scale);
    void move(int delta_x, int delta_y);
    void scale(int inlog, int center_x, int center_y);
    void updateCachedDepth();
    template <typename BASE>
    void pixelXtoRE(int x, MandelMath::number<BASE> *result);
    template <typename BASE>
    void pixelYtoIM(int y, MandelMath::number<BASE> *result);
  };
  struct Orbit
  {
    MandelEvaluator<MandelMath::number_any> evaluator;
    JuliaPointStore pointDataStore; //but we need it for long-term storage of results
    double first_mu_re, first_mu_im, first_mum_re, first_mum_im;
    Orbit(MandelMath::NumberType ntype, Orbit const *source);
    ~Orbit();
    constexpr static int LEN=JuliaPoint<MandelMath::number_any>::LEN;
  };
  struct PrecisionRecord
  {
    MandelMath::NumberType ntype;
    Orbit orbit;
    JuliaPoint<MandelMath::number_any> wtiPoint;
    Params params;
    Position position;
    //MandelMath::number<MandelMath::number_any> tmp_place;
    void *points; //array<double|float128|...>[width*height]
    int threadCount;
    std::variant<std::nullptr_t,
                 MandelEvaluator<double> **,
                 MandelEvaluator<__float128> **,
                 MandelEvaluator<MandelMath::dd_real> **,
                 MandelEvaluator<MandelMath::dq_real> **,
                 MandelEvaluator<MandelMath::real642> **> threads;
    PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, JuliaModel *doneReceiver);
    ~PrecisionRecord();
    //constexpr static int LEN=Params::LEN+JuliaPoint<MandelMath::number_a *>::LEN+Position::LEN+Orbit::LEN  +4   +1;
  } *precisionRecord;
};

#endif // JULIAMODEL_H
