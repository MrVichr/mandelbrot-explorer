#ifndef MANDELMODEL_H
#define MANDELMODEL_H

#include <QObject>
#include <QImage>
#include <QReadWriteLock>

#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class MandelModel: public QObject
{
  Q_OBJECT
protected:
  enum IndexIntoWorker
  {
    iiw_reImToPixel=10
  };
public:
  MandelModel();
  ~MandelModel();
  Q_INVOKABLE void startRunning();
  /*template <class WORKER_MULTI>
  void transformStore(WORKER_MULTI *old_worker, WORKER_MULTI *old_sworker, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex<WORKER_MULTI> *old_c,
                      WORKER_MULTI *new_worker, WORKER_MULTI *new_sworker, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex<WORKER_MULTI> *new_c,
                      int inlog, int new_step_log);*/
  void transformStore(void *old_points, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex<MandelMath::number_any> *old_c,
                      void *new_points, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex<MandelMath::number_any> *new_c,
                      int inlog, int new_step_log);
  Q_INVOKABLE void setView_double(double c_re, double c_im, double scale);
  void setView(MandelMath::complex<MandelMath::number_any> const &c, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  Q_INVOKABLE void pause(bool pause);
  void startNewEpoch();
  Q_INVOKABLE int writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, MandelMath::complex<MandelMath::number_any> const &z, MandelMath::number<MandelMath::number_any> *tmp);
  Q_INVOKABLE void paintOrbit(ShareableImageWrapper image, int x, int y);
  Q_INVOKABLE QString pixelXtoRE_str(int x);
  Q_INVOKABLE QString pixelYtoIM_str(int y);
  Q_INVOKABLE QString getTimes();
  Q_INVOKABLE QString getTextXY();
  Q_INVOKABLE QString getTextInfoGen();
  Q_INVOKABLE QString getTextInfoSpec();
  ShareableViewInfo getViewInfo();
  Q_PROPERTY(ShareableViewInfo viewInfo READ getViewInfo CONSTANT)// WRITE setViewInfo NOTIFY viewInfoChanged)
  Q_INVOKABLE ShareableViewInfo makeViewInfo(const QVariantMap &params);

  enum paintStyle
  {
    paintStyleKind=0,
    paintStyleCls=1,
    paintStyleExter=2,
    paintStyleExterAngle=3,
    paintStyleInter=4,
    paintStyleNear=5,
    paintStyleFZ=6,
    paintStyleFC=7
  };
  Q_ENUM(paintStyle);
  paintStyle _selectedPaintStyle;
  Q_PROPERTY(paintStyle selectedPaintStyle READ getselectedPaintStyle WRITE setselectedPaintStyle NOTIFY selectedPaintStyleChanged)
  paintStyle getselectedPaintStyle() { return _selectedPaintStyle; }
  void setselectedPaintStyle(paintStyle ps) { _selectedPaintStyle=ps; }
  int _threadsWorking;
  Q_PROPERTY(int threadsWorking READ getThreadsWorking CONSTANT)
  int getThreadsWorking() { return _threadsWorking; }
  Q_PROPERTY(int threadsMax READ getThreadCount CONSTANT)
  int getThreadCount() { return precisionRecord->threadCount; }

  enum precision
  {
    precisionDouble=0,
#if !ONLY_DOUBLE_WORKER
    precisionFloat128=1,
    precisionDDouble=2,
    precisionQDouble=3,
    precisionReal642=4
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
  void setextAngleZoom(int zoom) { _extAngleZoom=zoom; }

  QVector<int> periodToIndexCache;
  int periodToIndex(int period);
protected:
  QReadWriteLock threading_mutex;
  template <typename BASE>
  int giveWorkThreaded(MandelEvaluator<BASE> *me);
  template <typename BASE>
  int doneWorkThreaded(MandelEvaluator<BASE> *me, bool giveWork);
public slots:
  void doneWorkInThread(MandelEvaluatorThread *me);
  void selectedPrecisionChanged();
signals:
  void selectedPaintStyleChanged();
  void selectedPrecisionChange();
  void extAngleZoomChange();
  void triggerComputeThreaded(int epoch);
protected:
  //MandelMath::worker_multi::Allocator<MandelMath::worker_multi> *storeAllocator;
  //MandelMath::worker_multi *storeWorker; //pointStore
  MandelPointStore *pointStore;

  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  int nextGivenPointIndex;
  int effortBonus;
  //constexpr static int MAX_EFFORT=17;//131072 iters;
  static constexpr int MAX_EFFORT=22;//
  QElapsedTimer timerWriteToImage;

  struct Position
  {
    static constexpr int LEN=2;
    MandelMath::complex<MandelMath::number_any> center;
    int step_log;
    double step_size; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position(MandelMath::complex<MandelMath::number_any>::Scratchpad *spad, const Position *source);
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
    //MandelMath::worker_multi::Allocator pointAllocator;
    //MandelPointStore pointDataStore;
    //-> evaluator.currentData MandelPoint pointData;
    struct Bulb
    {
      bool valid_unused;
      MandelMath::complex<MandelMath::number_any> cb_unused;
      MandelMath::complex<MandelMath::number_any> rb_unused;
      MandelMath::complex<MandelMath::number_any> xc_unused;
      MandelMath::complex<MandelMath::number_any> baseZC_unused;
      MandelMath::complex<MandelMath::number_any> baseCC_unused;
      MandelMath::complex<MandelMath::number_any> baseFz;
      int foundMult_;
      bool is_card_unused;
      Bulb(MandelMath::complex<MandelMath::number_any>::Scratchpad *spad);
      ~Bulb();
      constexpr static int LEN=12;
    } bulb;
    Orbit(MandelMath::NumberType ntype);
    ~Orbit();
    constexpr static int LEN=Bulb::LEN;
  };
  struct PrecisionRecord
  {
    /*enum IndexIntoWorker
    {
      iiw_wtiPoint=5,
      iiw__LEN=0,
    };*/
    MandelMath::NumberType ntype;
    Orbit orbit;
    MandelPoint<MandelMath::number_any> wtiPoint;
    Position position;
    MandelMath::complex<MandelMath::number_any> lagu_c; //TODO: create struct Params and move there, like LaguerreModel
    MandelMath::complex<MandelMath::number_any> lagu_r; //same
    MandelMath::number<MandelMath::number_any> tmp_place;
    void *points; //array<double|float128|...>[width*height]
    int threadCount;
    std::variant<std::nullptr_t,
                 MandelEvaluator<double> **,
                 MandelEvaluator<__float128> **,
                 MandelEvaluator<MandelMath::dd_real> **,
                 MandelEvaluator<MandelMath::dq_real> **,
                 MandelEvaluator<MandelMath::real642> **> threads;
    PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, MandelModel *doneReceiver);
    ~PrecisionRecord();
    //constexpr static int LEN=ShareableViewInfo::LEN+MandelPoint_<MandelMath::number_a *>::LEN+Position::LEN+Orbit::LEN+5  +6;
      //setViewDouble=2 setView=2 updateCachedDepth=2 -> +6
  } *precisionRecord;
};

#endif // MANDELMODEL_H
