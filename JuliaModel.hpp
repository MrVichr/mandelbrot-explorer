#ifndef JULIAMODEL_H
#define JULIAMODEL_H

#include <QObject>
#include <QImage>
#include <QReadWriteLock>

#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class JuliaModel: public QObject
{
  Q_OBJECT
protected:
  enum IndexIntoWorker
  {
    iiw_reImToPixel=10
  };
public:
  JuliaModel();
  ~JuliaModel();
  Q_INVOKABLE void startRunning();
  void transformStore(void *old_points, JuliaPointStore *old_store, int old_width, int old_height, const MandelMath::complex<MandelMath::number_a *> *old_c,
                      void *new_points, JuliaPointStore *new_store, int new_width, int new_height, const MandelMath::complex<MandelMath::number_a *> *new_c,
                      int inlog, int new_step_log);
  void recomputeRoot(int max_effort);
  Q_INVOKABLE void setParams(ShareableViewInfo viewInfo);
  void setView(const MandelMath::complex<MandelMath::number_a *> *c, double scale);
  Q_INVOKABLE void drag(double delta_x, double delta_y);
  Q_INVOKABLE void zoom(double x, double y, int inlog);
  Q_INVOKABLE void setImageSize(int width, int height);
  Q_INVOKABLE void pause(bool pause);
  void startNewEpoch();
  //void giveWorkAll();
  Q_INVOKABLE int writeToImage(ShareableImageWrapper img);
  void reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex<MandelMath::number_a *> *c, MandelMath::number<MandelMath::number_a *> *tmp);
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
    MandelMath::complex<MandelMath::number_a *> c;
    MandelMath::complex<MandelMath::number_a *> root;
    Params(MandelMath::NumberType ntype, const Params *source);
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
    MandelMath::complex<MandelMath::number_a *> center;
    int step_log;
    double step_size; //TODO: should use special methods on number to add, mul and div by 2^-step_log
    int cached_center_re_mod; //(center/step) mod 32768
    int cached_center_im_mod;
    Position(MandelMath::NumberType ntype, const Position *source);
    ~Position();
    //void assign(Position *src);
    void setView(const MandelMath::complex<MandelMath::number_a *> *c, double scale);
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
    MandelEvaluator<MandelMath::number_a *> evaluator;
    JuliaPointStore pointDataStore; //but we need it for long-term storage of results
    JuliaPoint<MandelMath::number_a *> pointData;
    double first_mu_re, first_mu_im, first_mum_re, first_mum_im;
    Orbit(MandelMath::NumberType ntype);
    ~Orbit();
    constexpr static int LEN=JuliaPoint<MandelMath::number_a *>::LEN;
  };
  struct PrecisionRecord
  {
    MandelMath::NumberType ntype;
    JuliaPoint<MandelMath::number_a *> wtiPoint;
    Params params;
    Position position;
    Orbit orbit;
    MandelMath::number<MandelMath::number_a *> tmp_place;
    void *points; //array<double|float128|...>[width*height]
    int threadCount;
    MandelEvaluator<MandelMath::number_a *> **threads;
    PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, JuliaModel *doneReceiver);
    ~PrecisionRecord();
    //constexpr static int LEN=Params::LEN+JuliaPoint<MandelMath::number_a *>::LEN+Position::LEN+Orbit::LEN  +4   +1;
  } *precisionRecord;
};

#endif // JULIAMODEL_H
