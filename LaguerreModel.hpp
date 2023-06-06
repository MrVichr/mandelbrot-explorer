#ifndef LAGUERREMODEL_H
#define LAGUERREMODEL_H

#include <QObject>
#include <QImage>
#include <QReadWriteLock>

#include "atomic_update_rect.hpp"
#include "ShareableImageWrapper.hpp"
#include "MandelEvaluator.hpp"

class LaguerreModel: public QObject
{
  Q_OBJECT
protected:
  enum IndexIntoWorker
  {
    iiw_reImToPixel=10
  };
public:
  LaguerreModel();
  ~LaguerreModel();
  Q_INVOKABLE void startRunning();
  /*template <class WORKER_MULTI>
  void transformStore(WORKER_MULTI *old_worker, WORKER_MULTI *old_sworker, LaguerrePointStore *old_store, int old_width, int old_height, const MandelMath::complex<WORKER_MULTI> *old_c,
                      WORKER_MULTI *new_worker, WORKER_MULTI *new_sworker, LaguerrePointStore *new_store, int new_width, int new_height, const MandelMath::complex<WORKER_MULTI> *new_c,
                      int inlog, int new_step_log);*/
  void transformStore(void *old_points, LaguerrePointStore *old_store, int old_width, int old_height, MandelMath::complex<MandelMath::number_any> const &old_c,
                      void *new_points, LaguerrePointStore *new_store, int new_width, int new_height, MandelMath::complex<MandelMath::number_any> const &new_c,
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
    static constexpr int LEN=5;
    int period;
    int nth_fz;
    MandelMath::complex<MandelMath::number_any> c;
    MandelMath::complex<MandelMath::number_any> root;
    MandelMath::number<MandelMath::number_any> nth_fz_limit;
    Params(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Params *source);
  };


  enum paintStyle
  {
    paintStyleCls=0,
    paintStyleNthFz=1,
    paintStyleNthFz1=2,
    paintStyleRoots=3,
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

protected:
  QReadWriteLock threading_mutex;
  template <typename BASE>
  int giveWorkThreaded(MandelEvaluator<BASE> *me);
  template <typename BASE>
  int doneWorkThreaded(MandelEvaluator<BASE> *me, int result, bool giveWork);

  AtomicUpdateRect image_dirty;
  void invalidateMainImage()
  {
    image_dirty.invalidate();
  }
public slots:
  void doneWorkInThread(MandelEvaluatorThread *me);
  void selectedPrecisionChanged();
signals:
  void selectedPaintStyleChanged();
  void selectedPrecisionChange();
  void triggerLaguerreThreaded(int epoch);
protected:
  //MandelMath::worker_multi::Allocator<MandelMath::worker_multi> *storeAllocator;
  //MandelMath::worker_multi *storeWorker; //pointStore
  LaguerrePointStore *pointStore;
  //constexpr static int MAX_ZOOM_IN_DOUBLE=55;//53;
  //MandelMath::number_store::DbgType currentMath;
  int epoch;
  int imageWidth;
  int imageHeight;
  int nextGivenPointIndex;
  int effortBonus;
  //constexpr static int MAX_EFFORT=17;//131072 iters;
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
    //MandelMath::worker_multi::Allocator pointAllocator;
    //-> evaluator.tmpLaguerreStore
    LaguerrePointStore pointDataStore; //but we need it for long-term storage of results
    //-> evaluator.tmpLaguerrePoint ? LaguerrePoint<MandelMath::worker_multi> pointData;
    LaguerrePoint<MandelMath::number_any> pointData;
    double first_mu_re, first_mu_im, first_mum_re, first_mum_im;
    Orbit(MandelMath::NumberType ntype, Orbit const *source);
    ~Orbit();
    constexpr static int LEN=LaguerrePoint<MandelMath::number_any>::LEN;
  };
  struct PrecisionRecord
  {
    MandelMath::NumberType ntype;
    Orbit orbit;
    LaguerrePoint<MandelMath::number_any> wtiPoint;
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
    PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, LaguerreModel *doneReceiver);
    ~PrecisionRecord();
    //constexpr static int LEN=Params::LEN+LaguerrePoint<MandelMath::number_a *>::LEN+Position::LEN+Orbit::LEN  +4   +1;
  } *precisionRecord;
};

#endif // LAGUERREMODEL_H
