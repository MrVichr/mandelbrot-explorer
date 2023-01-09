#ifndef MANDELEVALUATOR_HPP
#define MANDELEVALUATOR_HPP

#include <QThread>

#include "MandelMath.hpp"
//#include "double_double.hpp"

enum NewtonNaiveChoice { nc03, nc05, nc08, ncWide, nc100, nc90, nc80, nc60, ncClose };

struct LaguerrePointStore
{
  enum State { stUnknown, stWorking, stResolved, stFail };
  std::atomic<State> state;
  double firstM;
  double firstStep_re, firstStep_im;
  NewtonNaiveChoice naiveChoice;
  int iter;
  LaguerrePointStore();
  void assign(const LaguerrePointStore *src);
};

template <class BASE>
struct LaguerrePoint
{
  enum IndexIntoWorker
  {
    iiw__BASE=0,
    iiw_r=iiw__BASE+0,
    iiw_fz_r=iiw__BASE+2,
    iiw_nth_fz=iiw__BASE+4,
    iiw__END=iiw__BASE+6
  };
  static constexpr size_t LEN=iiw__END-iiw__BASE;
  LaguerrePointStore *store;
  LaguerrePoint(LaguerrePointStore *store, MandelMath::NumberType ntype);
  MandelMath::complex<BASE> r;
  MandelMath::complex<BASE> fz_r;
  MandelMath::complex<BASE> nth_fz;
  LaguerrePoint &operator =(LaguerrePoint &src) = delete;
  //void assign(const LaguerrePoint<WORKER_MULTI, IIW_SRC> &src);
  void readFrom(void *storage, int index);
  void writeTo(void *storage, int index);
  void zero(const MandelMath::complex<BASE> *c);
};

struct MandelPointStore
{
  enum WorkState { stIdle, stWorking, stDone };
  std::atomic<WorkState> wstate;
  enum ResultState { stUnknown, stOutside, stOutAngle, stBoundary, stMisiur, stDiverge, stPeriod2, stPeriod3, stMaxIter };
  ResultState rstate;
  bool has_fc_r;
  int lookper_startiter, lookper_prevGuess, lookper_lastGuess;
  bool lookper_nearr_dist_touched; //check if equal but only once; in theory should only happen at dist=0
  int near0iter_1; //actual+1
  int nearziter_0;
  int period;
  int surehand;
  int iter, newton_iter;
  double exterior_hits, exterior_avoids; //upper and lower bound
  struct
  {
    double hits;
    double hits_re, hits_im;
    void zero() {hits=-1; hits_re=0; hits_im=0; }
  } interior;

  MandelPointStore();
  void assign(const MandelPointStore *src);
};

template <class BASE>
struct MandelPoint
{
  //static constexpr size_t LEN=20;
  enum IndexIntoWorker
  {
    iiw__BASE=0,
    iiw_f_=iiw__BASE+0,
    iiw_fc_c=iiw__BASE+2,
    iiw_fz_r=iiw__BASE+4,
    iiw_fz_c_mag=iiw__BASE+6,
    iiw_sure_fz_mag=iiw__BASE+7,
    iiw_sure_startf=iiw__BASE+8,
    iiw_lookper_startf=iiw__BASE+10,
    iiw_lookper_nearr=iiw__BASE+12,
    iiw_lookper_nearr_dist=iiw__BASE+14,
    iiw_lookper_totalFzmag=iiw__BASE+15,
    iiw_near0m=iiw__BASE+16,
    iiw_nearzm=iiw__BASE+17,
    iiw_root=iiw__BASE+18,
    iiw_extangle=iiw__BASE+20,
    iiw__END=iiw__BASE+21
  };
  static constexpr size_t LEN=iiw__END-iiw__BASE;
  MandelPointStore *store;
  MandelPoint(MandelPointStore *store, MandelMath::NumberType ntype);
  MandelMath::complex<BASE> f;
  MandelMath::complex<BASE> fc_c; //fc_c, or fz_r if stPeriod2 or stPeriod3
  MandelMath::complex<BASE> fz_r;
  MandelMath::number<BASE> fz_c_mag; //since the beginning
  MandelMath::number<BASE> sure_fz_mag; //since near0
  MandelMath::complex<BASE> sure_startf;
  MandelMath::complex<BASE> lookper_startf;
  MandelMath::complex<BASE> lookper_nearr;
  MandelMath::number<BASE> lookper_nearr_dist;
  MandelMath::number<BASE> lookper_totalFzmag;
  MandelMath::number<BASE> near0m;
  MandelMath::number<BASE> nearzm;
  MandelMath::complex<BASE> root;
  MandelMath::number<BASE> extangle;
  MandelPoint &operator =(MandelPoint &src) = delete;
  //template <int IIW_SRC>
  //void assign<WORKER_MULTI, IIW_OFFSET, IIW_SRC>(const MandelPoint<WORKER_MULTI, IIW_SRC> &src);
  //template <int IIW_SRC>
  //void assign(const MandelPoint<WORKER_MULTI, IIW_SRC> &src);
  void readFrom(void *storage, int index);//BASE src[LEN]);
  void writeTo(void *storage, int index);
  void zero(const MandelMath::complex<BASE> *first_z);
};

struct JuliaPointStore
{
  enum WorkState { stIdle, stWorking, stDone };
  std::atomic<WorkState> wstate;
  enum ResultState { stUnknown, stOutside, stOutAngle, stBoundary, stMisiur, stDiverge, stPeriod2, stPeriod3, stMaxIter };
  ResultState rstate;
  //bool has_fc_r;
  int lookper_startiter;//, lookper_prevGuess, lookper_lastGuess;
  //bool lookper_nearr_dist_touched; //check if equal but only once; in theory should only happen at dist=0
  int near0iter_1; //actual+1
  int nearziter_0;
  //int period;
  //int surehand;
  int iter;
  struct
  {
    double verify_hits, verify_avoids; //official estimate in one go
    double none_hits, none_avoids; //upper and lower bound; should==verify
    double always_hits, always_avoids;
    double condi_hits, condi_avoids;
    double blend_hits_, blend_avoids;
    void zero(double x) {verify_hits=x; verify_avoids=x; none_hits=x; none_avoids=x; always_hits=x; always_avoids=x; condi_hits=x; condi_avoids=x; blend_hits_=x; blend_avoids=x; }
  } exterior;
  struct
  {
    double hits;
    double hits_re, hits_im;
    double phi_re, phi_im;
    double phi1_re, phi1_im;
    double phi2_re, phi2_im;
    double f1_re, f1_im;
    int first_under_1;
    int cycles_until_root;
    void zero() {hits=-1; hits_re=0; hits_im=0; phi_re=0; phi_im=0; phi1_re=0; phi1_im=0; phi_re=0; phi2_im=0; f1_re=0; f1_im=0; first_under_1=0; cycles_until_root=0; }
  } interior;

  JuliaPointStore();
  void assign(const JuliaPointStore *src);
};

template <class BASE>
struct JuliaPoint
{
  //static constexpr size_t LEN=20;
  enum IndexIntoWorker
  {
    iiw__BASE=0,
    iiw_f_=iiw__BASE+0,
    iiw_fc_c=iiw__BASE+2,
    iiw_fz_c_mag=iiw__BASE+4,
    iiw_lookper_distr=iiw__BASE+5,
    iiw_lookper_fz=iiw__BASE+7,
    iiw_near0m=iiw__BASE+9,
    iiw_nearzm=iiw__BASE+10,
    iiw_near0fzm=iiw__BASE+11,
    iiw_since0fzm=iiw__BASE+12,
    iiw_extangle=iiw__BASE+13,
    iiw__END=iiw__BASE+14
  };
  static constexpr size_t LEN=iiw__END-iiw__BASE;
  JuliaPointStore *store;
  JuliaPoint(JuliaPointStore *store, MandelMath::NumberType ntype);
  MandelMath::complex<BASE> f;
  MandelMath::complex<BASE> fz_z;
  MandelMath::number<BASE> fz_z_mag; //since the beginning
  MandelMath::complex<BASE> lookper_distr;
  MandelMath::complex<BASE> lookper_fz;
  MandelMath::number<BASE> nearzm;
  MandelMath::number<BASE> near0fzm; //from start to near0
  MandelMath::number<BASE> near0m;
  MandelMath::number<BASE> since0fzm; //from near0+1 to end
  MandelMath::number<BASE> extangle;
  JuliaPoint &operator=(JuliaPoint &src) = delete;
  //template <int IIW_SRC>
  //void assign<WORKER_MULTI, IIW_OFFSET, IIW_SRC>(const MandelPoint<WORKER_MULTI, IIW_SRC> &src);
  //template <int IIW_SRC>
  //void assign(const MandelPoint<WORKER_MULTI, IIW_SRC> &src);
  void readFrom(void *storage, int index);//BASE src[LEN]);
  void writeTo(void *storage, int index);
  void zero(const MandelMath::complex<BASE> *first_z);
};

class ShareableViewInfo: public QObject
{
  Q_OBJECT
protected:
  int refcount_unused;
public:
  //static constexpr int LEN=5;
  ShareableViewInfo(): view(), c(), nth_fz_limit() { } //Qt uses this and operator= instead of copy constructor :-/
  /*template <typename BASE>
  ShareableViewInfo(): c(new MandelMath::number<BASE>(), new MandelMath::number<BASE>()),
                       root(new MandelMath::number<BASE>(), new MandelMath::number<BASE>()),
                       nth_fz_limit(new MandelMath::number<BASE>()), scale(1), period(0), nth_fz(0) { }*/
  ShareableViewInfo(MandelMath::NumberType ntype): view(ntype), c(ntype), nth_fz_limit(ntype), scale(1), nth_fz(0), max_root_effort(3) { }
  ShareableViewInfo(ShareableViewInfo &src);
  ShareableViewInfo(const ShareableViewInfo &src);
  ShareableViewInfo(ShareableViewInfo &&src); //important
  //~ShareableViewInfo();
  ShareableViewInfo &operator=(ShareableViewInfo &src);
  ShareableViewInfo &operator=(ShareableViewInfo &&src);
  MandelMath::complex<MandelMath::number_a *> view;
  MandelMath::complex<MandelMath::number_a *> c;
  MandelMath::number<MandelMath::number_a *> nth_fz_limit;
  double scale;
  int nth_fz;
  int max_root_effort;
};
Q_DECLARE_METATYPE(ShareableViewInfo);

template <class BASE>
class LaguerreStep
{
protected:
  /*enum IndexIntoWorker
  {
    iiw__BASE=133,
    iiw_step=iiw__BASE+0,
    iiw_s1=iiw__BASE+2,
    iiw_s2=iiw__BASE+4,
    iiw_tmp1=iiw__BASE+6,
    iiw_tmp2=iiw__BASE+8,
    iiw_laguG=iiw__BASE+9,
    iiw_laguG2=iiw__BASE+11,
    iiw_laguH=iiw__BASE+13,
    iiw_laguX=iiw__BASE+15,
    iiw_fzzf=iiw__BASE+17,
    iiw__END=iiw__BASE+19
  };*/
public:
  //static constexpr size_t LEN=19;
  LaguerreStep(MandelMath::NumberType ntype);
  //~LaguerreStep();
  //void switchType(MandelMath::number_worker *worker);
  //do one Laguerre step
  bool eval(typename MandelMath::complex<BASE>::Scratchpad *tmp,
            int lg2_degree, const MandelMath::complex<BASE> *const f,
                            const MandelMath::complex<BASE> *const f_z,
                            const MandelMath::complex<BASE> *const f_zz); //->step_re, step_im

  //results
  MandelMath::complex<BASE> step;
  struct
  {
    double mum_re, mum_im;
    double mu_re, mu_im;
    int lastm;
  } dbg;

  //temporary
protected:
  MandelMath::complex<BASE> s1; //TODO: unused
  MandelMath::complex<BASE> s2;
  MandelMath::complex<BASE> tmp1;
  MandelMath::number<BASE> tmp2;
  MandelMath::complex<BASE> laguG;
  MandelMath::complex<BASE> laguG2;
  MandelMath::complex<BASE> laguH;
  MandelMath::complex<BASE> laguX;
  MandelMath::complex<BASE> fzzf;
};

template <typename BASE>
class MandelLoopEvaluator
{
public:
  /*enum IndexIntoWorker
  {
    iiw__BASE=75,
    iiw_help_c=iiw__BASE+0,
    iiw_f=iiw__BASE+2,
    iiw_f_z=iiw__BASE+4,
    iiw_f_c=iiw__BASE+6,
    iiw_f_zz=iiw__BASE+8,
    iiw_f_zc=iiw__BASE+10,
    iiw_f_cc=iiw__BASE+12,
    iiw_f_zzc=iiw__BASE+14,
    iiw_first_multi=iiw__BASE+16,
    iiw_sumA=iiw__BASE+18,
    iiw_fz_mag=iiw__BASE+20,
    iiw_fz_minmag=iiw__BASE+21,
    iiw_s1=iiw__BASE+22,
    iiw_s2=iiw__BASE+24,
    iiw__END=iiw__BASE+26
  };*/
  //static constexpr size_t LEN=26;
  MandelLoopEvaluator(MandelMath::NumberType ntype);
  //~MandelLoopEvaluator();
  /*union Place
  {
    MandelMath::dd_real (*dd)[LEN];
    MandelMath::multiprec (*multi)[LEN];
    Place(): dd(nullptr) { }
  } place;*/
  //void switchType(MandelMath::number_worker *worker);
  //eval F_c(c)-c
  bool evalg(typename MandelMath::complex<BASE>::Scratchpad *tmp,
             int period, const MandelMath::complex<BASE> *c); //->f, f_c, f_cc
  //eval F_c(z) -z if minusZ
  bool eval2(typename MandelMath::complex<BASE>::Scratchpad *tmp,
             int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z, bool minusZ, bool doInit=true);     //->f, f_z, f_zz, f_c, f_zc, f_cc
  bool eval2_mag(typename MandelMath::complex<BASE>::Scratchpad *tmp,
                 int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z); //->f, f_z, f_zz, f_c, f_zc, f_z_mag
  bool eval_zz(typename MandelMath::complex<BASE>::Scratchpad *tmp,
               int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z, bool minusZ, bool doInit=true);   //->f, f_z, f_zz
  bool eval2zzc(typename MandelMath::complex<BASE>::Scratchpad *tmp,
                int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z);  //->f, f_z, f_zz, f_c, f_zc, f_cc, f_zzc
  bool eval_multi(typename MandelMath::complex<BASE>::Scratchpad *tmp,
                  int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z, const MandelMath::complex<BASE> *f_z_target, double dist_tolerance); //->f, f_z, f_zz, multi, first_multi
  bool eval_near0(typename MandelMath::complex<BASE>::Scratchpad *tmp,
                  int period, const MandelMath::complex<BASE> *c); //->near0iter_1
  bool eval_ext_mandMJ(typename MandelMath::complex<BASE>::Scratchpad *tmp, double mandel,
                    const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *z, int iter); //

  //inputs
  //MandelMath::complex<BASE> help_c; //for use by caller if needed
    // TODO: back into params.c now that we have first_z ?
  //MandelMath::number_store z_re, z_im;
  //MandelMath::number_store c_re, c_im;

  //results
  MandelMath::complex<BASE> f;
  MandelMath::complex<BASE> f_z;
  MandelMath::complex<BASE> f_c;
  MandelMath::complex<BASE> f_zz;
  MandelMath::complex<BASE> f_zc;
  MandelMath::complex<BASE> f_cc;
  MandelMath::complex<BASE> f_zzc;
  int multi;
  MandelMath::complex<BASE> first_multi;
  int near0iter_1;
  MandelMath::complex<BASE> sumA;
  MandelMath::number<BASE> f_z_mag;
  MandelMath::number<BASE> f_z_minmag; //also near0mag

  //temporary
protected:
  MandelMath::complex<BASE> s1;
  MandelMath::complex<BASE> s2;
};

template <class WORKER_MULTI>
class MandelEvaluator;

class MandelEvaluatorThread: public QThread
{
  Q_OBJECT
protected:
  MandelEvaluator<MandelMath::number_a *> *owner;
public:
  MandelEvaluatorThread(MandelEvaluator<MandelMath::number_a *> *owner);
  void syncMandel();
  void syncJulia(int juliaPeriod);
  int syncLaguerre(int period, MandelMath::complex<MandelMath::number_a *> *r, const bool fastHoming);
public slots:
  void doMandelThreaded(int epoch);
  void doJuliaThreaded(int epoch, int juliaPeriod);
  void doLaguerreThreaded(int epoch, int laguerrePeriod);
signals:
  void doneMandelThreaded(MandelEvaluatorThread *me);
  void doneJuliaThreaded(MandelEvaluatorThread *me);
  void doneLaguerreThreaded(MandelEvaluatorThread *me);
};

template <class BASE>
class MandelEvaluator //templates cannot be Q_OBJECT (or QThread) and cannot receive signals
{
protected:
  //MandelMath::worker_multi::Allocator self_allocator;
public:
  constexpr static double LARGE_FLOAT2=1e60;
  constexpr static double MAGIC_MIN_SHRINK=1.5;
  constexpr static int MAX_PERIOD=8000;
  /*enum IndexIntoWorker
  {
    iiw_currentData=2,
    iiw_currentParams=28,
    iiwLaguerrePoint=100,
  };*/
  int threaded_errorcode;
  int busyEpoch;
  MandelEvaluatorThread thread;
  //WORKER_MULTI *currentWorker;
  MandelMath::NumberType ntype; //typeEmpty for number * ?
  MandelEvaluator(MandelMath::NumberType ntype, bool dontRun);
  virtual ~MandelEvaluator(); //must be virtual, or would have to typecast in delete threads[i]
  void startRunning();
/*#if NUMBER_DOUBLE_EXISTS
  static void simple_double(double cr, double ci, MandelPoint<MandelMath::worker_multi> &data, int maxiter);
#endif //NUMBER_DOUBLE_EXISTS
  static void simple_ddouble(MandelMath::dd_real *cr, MandelMath::dd_real *ci, MandelPoint &data, int maxiter);
  static void simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter);*/

  //MandelMath::number_worker::Type currentType;
  //MandelMath::worker_multi::Allocator *currentAllocator;
  //void switchType(MandelMath::number_worker *worker);
  int workIfEpoch;
  int pointsComputed;
  qint64 totalNewtonIterations;
  QElapsedTimer timeThreaded;
  qint64 timeThreadedTotal;
  QElapsedTimer timeOuter;
  qint64 timeOuterTotal;
  QElapsedTimer timeInner;
  qint64 timeInnerTotal;
  QElapsedTimer timeInvoke;
  qint64 timeInvokePostTotal;
  qint64 timeInvokeSwitchTotal;

  struct ComputeParams
  {
    /*enum IndexIntoWorker
    {
      iiw__BASE=0,
      iiw_c=iiw__BASE+0,
      iiw_r=iiw__BASE+2,
      iiw_z=iiw__BASE+4,
      iiw__END=iiw__BASE+6
    };
    static constexpr size_t LEN=iiw__END-iiw__BASE;*/
    MandelMath::complex<BASE> c;
    MandelMath::complex<BASE> juliaRoot;
    double patchSizeExterior;
    MandelMath::number<BASE> juliaBailout_; //not quite 4
    MandelMath::complex<BASE> juliaAlpha; // |alpha|<1
    MandelMath::complex<BASE> first_z;
    int epoch;
    int pixelIndex;
    int maxiter;
    bool breakOnNewNearest;
    bool want_fc_r;
    bool want_extangle;
    int nth_fz;
    ComputeParams(MandelMath::NumberType ntype);
  } currentParams;
  LaguerrePointStore laguerreStore;
  LaguerrePoint<BASE> laguerreData;
  MandelPointStore mandelDataStore;
  MandelPoint<BASE> mandelData;
  JuliaPointStore juliaStore;
  JuliaPoint<BASE> juliaData;
  typename MandelMath::complex<BASE>::Scratchpad tmp;

  void syncMandelSplit();
  void syncJuliaSplit(int juliaPeriod);
  MandelLoopEvaluator<BASE> loope;
  int laguerre(int period, const MandelMath::complex<BASE> *c, MandelMath::complex<BASE> *r, const bool fastHoming);
  struct NewtRes
  {
    //static constexpr int LEN=8;
    int cyclesNeeded;
    //now firstMum_re double firstM;
    MandelMath::complex<BASE> fz_r; //TODO: would be nice to move into newt and move newton() to newt
    MandelMath::complex<BASE> nth_fz;
    MandelMath::complex<BASE> first_guess_lagu;
    MandelMath::complex<BASE> first_guess_newt;
    double first_fejer_re, first_fejer_im;
    double first_naive1_re, first_naive1_im, first_naive2_re, first_naive2_im, first_naive_re, first_naive_im;
    NewtonNaiveChoice naiveChoice;
    double first_neumaier1_re, first_neumaier1_im; //circle enclosing 1 but not 2 roots (approx)
    double first_neumaier2_re, first_neumaier2_im;   //circle enclosing 2 but not 3 roots (approx) (never valid without  f''')
    double ostrowski_r1c, ostrowski_r1x, ostrowski_r2c;
    //double first_lagum_re, first_lagum_im;
    double first_lagu1_re, first_lagu1_im, first_lagu1o_re, first_lagu1o_im;
    double firstMu_re, firstMu_im, firstMum_re, firstMum_im;
    double accy_tostop, accy_multiplier, accy_noise; //in units of eps2()
    NewtRes(MandelMath::NumberType ntype);
  } newtres;
protected:
  struct Eval
  {
    /*enum IndexIntoWorker
    {
      iiw_fz_mag=42,
      iiw_near0fmag=43,
    };*/
    //static constexpr int LEN=3;
    MandelMath::complex<BASE> fz_r;
    MandelMath::number<BASE> fz_mag;
    //MandelMath::number<BASE> near0fmag;
    Eval(MandelMath::NumberType ntype);
  } eval;
  struct Newt
  {
    /*enum IndexIntoWorker
    {
      iiw_tmp2=68,
    };
    static constexpr int LEN=25;*/
    MandelMath::complex<BASE> bestr;
    MandelMath::complex<BASE> f_r;
    //MandelMath::number_store fz_r_re, fz_r_im;
    MandelMath::complex<BASE> fzz_r;
    MandelMath::complex<BASE> tmp1;
    //MandelMath::number_store fzfix_re, fzfix_im;
    MandelMath::complex<BASE> laguH;
    MandelMath::complex<BASE> laguG;
    MandelMath::complex<BASE> laguG2;
    MandelMath::complex<BASE> laguX;
    MandelMath::complex<BASE> newtX;
    MandelMath::complex<BASE> prevR;
    MandelMath::complex<BASE> prevGz;
    MandelMath::complex<BASE> fzzf;
    MandelMath::number<BASE> tmp2;
    Newt(MandelMath::NumberType ntype);
  } newt;
  struct InteriorInfo
  {
    /*enum IndexIntoWorker
    {
      iiw_inte_abs=71,
    };*/
    //static constexpr int LEN=6;
    MandelMath::complex<BASE> inte;
    MandelMath::number<BASE> inte_abs;
    MandelMath::complex<BASE> fz;
    MandelMath::number<BASE> fz_mag;
    MandelMath::complex<BASE> alpha;
    MandelMath::complex<BASE> alphak;
    MandelMath::complex<BASE> zoom;
    InteriorInfo(MandelMath::NumberType ntype);
  } interior;
public:
  struct Bulb
  {
    MandelLoopEvaluator<BASE> *loope;
    MandelMath::complex<BASE> baseZC;
    MandelMath::complex<BASE> baseCC;
    MandelMath::complex<BASE> t1;
    MandelMath::complex<BASE> t2;
    MandelMath::complex<BASE> t3;
    MandelMath::complex<BASE> deltac;
    MandelMath::complex<BASE> deltar;
    MandelMath::complex<BASE> first_cb;
    double firstStep;
    MandelMath::complex<BASE> dbg_first_rb;
    MandelMath::complex<BASE> prev_rb;
    MandelMath::complex<BASE> target_f_z;

    MandelMath::complex<BASE> res_cb;
    MandelMath::complex<BASE> res_rb;
    MandelMath::complex<BASE> res_xc;
    MandelMath::complex<BASE> res_baseZC;
    MandelMath::complex<BASE> res_baseCC;
    bool res_card;
    int res_foundMult;
    bool res_valid;
    MandelMath::complex<BASE> res_baseFz;
    //MandelMath::number_store test_x0_re, test_x0_im;
    //MandelMath::number_store test_xn_re, test_xn_im;
    LaguerreStep<BASE> lagu;
    Bulb(MandelMath::NumberType ntype, MandelLoopEvaluator<BASE> *loope);
    ~Bulb() { }
  protected:
    //void fixRnearBase(MandelMath::complex_place<WORKER_MULTI> *r, const MandelMath::complex_place<WORKER_MULTI> *c, int period, int *mult);
  public:
    void findBulbBase(typename MandelMath::complex<BASE>::Scratchpad *tmp, int period2, const MandelMath::complex<BASE> *c);//, MandelMath::complex<WORKER_MULTI> *cb, MandelMath::complex<WORKER_MULTI> *rb, MandelMath::complex<WORKER_MULTI> *xc, MandelMath::complex<WORKER_MULTI> *baseZC, MandelMath::complex<WORKER_MULTI> *baseCC, bool *is_card, int *foundMult, MandelMath::complex<WORKER_MULTI> *baseFz);
    //static constexpr int LEN=MandelLoopEvaluator<BASE>::LEN+32+LaguerreStep<BASE>::LEN;
  } bulb;
  struct ExtAngle
  {
    static constexpr double SAFE_RADIUS=6;
    using number=MandelMath::number<BASE>;
    using complex=MandelMath::complex<BASE>;
    static constexpr double SPECIAL_VALUE_EDGE=10;
    static constexpr double SPECIAL_VALUE_DEEP=11;
    MandelEvaluator<BASE> *owner;
    MandelMath::NumberType ntype;
    MandelLoopEvaluator<BASE> *loope;
    LaguerreStep<BASE> *lagus;
    MandelMath::complex<BASE> z;
    MandelMath::complex<BASE> r;
    MandelMath::number<BASE> angleC;
    MandelMath::number<BASE> angle;
    MandelMath::complex<BASE> x;
    MandelMath::complex<BASE> target;
    MandelMath::complex<BASE> dldlz;
    MandelMath::complex<BASE> d2ldlz2;
    struct Dbg
    {
      int first_guess_valid;
      complex first_guess_0, first_guess_1_, first_guess_2;
      int last_guess_valid;
      complex last_guess_0, last_guess_1_, last_guess_2;
      Dbg(MandelMath::NumberType ntype): first_guess_valid(0), first_guess_0(ntype), first_guess_1_(ntype), first_guess_2(ntype),
                                         last_guess_valid(0), last_guess_0(ntype), last_guess_1_(ntype), last_guess_2(ntype) {}
    } dbg;
    ExtAngle(MandelMath::NumberType ntype, MandelLoopEvaluator<BASE> *loope, LaguerreStep<BASE> *lagus, MandelEvaluator<BASE> *owner);
    ~ExtAngle();
    void computeMJ(number *result, bool mandel, int iter, const complex *c, const complex *z, typename complex::Scratchpad *tmp);
  } extangle;
  //static constexpr size_t LEN=10+ ComputeParams::LEN+MandelPoint_<BASE>::LEN+LaguerrePoint<BASE>::LEN+JuliaPoint<BASE>::LEN+NewtRes::LEN+Eval::LEN+Newt::LEN+InteriorInfo::LEN+Bulb::LEN;
public:
  struct
  {
    std::function<int (MandelEvaluator *me)> give;
    std::function<int (MandelEvaluator *me, int result, bool giveWork)> doneLaguerre;
    std::function<int (MandelEvaluator *me, bool giveWork)> doneMandel;
    std::function<int (MandelEvaluator *me, bool giveWork)> doneJulia;
  } threaded;
protected:

  protected:
  int periodCheck(int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *root_seed, bool exactMatch);
  int estimateInterior(int period, const MandelMath::complex<BASE> *c, const MandelMath::complex<BASE> *root);//, InteriorInfo *interior);
  void mandel_until_bailout();
  void julia_until_bailout();
  void evaluateMandel();
  void doMandelThreadedSplit(int epoch);
  void evaluateJulia(int juliaPeriod);
  void doJuliaThreadedSplit(int epoch, int juliaPeriod);
  void doLaguerreThreadedSplit(int epoch, int laguerrePeriod);
  friend MandelEvaluatorThread;
/*protected slots:
  void doCompute();
  //void doNewton();
signals:
  void doneCompute(MandelEvaluator *me);
  void doneNewton(MandelEvaluator *me, int result);*/
};

#endif // MANDELEVALUATOR_HPP
