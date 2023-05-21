#include "MandelEvaluator.hpp"
#define _USE_MATH_DEFINES //some magic
#include <cmath>
//C++20 also don't have #include <numbers>
//has gcd since C++17 which I apparently don't have #include <numeric>

#define USE_GCD_FOR_CHECKPERIOD 0
#define CLEVER_FIX 0
#define THREADED_DONE_GIVE_WORK 0
//TODO: delete surehand
#define SUREHAND_CHOICE 0 //0..old only, 1..surehand only, 2..both
  //1 tries newton fewer times
  //0 tries newton with smaller periods that are sometimes correct, thanks to restarts

void nop(); //defined in MandelMath

LaguerrePointStore::LaguerrePointStore(): state(State::stUnknown),
  firstM(0), firstStep_re(0), firstStep_im(0), iter(0)
{

}

void LaguerrePointStore::assign(const LaguerrePointStore *src)
{
  working_assert(src!=nullptr);
  //TODO: if (src==nullptr) see LaguerrePoint::zero
  state=src->state.load();
  firstM=src->firstM;
  firstStep_re=src->firstStep_re;
  firstStep_im=src->firstStep_im;
  iter=src->iter;
}


template<typename BASE>
LaguerrePoint<BASE>::LaguerrePoint(LaguerrePointStore *store, MandelMath::number<BASE>::Scratchpad *spad):
  store(store), r(spad), fz_r(spad), nth_fz(spad)
{
}

template<typename BASE>
void LaguerrePoint<BASE>::readFrom(void *storage, int index)
{
  //store->assign(src.store);
  r.readFrom(storage, index+iiw_r);
  fz_r.readFrom(storage, index+iiw_fz_r);
  nth_fz.readFrom(storage, index+iiw_nth_fz);
}

template<class BASE>
void LaguerrePoint<BASE>::writeTo(void *storage, int index)
{
  r.writeTo(storage, index+iiw_r);
  fz_r.writeTo(storage, index+iiw_fz_r);
  nth_fz.writeTo(storage, index+iiw_nth_fz);
}

template<class BASE>
void LaguerrePoint<BASE>::zero(const MandelMath::complex<BASE> &c)
{
  store->state=LaguerrePointStore::State::stUnknown;
  store->firstM=0;
  store->firstStep_re=0;
  store->firstStep_im=0;
  store->iter=0;
  r.assign_across(c);
  fz_r.zero(1, 0);
  nth_fz.zero();
}

template <typename BASE>
ComputeMandelParams<BASE>::ComputeMandelParams(MandelMath::number<BASE>::Scratchpad *spad):
  c(spad), bailout(spad), first_z(spad)
{
  bailout.zero(4);
}

template<typename BASE> template<typename OTHER_BASE>
void ComputeMandelParams<BASE>::assign_across(const ComputeMandelParams<OTHER_BASE> &src)
{
  c.assign_across(src.c);
  bailout.assign_across(src.bailout);
  //no need to copy first_z
  //first_z_.assign_across(&src.first_z_);
}

MandelPointStore::MandelPointStore(): wstate(WorkState::stIdle), rstate(ResultState::stUnknown), iter(0)
{

}

void MandelPointStore::assign(const MandelPointStore *src)
{
  wstate.store(src->wstate.load());
  rstate=src->rstate;
  iter=src->iter;
  has_fc_r=src->has_fc_r;
  lookper_startiter=src->lookper_startiter;
  lookper_prevGuess=src->lookper_prevGuess;
  lookper_lastGuess=src->lookper_lastGuess;
  lookper_nearr_dist_touched=src->lookper_nearr_dist_touched;
  near0iter_1=src->near0iter_1;
  nearziter_0=src->nearziter_0;
  newton_iter=src->newton_iter;
  period=src->period;
  surehand=src->surehand;
  exterior_hits=src->exterior_hits;
  exterior_avoids=src->exterior_avoids;
  interior=src->interior;
}

template<typename BASE>
MandelPoint<BASE>::MandelPoint(MandelPointStore *store, MandelMath::number<BASE>::Scratchpad *spad):
  store(store), f(spad), fc_c(spad), fz_r(spad), fz_c_mag(spad),
  sure_fz_mag(spad), sure_startf(spad), lookper_startf(spad), lookper_nearr(spad),
  lookper_nearr_dist(spad), lookper_totalFzmag(spad),
  near0m(spad), nearzm(spad), root(spad), extangle(spad)
{
}

/*
what's the problem here?
template<int IIW_OFFSET>
MandelPoint<worker_multi, IIW_OFFSET>::MandelPoint(MandelPointStore *store, MandelMath::worker_multi::Allocator<worker_multi> *allocator):
  self_allocator(allocator, LEN),
  store(store), f(&self_allocator), fc_c(&self_allocator), fz_r(&self_allocator), fz_c_mag(&self_allocator),
  sure_fz_mag(&self_allocator), sure_startf(&self_allocator), lookper_startf_(&self_allocator), lookper_nearr(&self_allocator),
  lookper_nearr_dist(&self_allocator), lookper_totalFzmag(&self_allocator),
  near0f(&self_allocator), root(&self_allocator)
{
  //reset();
  //should be overwritten before read:
  working_assert(self_allocator.checkFill());
}
*/

template <typename BASE>
void MandelPoint<BASE>::readFrom(void *storage, int index)
{
  //store->assign(src.store);
  f.readFrom(storage, index+iiw_f_);
  fc_c.readFrom(storage, index+iiw_fc_c);
  fz_r.readFrom(storage, index+iiw_fz_r);
  fz_c_mag.readFrom(storage, index+iiw_fz_c_mag);
  sure_fz_mag.readFrom(storage, index+iiw_sure_fz_mag);
  sure_startf.readFrom(storage, index+iiw_sure_startf);
  lookper_startf.readFrom(storage, index+iiw_lookper_startf);
  lookper_nearr.readFrom(storage, index+iiw_lookper_nearr);
  lookper_nearr_dist.readFrom(storage, index+iiw_lookper_nearr_dist);
  lookper_totalFzmag.readFrom(storage, index+iiw_lookper_totalFzmag);
  //MandelMath::number_instance<WORKER_MULTI>(&lookper_nearr_dist).assign(src.lookper_nearr_dist.ptr);
  //MandelMath::number_instance<WORKER_MULTI>(&lookper_totalFzmag).assign(src.lookper_totalFzmag.ptr);
  near0m.readFrom(storage, index+iiw_near0m);
  nearzm.readFrom(storage, index+iiw_nearzm);
  root.readFrom(storage, index+iiw_root);
  extangle.readFrom(storage, index+iiw_extangle);
}

template <typename BASE>
void MandelPoint<BASE>::writeTo(void *storage, int index)
{
  //store->assign(src.store);
  f.writeTo(storage, index+iiw_f_);
  fc_c.writeTo(storage, index+iiw_fc_c);
  fz_r.writeTo(storage, index+iiw_fz_r);
  fz_c_mag.writeTo(storage, index+iiw_fz_c_mag);
  sure_fz_mag.writeTo(storage, index+iiw_sure_fz_mag);
  sure_startf.writeTo(storage, index+iiw_sure_startf);
  lookper_startf.writeTo(storage, index+iiw_lookper_startf);
  lookper_nearr.writeTo(storage, index+iiw_lookper_nearr);
  lookper_nearr_dist.writeTo(storage, index+iiw_lookper_nearr_dist);
  lookper_totalFzmag.writeTo(storage, index+iiw_lookper_totalFzmag);
  //MandelMath::number_instance<WORKER_MULTI>(&lookper_nearr_dist).assign(src.lookper_nearr_dist.ptr);
  //MandelMath::number_instance<WORKER_MULTI>(&lookper_totalFzmag).assign(src.lookper_totalFzmag.ptr);
  near0m.writeTo(storage, index+iiw_near0m);
  nearzm.writeTo(storage, index+iiw_nearzm);
  root.writeTo(storage, index+iiw_root);
  extangle.writeTo(storage, index+iiw_extangle);
}

template <typename BASE>
void MandelPoint<BASE>::zero(const MandelMath::complex<BASE> &first_z)
{
  f.assign_across(first_z);

  root.im.assign(f.im); //don't have "tmp" here so do it the long way... for now
  root.im.sqr();
  near0m.assign(f.re);
  near0m.sqr();
  near0m.add(root.im);

  fc_c.zero(1, 0);
  fz_r.zero(1, 0);
  fz_c_mag.zero(1);
  sure_fz_mag.zero(1);
  sure_startf.assign(f);
  store->lookper_prevGuess=0;
  store->lookper_lastGuess=0;
  //lookper resets at first iter
  store->near0iter_1=1;
  store->nearziter_0=1;

  //near0m.assign(f.getMag_tmp(tmp));

  //nearzm_.assign(near0m_); //could start with infinity but we should eventually iterate closer than 0... I think
  nearzm.zero(20);
  store->period=0;
  store->surehand=0;
  root.zero(0, 0);
  extangle.zero(0);

  store->wstate=MandelPointStore::WorkState::stIdle;
  store->rstate=MandelPointStore::ResultState::stUnknown;
  store->iter=0;
  store->newton_iter=0;
  store->exterior_avoids=-1;
  store->exterior_hits=-1;
  store->interior.zero();
  store->has_fc_r=false;
  /*
    real exterior:=0
    real interior:=0
    initwinding(c)
    complex interiorComplex:=0
    int period:=0
    complex root:=0
  */
}

template <typename BASE>
ComputeJuliaParams<BASE>::ComputeJuliaParams(MandelMath::number<BASE>::Scratchpad *spad):
  period(1), //must agree with c==0+0i
  c(spad), root(spad), root0(spad), patchSizeExterior(0), bailout(spad), alphaShort(spad), alpha(spad), first_z(spad)
{
  bailout.zero(4);
}

template<typename BASE> template<typename OTHER_BASE>
void ComputeJuliaParams<BASE>::assign_across(const ComputeJuliaParams<OTHER_BASE> &src)
{
  period=src.period;
  c.assign_across(src.c);
  root.assign_across(src.root);
  root0.assign_across(src.root0);
  alphaShort.assign_across(src.alphaShort);
  alpha.assign_across(src.alpha);
  patchSizeExterior=src.patchSizeExterior;
  bailout.assign_across(src.bailout);
  //no need to copy first_z
}

JuliaPointStore::JuliaPointStore(): wstate(WorkState::stIdle), rstate(ResultState::stUnknown), nearmr_difficult(false), iter(0)
{

}

void JuliaPointStore::assign(const JuliaPointStore *src)
{
  wstate.store(src->wstate.load());
  rstate=src->rstate;
  iter=src->iter;
  //has_fc_r=src->has_fc_r;
  lookper_startiter=src->lookper_startiter;
  //lookper_prevGuess=src->lookper_prevGuess;
  //lookper_lastGuess=src->lookper_lastGuess;
  //lookper_nearr_dist_touched=src->lookper_nearr_dist_touched;
  near0iter_1=src->near0iter_1;
  nearziter_0=src->nearziter_0;
  bigfzfzziter=src->bigfzfzziter;
  nearmriter=src->nearmriter;
  nearmr_difficult=src->nearmr_difficult;
  //newton_iter=src->newton_iter;
  //period=src->period;
  //surehand=src->surehand;
  exterior=src->exterior;
  interior=src->interior;
}

template<typename BASE>
JuliaPoint<BASE>::JuliaPoint(JuliaPointStore *store, MandelMath::number<BASE>::Scratchpad *spad):
  store(store), f(spad), fz_z(spad), fzz_z(spad), fz_z_mag(spad),
  //sure_fz_mag(ntype), sure_startf(ntype), lookper_startf(ntype),
  lookper_distr(spad), lookper_fz(spad),
  //lookper_nearr_dist(ntype), lookper_totalFzmag(ntype),
  nearzm(spad), near0fzm(spad), near0m(spad), since0fzm(spad), bigfzfzzm(spad), shrinkfactor(spad),
  nearmrm(spad), nearmr_f(spad), nearmr_fz(spad), extangle(spad), nearmr(*this)
{
}

/*
what's the problem here?
template<int IIW_OFFSET>
MandelPoint<worker_multi, IIW_OFFSET>::MandelPoint(MandelPointStore *store, MandelMath::worker_multi::Allocator<worker_multi> *allocator):
  self_allocator(allocator, LEN),
  store(store), f(&self_allocator), fc_c(&self_allocator), fz_r(&self_allocator), fz_c_mag(&self_allocator),
  sure_fz_mag(&self_allocator), sure_startf(&self_allocator), lookper_startf_(&self_allocator), lookper_nearr(&self_allocator),
  lookper_nearr_dist(&self_allocator), lookper_totalFzmag(&self_allocator),
  near0f(&self_allocator), root(&self_allocator)
{
  //reset();
  //should be overwritten before read:
  working_assert(self_allocator.checkFill());
}
*/

template <typename BASE>
void JuliaPoint<BASE>::readFrom(void *storage, int index)
{
  //store->assign(src.store);
  f.readFrom(storage, index+iiw_f_);
  fz_z.readFrom(storage, index+iiw_fz_z);
  fzz_z.readFrom(storage, index+iiw_fzz_z);
  fz_z_mag.readFrom(storage, index+iiw_fz_z_mag);
  lookper_distr.readFrom(storage, index+iiw_lookper_distr);
  lookper_fz.readFrom(storage, index+iiw_lookper_fz);
  near0m.readFrom(storage, index+iiw_near0m);
  nearzm.readFrom(storage, index+iiw_nearzm);
  near0fzm.readFrom(storage, index+iiw_near0fzm);
  since0fzm.readFrom(storage, index+iiw_since0fzm);
  bigfzfzzm.readFrom(storage, index+iiw_bigfzfzzm);
  shrinkfactor.readFrom(storage, index+iiw_shrinkfactor);
  nearmrm.readFrom(storage, index+iiw_nearmrm);
  nearmr_f.readFrom(storage, index+iiw_nearmr_f);
  nearmr_fz.readFrom(storage, index+iiw_nearmr_fz_);
  extangle.readFrom(storage, index+iiw_extangle);
}

template <typename BASE>
void JuliaPoint<BASE>::writeTo(void *storage, int index)
{
  //store->assign(src.store);
  f.writeTo(storage, index+iiw_f_);
  fz_z.writeTo(storage, index+iiw_fz_z);
  fzz_z.writeTo(storage, index+iiw_fzz_z);
  fz_z_mag.writeTo(storage, index+iiw_fz_z_mag);
  lookper_distr.writeTo(storage, index+iiw_lookper_distr);
  lookper_fz.writeTo(storage, index+iiw_lookper_fz);
  near0m.writeTo(storage, index+iiw_near0m);
  nearzm.writeTo(storage, index+iiw_nearzm);
  near0fzm.writeTo(storage, index+iiw_near0fzm);
  since0fzm.writeTo(storage, index+iiw_since0fzm);
  bigfzfzzm.writeTo(storage, index+iiw_bigfzfzzm);
  shrinkfactor.writeTo(storage, index+iiw_shrinkfactor);
  nearmrm.writeTo(storage, index+iiw_nearmrm);
  nearmr_f.writeTo(storage, index+iiw_nearmr_f);
  nearmr_fz.writeTo(storage, index+iiw_nearmr_fz_);
  extangle.writeTo(storage, index+iiw_extangle);
}

template <typename BASE>
void JuliaPoint<BASE>::zero(const MandelMath::complex<BASE> &first_z)
{
  f.assign_across(first_z);

  nearzm.assign(f.im); //don't have "tmp" here so do it the long way... for now
  nearzm.sqr();
  near0m.assign(f.re);
  near0m.sqr();
  near0m.add(nearzm);
  nearzm.zero(20);//assign(near0m);
  near0fzm.zero(1);
  since0fzm.zero(1);
  bigfzfzzm.zero(0);
  shrinkfactor.zero(1, 0);

  fz_z.zero(1, 0);
  fz_z_mag.zero(1);
  fzz_z.zero(0, 0);
  //lookper resets at first iter
  //lookper_distr
  //lookper_fz
  store->near0iter_1=1;
  store->nearziter_0=1;
  store->bigfzfzziter=0;
  store->nearmriter=0;
  store->nearmr_difficult=false;

  //near0m.assign(f.getMag_tmp(tmp));

  //nearzm_.assign(near0m_); //could start with infinity but we should eventually iterate closer than 0... I think
  //nearzm.zero(20);
  nearmrm.zero(20); // //overwrite at iter 0 (we don't have -root0 here)
  nearmr_f.assign(first_z);
  nearmr_fz.zero(1, 0);
  extangle.zero(0);

  store->wstate=JuliaPointStore::WorkState::stIdle;
  store->rstate=JuliaPointStore::ResultState::stUnknown;
  store->iter=0;
  //store->newton_iter=0;
  store->exterior.zero(-1);
  store->interior.zero();
  //store->has_fc_r=false;
  /*
    real exterior:=0
    real interior:=0
    initwinding(c)
    complex interiorComplex:=0
    int period:=0
    complex root:=0
  */
}

template<typename BASE>
void JuliaPoint<BASE>::NearMR::tap(ComputeJuliaParams<BASE> &params, MandelMath::complex<BASE> *tmpx)
{
#if 1 //c blocked by r
  tmpx->re.assign(*owner.f.dist2_tmp(params.root));
  tmpx->im.assign(*owner.f.dist2_tmp(params.c));
  bool is_close_r=tmpx->re.isle(owner.nearmrm);
  bool is_close_c=tmpx->im.isle(owner.nearmrm);
  assert(owner.store!=nullptr);
  if (!is_close_r)
  {
    if (is_close_c) //closer to c but not r
    {
      owner.store->nearmr_difficult=false;
      owner.nearmrm.assign(tmpx->im);
      owner.store->nearmriter=owner.store->iter+1;
      owner.nearmr_f.assign(owner.f);
      owner.nearmr_fz.assign(owner.fz_z);
    }
    //else f far from both, just ignore
  }
  else if (!is_close_c) //closer to r but not c
  { //settled on getting closer and closer to r
    //owner.store->nearmr_difficult=false;
    owner.nearmrm.assign(tmpx->re);
    //owner.store->nearmriter_=owner.store->iter+2;
  }
  else
  { //closer to both
    if (tmpx->re.isle(tmpx->im))
    { //closer to r than c
      /*if (owner.nearmr_prev.mulreT_tmp(&params.root0, scr)->isl0()) //came from closer to -r0
      {
        owner.store->nearmr_difficult=false;
        owner.nearmrm_.assign(tmpx->re);
        owner.store->nearmriter_=owner.store->iter+1;
        owner.nearmr_fz.assign(&owner.fz_z);
      }
      else
      {
        owner.store->nearmr_difficult=false;//true;
        owner.nearmrm_.assign(tmpx->re);
        owner.store->nearmriter_=owner.store->iter+1;
        owner.nearmr_fz.assign(&owner.fz_z);
      }*/
      owner.store->nearmr_difficult=false;//true;
      owner.nearmrm.assign(tmpx->re);
      owner.store->nearmriter=owner.store->iter+1;
      owner.nearmr_f.assign(owner.f);
      owner.nearmr_fz.assign(owner.fz_z);
    }
    else
    { //closer to c than r
      owner.store->nearmr_difficult=false;
      owner.nearmrm.assign(tmpx->im);
      owner.store->nearmriter=owner.store->iter+1;
      owner.nearmr_f.assign(owner.f);
      owner.nearmr_fz.assign(owner.fz_z);
    }
  }
  //owner.nearmr_prev.assign(&owner.f);
#else
#if 0 //-r0 blocked by +r0
  tmpx->from_pmdist(owner.f, params.root0, scr);
#elif 1 //0 blocked by +r0
  tmpx->re.assign(*owner.f.dist2_tmp(&params.root0, scr));
  tmpx->im.assign(*owner.f.getMag_tmp(scr));
#endif
  bool is_close_p=tmpx->re.isle(owner.nearmrm);
  bool is_close_m=tmpx->im.isle(owner.nearmrm);
  assert(owner.store!=nullptr);
  if (!is_close_p)
  {
    if (is_close_m) //closer to 0 but not r
    {
      owner.store->nearmr_difficult=false;
      owner.nearmrm.assign(tmpx->im);
      owner.store->nearmriter_=owner.store->iter+2;
    }
    //else f far from both, just ignore
  }
  else if (!is_close_m) //closer to r but not 0
  { //settled on getting closer and closer to r
    //owner.store->nearmr_difficult=false;
    owner.nearmrm.assign(tmpx->re);
    //owner.store->nearmriter_=owner.store->iter+2;
  }
  else
  { //closer to both
    if (tmpx->re.isle(tmpx->im))
    { //closer to r than 0
      owner.store->nearmr_difficult=false;//true;
      owner.nearmrm.assign(tmpx->re);
      //owner.store->nearmriter_=owner.store->iter+2;
    }
    else
    { //closer to 0 than r
      owner.store->nearmr_difficult=false;
      owner.nearmrm.assign(tmpx->im);
      owner.store->nearmriter_=owner.store->iter+2;
    }
  }
#endif
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &src): ShareableViewInfo((ShareableViewInfo const &)src)
{
  //view.assign(src.view);
  //c.assign(src.c);
  //nth_fz_limit.assign(src.nth_fz_limit);
  //dbgPoint();
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo const &src): QObject(),
  spad(src.spad), view(src.view), c(src.c), nth_fz_limit(src.nth_fz_limit),
  scale(src.scale), nth_fz(src.nth_fz), max_root_effort(src.max_root_effort)
{ //why do you need this?
  //doesn't "need" but calls it anyway dbgPoint();
  //dbgPoint();
}

ShareableViewInfo::ShareableViewInfo(ShareableViewInfo &&src): QObject(),
  spad(src.spad), view(std::move(src.view)), c(std::move(src.c)), nth_fz_limit(std::move(src.nth_fz_limit)),
  scale(std::move(src.scale)), nth_fz(std::move(src.nth_fz)), max_root_effort(std::move(src.max_root_effort))
{
  //view.assign(src.view);
  //c.assign(src.c);
  //nth_fz_limit.assign(src.nth_fz_limit);
  dbgPoint();
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &src)
{
  view.constructLateBecauseQtIsAwesome(src.view);
  //view.re.constructLateBecauseQtIsAwesome(src.view.re);
  //view.im.constructLateBecauseQtIsAwesome(src.view.im);
  //view.assign(src.view);//view.assign_across(src.view); //really across? I think the "across" happens later
  //c.re.constructLateBecauseQtIsAwesome(src.c.re);
  //c.im.constructLateBecauseQtIsAwesome(src.c.im);
  c.constructLateBecauseQtIsAwesome(src.c);//c.assign_across(src.c);
  nth_fz_limit.constructLateBecauseQtIsAwesome(src.nth_fz_limit);
  //nth_fz_limit.assign(src.nth_fz_limit);
  scale=src.scale;
  nth_fz=src.nth_fz;
  max_root_effort=src.max_root_effort;
  return *this;
}

ShareableViewInfo &ShareableViewInfo::operator=(ShareableViewInfo &&src)
{
  return operator=((ShareableViewInfo &)src);
}

template<typename BASE>
LaguerreStep<BASE>::LaguerreStep(MandelMath::number<BASE>::Scratchpad *spad):
  step(spad), s1(spad), s2(spad), tmp1(spad), tmp2(spad),
  laguG(spad), laguG2(spad), laguH(spad), laguX(spad), fzzf(spad)
{
}

template <typename BASE>
bool LaguerreStep<BASE>::eval(int lg2_degree, MandelMath::complex<BASE> const &f,
                                              MandelMath::complex<BASE> const &f_z,
                                              MandelMath::complex<BASE> const &f_zz)
{
  if (f.is0())
  {
    step.zero(0, 0);
    return true;
  };

  double order1;
  int maxm;
  if (lg2_degree<5)
  {
    maxm=1+lg2_degree; //in theory up to n-1 but for Mandelbrot only 1+period = 1+log2(n)
    order1=ldexp(1, -lg2_degree);
  }
  else if (lg2_degree<1024)
  {
    maxm=15;
    order1=ldexp(1, -lg2_degree);
  }
  else
  {
    maxm=15;
    order1=0;
  }

  //1/f should be fine, or we'd be at the root
  tmp1.assign(f);
  tmp1.recip();    //1/f
  laguG.assign(f_z);
  laguG.mul(tmp1); //laguG = f'/f
  fzzf.assign(f_zz);
  fzzf.mul(tmp1); //f''/f


    // laguH=fzf^2-fzzf
    // m=Round( Re(G^2*H^T)/mag(H) )
    // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
    laguG2.assign(laguG);
    laguG2.sqr();    //G^2
    laguH.assign(laguG2);
    laguH.sub(fzzf); //G^2-fzzf = f'^2/f^2-f''*f/f^2 = laguH
    //currentWorker->assign(tmp1.re_s, laguG2.re_s);
    //currentWorker->assign(tmp1.im_s, laguG2.im_s);
    int m=1;
    {
      /*double G2HT_re=currentWorker->toDouble(laguG2.mulreT(&laguH));
      double H_mag=currentWorker->toDouble(laguH.getMagTmp());
      //turns out that if mu=m then mu=m=G^2/H
      //1.5*mag(H)>Re(G^2*H^T) ... m=1
      //300*mag(H)<Re(G^2*H^T) ... m=300
      //mag(H)<Re(G^2*H^T)*order1 ... m=1/order1
      if (1.5*H_mag>=G2HT_re) //>= so we don't divide G^2/H if H=G=0
        m=1;
      else if ((clever.mult+0.5)*H_mag<=G2HT_re)
        m=1; //best practice is to use m=1 if H=0   clever.mult;
      else if (H_mag*(maxm-0.5)<G2HT_re)
        m=maxm;
      else
        m=qRound(G2HT_re/H_mag);*/

      //solve for m=mu:   m=G^2/H
      double mum_re=1, mum_im=0; //better than mu?
      double h_re=laguH.re.toDouble();
      double h_im=laguH.im.toDouble();
      double h_mag=h_re*h_re+h_im*h_im;
      double g2_re=laguG2.re.toDouble();
      double g2_im=laguG2.im.toDouble();
      if (h_mag>0.01)
      { //h_mag ok
        mum_re=(g2_re*h_re+g2_im*h_im)/h_mag;
        mum_im=(g2_im*h_re-g2_re*h_im)/h_mag;
      };
      dbg.mum_re=mum_re;
      dbg.mum_im=mum_im;
    }

    //m= some func of mu where mu is solution of ((1-1/n)*H/G^2-1/n) mu^2 + 2*mu/n -1=0
    //with m as input:                           ((1-m/n)*H/G^2-1/n) mu^2 + m/n 2*mu -m = 0
    double G2_mag=laguG2.getMag_double();
    if (G2_mag<0.01)
    { //G2_mag bad
      m=1;
      dbg.mu_re=1;
      dbg.mu_im=0;
    }
    else
    {
      laguX.assign(laguG2);
      laguX.im.chs();
      laguX.mul(laguH);
      double a_re=laguX.re.toDouble()/G2_mag*(1-order1)-order1;
      double a_im=laguX.im.toDouble()/G2_mag*(1-order1);
      double mu_re, mu_im;
      MandelMath::complex_double_quadratic(&mu_re, &mu_im, a_re, a_im, order1, 0, -1, 0);
      dbg.mu_re=mu_re;
      dbg.mu_im=mu_im;
      if (!(mu_re>=1.3)) //also m=1 if mu_re is NaN    (mu_re<1.3)
        m=1;
      else {/*if (abs(mu_im)>mu_re/2)
        m=1;
      else
      {
        double mu_mag=mu_re*mu_re+mu_im*mu_im;
        m=qRound(sqrt(mu_mag)); //or just round mu_re?
        */
        m=qRound(mu_re);
        if (m>maxm)
          m=maxm;
      }
    }
#if 0
    if (newtonCycle==0)
    {
      //Fejer bound: smaller solution x of
      //fzz/(n-1) x^2+2 fz x + n f=0
      //x=y*n
      //fzz*n/(n-1) y^2+2 fz y + f=0

      double r_re=currentWorker->toDouble(r->re_s);
      double r_im=currentWorker->toDouble(r->im_s);
      //numbers are small but don't need precision so let's do it in double
      double a_re=currentWorker->toDouble(fzz_r.re_s)/(1-order1);
      double a_im=currentWorker->toDouble(fzz_r.im_s)/(1-order1);
      double fz_re=currentWorker->toDouble(fz_r.re_s);
      double fz_im=currentWorker->toDouble(fz_r.im_s);
      double f_re=currentWorker->toDouble(f_r.re_s);
      double f_im=currentWorker->toDouble(f_r.im_s);
      MandelMath::complex_double_quadratic(
            &newtres.first_fejer_re, &newtres.first_fejer_im,
            a_re, a_im,
            fz_re, fz_im,
            f_re, f_im);
      newtres.first_fejer_re=r_re+ldexp(newtres.first_fejer_re, period);
      newtres.first_fejer_im=r_im+ldexp(newtres.first_fejer_im, period);

      //Batra's bound https://www.tuhh.de/ti3/paper/rump/Ru03c.pdf theorem 3.8
        //but only for real coefficients
      //|fz r|-|f + fzz/2 r^2|=0, find r
      //sqrt(fz fz^T) r=sqrt((f + fzz/2 r^2)(f^T + fzz^T/2 r^2))
      //sqrt(fz fz^T) r=sqrt((|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4))
      //(a+bi)(c-di)+(a-bi)(c+di)=2ac+2bd=2 Re(f fzz^T)
      //|fz|^2 r^2=|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4
      //0=|f|^2+ (Re(f^T fzz)-|fz|^2) rr + |fzz|^2/4 rr^2    r=sqrt(rr)

      /*MandelMath::complex_double_quadratic(&newtres.first_batra, &a_im,
          currentWorker->toDouble(fzz_r.getMagTmp())/4, 0,
          (currentWorker->toDouble(f_r.mulreT(&fzz_r))-currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
          currentWorker->toDouble(f_r.getMagTmp()), 0);
      if (newtres.first_batra>=0)
        newtres.first_batra=sqrt(newtres.first_batra);*/

      //https://ur.booksc.eu/book/5736333/a5b588
      //ZAMM - Journal of Applied Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik und Mechanik
      //1988 Vol. 68; Iss. 6
      //Dr. A. Neumaier: An Existence Test for Root Clusters and Multiple Roots
      //fi(c, r, alpha)=r abs(re((f(c+r e^ialpha)-f(c))/(c+r e^ialpha)))-abs(f(c))
      //  addition from https://ur.booksc.eu/book/5736333/a5b588 remark 3:
      //  f needs to be divided (or rotated) by f' first to make f' real
      //for all alpha, which r makes fi==0 ?
      //abs(re(f'*r+f''/2 r^2 e^ialpha))=abs(f)
      //for max re(f'*r+f''/2 r^2 e^ialpha), we need max re(f'+f''/2 r e^ialpha) because r is real
      //f'' e^ialpha=real
      //e^ialpha=f''^T/sqrt(f'' f''^T)=sqrt(f''^T/f'')
      //abs(re(f'*r+ r^2/2 sqrt(f'' f''^T)))-abs(f)=0
      //r*abs(re(f'))+ r^2/2 sqrt(f'' f''^T)-abs(f)=0
      /*if (currentWorker->isle0(fz_r.re_s))
        MandelMath::complex_double_quadratic(&newtres.first_batra, &a_im,
            +sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
            //currentWorker->toDouble(fz_r.re_s)/2, 0,
            sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
            +sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);
      else*/
      MandelMath::complex_double_quadratic(&newtres.first_neumaier1_re_, &newtres.first_neumaier1_im_,
          -sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
          sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
          -sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);

      /* Neumaier for k=2:
      Re(f''(z)/2) > |f| r^-2 + |f'| r^-1      r real, r>0
      f''(z)=f''+(z-z0)f'''+...
      |f''+r|f'''|+...|/2 r^2 > |f| + |f'| r
      ...+|f'''|r^3/2+|f''| r^2/2 > |f| + |f'|r
                      |f''| r^2/2 > |f| + |f'|r
      gives always r<0 but that's the wrong root
      |f''| r^2/2 - |f'|r - |f| =0
      2*|f'|/|f''|+-sqrt(4*|f'|^2/|f''|^2-8*|f|/|f''|)
      |f'|/|f''|+-sqrt(|f'|^2/|f''|^2+2*|f|/|f''|)
      works if |f'|^2/|f''|+2*|f|>0 i.e. always
      but r1 always<0, r2>2*|f'|/|f''|
      r2=(|f'|+sqrt(|f'|^2+2*|f|*|f''|))/|f''|

      test x(x-1) at 2+i
      f=1+3i f'=2x-1=3+2i f''=2
      r2=(|3+2*I|+sqrt(|3+2*I|^2+2*2*|1+3*I|))/2
      4.33502318885498454
      correct is 2.236
      */
      double fm=sqrt(currentWorker->toDouble(f_r.getMagTmp()));
      double fzm=sqrt(currentWorker->toDouble(fz_r.getMagTmp()));
      double fzzm=sqrt(currentWorker->toDouble(fzz_r.getMagTmp()));
      newtres.first_neumaier2_re=(fzm + sqrt(fzm*fzm+2*fm*fzzm))/fzzm;
      newtres.first_neumaier2_im=0;

      /* naive: approximate f with c(x-a)^m
      m=f'^2/(f'^2-f f'') = f'^2/f^2/(f'^2/f^2-f''/f)=G^2/H
      x-a=m/(f'/f)=m/G=G/H    looks good if |m_im|<|m_re|
      m*(x-a)=G^3/H^2

      trouble: singularities when f f''=f'^2 -> m=infinity, iteration jumps too far
                                  f'=0 -> m=0, m/(f/f') jumps too little
      */
      /*double g_re=currentWorker->toDouble(laguG.re_s);
      double g_im=currentWorker->toDouble(laguG.im_s);
      double g_mag=g_re*g_re+g_im*g_im;
      if (1e6*H_mag<=g_mag*g_mag)
      {
        newtres.first_naive_re=currentWorker->toDouble(r->re_s);
        newtres.first_naive_im=currentWorker->toDouble(r->im_s);
      }
      else
      {
        double g2_re=currentWorker->toDouble(laguG2.re_s);
        double g2_im=currentWorker->toDouble(laguG2.im_s);
        double h_re=currentWorker->toDouble(laguH.re_s);
        double h_im=currentWorker->toDouble(laguH.im_s);

        double m_re=(g2_re*h_re+g2_im*h_im)/H_mag;
        double m_im=(g2_im*h_re-g2_re*h_im)/H_mag;
        //couldn't find smooth function that:
        //1->1 2->2 3->3... 0->1 -1->1 i->1 -i->1
        //esp. since we need to have 1->1 exact and in neigborhood too
        if ((m_re<abs(m_im)*2))
        //if ((m_re<0.9) || (m_re<abs(m_im)*2)) //for m~0, we need something like sqrt(m): m is too small, 1 is too large
        {
          m_re=1;
          m_im=0;
        };
        newtres.first_naive_re=currentWorker->toDouble(r->re_s)-(m_re*g_re+m_im*g_im)/g_mag;
        newtres.first_naive_im=currentWorker->toDouble(r->im_s)-(m_im*g_re-m_re*g_im)/g_mag;
      }*/

      /* even naiver: show the 2 roots of c(x-a)(x-b) that have the same f, f', f''
      w.l.o.g. x=0
      c(x^2-(a+b)x+ab)=f''x^2/2+f'x+f
      cx^2-c(a+b)x+cab=f''x^2/2+f'x+f
      -f'/f''+-sqrt(f'^2/f''^2-2*f/f'')

      if x1 close to x2 (relative to x), use (x1+x2)/2 else use x1
      at |x1|=|x2|, 90 degrees..mult~2, use (x1+x2)/2
      at |x1|=|x2|, 60 degrees..mult~1, use x1
      at |x1|=0.8|x2|, 80% weight from x1
      at |x1|=0.5|x2|, 90% weight from x1
      at |x1|=0.3|x2|, use x1
      when x1~x2, correct guess is actually around 0.7 x1
      */
      a_re=currentWorker->toDouble(fzz_r.re_s)/2;
      a_im=currentWorker->toDouble(fzz_r.im_s)/2;
      MandelMath::complex_double_quadratic2(&newtres.first_naive1_re_, &newtres.first_naive1_im,
                                            &newtres.first_naive2_re, &newtres.first_naive2_im,
                                            a_re, a_im, fz_re/2, fz_im/2, f_re, f_im);
      double n2_rmag=1/(newtres.first_naive2_re*newtres.first_naive2_re+newtres.first_naive2_im*newtres.first_naive2_im);
      //d=naive1/naive2
      double d_re=(newtres.first_naive1_re_*newtres.first_naive2_re+newtres.first_naive1_im*newtres.first_naive2_im)*n2_rmag;
      double d_im=(newtres.first_naive1_im*newtres.first_naive2_re-newtres.first_naive1_re_*newtres.first_naive2_im)*n2_rmag;
      double d_mag=(d_re*d_re+d_im*d_im);
      double w1=1, w2=0;
      if (d_re<-0.5) //angle>120deg, even if close in magnitude
      { w1=1; w2=0; newtres.naiveChoice=NewtonNaiveChoice::ncWide; }
      else if (d_mag<0.3*0.3)
      { w1=1; w2=0; newtres.naiveChoice=NewtonNaiveChoice::nc03; }
      else if (d_mag<0.5*0.5)
      { w1=0.9; w2=0.1; newtres.naiveChoice=NewtonNaiveChoice::nc05; }
      else if (d_mag<0.8*0.8)
      { w1=0.8; w2=0.2; newtres.naiveChoice=NewtonNaiveChoice::nc08; } //or just 1;0
      else if (d_re<-0.1)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        w1=1; w2=0.0;
        newtres.naiveChoice=NewtonNaiveChoice::nc100;
      }
      else if (d_re<0)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        //don't trust M here
        w1=1; w2=0.0;
        newtres.naiveChoice=NewtonNaiveChoice::nc90_;
      }
      else if (d_re<0.1)
      {
        //most problematic case, some times best is [1,0] other times [1,1]
        //don't trust M here
        w1=1; w2=0.0;
        newtres.naiveChoice=NewtonNaiveChoice::nc80;
      }
      else if (d_re<0.5)
      {
        //can (try) use M here
        w1=1; w2=0;
        newtres.naiveChoice=NewtonNaiveChoice::nc60;
      }
      else
      {
        //can (try) use M here
        w1=newtres.firstMum_re_-1; w2=0;
        newtres.naiveChoice=NewtonNaiveChoice::ncClose;
      }
      newtres.first_naive_re=w1*newtres.first_naive1_re_+w2*newtres.first_naive2_re;
      newtres.first_naive_im=w1*newtres.first_naive1_im+w2*newtres.first_naive2_im;
      newtres.first_naive1_re_=r_re+newtres.first_naive1_re_;
      newtres.first_naive1_im=r_im+newtres.first_naive1_im;
      newtres.first_naive2_re=r_re+newtres.first_naive2_re;
      newtres.first_naive2_im=r_im+newtres.first_naive2_im;
      newtres.first_naive_re=r_re+newtres.first_naive_re;
      newtres.first_naive_im=r_im+newtres.first_naive_im;

      //Laguerre is the solution of
      //   c=-n  b=f'/f  a=f'^2/f^2*(1-n/m+1/m)-f''/f*(1-n/m)=H*(1-n/m)+G^2/m
      //   G=f'/f   H=G^2-f''/f
      //>> a=H*(1-m/n)-G^2/n  b=m*G/n  c=-m    ok
      /*
      a_re=currentWorker->toDouble(laguH.re_s)*(1-m*order1)-currentWorker->toDouble(laguG2.re_s)*order1;
      a_im=currentWorker->toDouble(laguH.im_s)*(1-m*order1)-currentWorker->toDouble(laguG2.im_s)*order1;
      double b_re=currentWorker->toDouble(laguG.re_s)*m*order1;
      double b_im=currentWorker->toDouble(laguG.im_s)*m*order1;
      MandelMath::complex_double_quadratic(
            &newtres.first_lagum_re, &newtres.first_lagum_im,
            a_re, a_im,
            b_re, b_im,
            -m, 0);
      newtres.first_lagum_re=currentWorker->toDouble(r->re_s)-newtres.first_lagum_re;
      newtres.first_lagum_im=currentWorker->toDouble(r->im_s)-newtres.first_lagum_im;
      */
      a_re=currentWorker->toDouble(laguH.re_s)*(1-order1)-currentWorker->toDouble(laguG2.re_s)*order1;
      a_im=currentWorker->toDouble(laguH.im_s)*(1-order1)-currentWorker->toDouble(laguG2.im_s)*order1;
      double b_re=currentWorker->toDouble(laguG.re_s)*order1;
      double b_im=currentWorker->toDouble(laguG.im_s)*order1;
      MandelMath::complex_double_quadratic2(
            &newtres.first_lagu1_re, &newtres.first_lagu1_im,
            &newtres.first_lagu1o_re, &newtres.first_lagu1o_im,
            a_re, a_im,
            b_re, b_im,
            -1, 0);
      newtres.first_lagu1_re=r_re-newtres.first_lagu1_re;
      newtres.first_lagu1_im=r_im-newtres.first_lagu1_im;
      newtres.first_lagu1o_re=r_re-newtres.first_lagu1o_re;
      newtres.first_lagu1o_im=r_im-newtres.first_lagu1o_im;
    };
#endif
  dbg.lastm=m;
  bool lagu_valid=false;
  bool newt_valid=false;
  if (order1>=0)
  {
    // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
    // all but last few cycles can be done just in double precision
    //   but the cost of this compared to evaluation of f,f',f'' is negligible
    laguX.assign(laguG2);
    laguX.lshift(-lg2_degree); //G^2/n
    laguX.rsub(laguH); //H-G^2/n
    tmp2.zero(m);
    tmp2.recip();
    tmp2.add_double(-order1); //1/m-1/n
    laguX.mul(tmp2); //(1/m-1/n)*(H-G^2/n)
    laguX.sqrt();
    //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
    if (laguX.mulreT_tmp(laguG)->isl0())
    {
      laguX.chs();
    };
    laguG.lshift(-lg2_degree); //G/n
    laguX.add(laguG);
    //if 1/n~0: a=1/(0 +- sqrt( (1/m)*(H) )), m can still be 1..max
    //   fine if H!=0:       a=1/( sqrt( (1/m)*(H) )), m can still be 1..max
    //   if H==0: 1/G/(1/n + sqrt( (1/300-1/n)*(-1/n) ))=1/G* -i*sqrt(300*n)
    //   if H=G=0: 1/0
    //if G=0: a=1/(+- sqrt( (1/m-1/n)*(H) ))     m=1
    //   fine if H!=0: a=+-(sqrt(n/(n-1))*sqrt(f/-f''))       x^2+9 at 0: f=9 f''=2 -> +-3i
    //   if H=0: a=1/0
    //if H=0: a=1/G*m*(1 - i*sqrt(n/m-1))  m~n -> a=n/G;  m~300 -> a=-i/G*sqrt(n*300)
    //        a=1/G*m*n*(1/n - i*sqrt(1/m/n-1/n^2))
    double X_mag=laguX.getMag_tmp()->toDouble();
    if (X_mag>=1e-60)
    {
      laguX.recip_prepared();
      lagu_valid=true;
    };
    //else
    //we should move the guess a little and try again, but
    //  we can leave this to the caller
    //return 0;
  };
  if (!f_z.is0()) //gz_r_mag!=0)
  {
    //newton near multiroot:
    //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
    //Newton corrected for multiroot = f/f'*m
    //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
    step.assign(f_z);
    step.recip();
    step.mul(f); //f/f'
    if (m!=1)
    {
      step.mul_double(m);
    };
    newt_valid=true;
  };
#if 0
  if (newtonCycle==0)
  {
    currentWorker->assign(&newtres.first_guess_newt_re, r->re_s);
    currentWorker->assign(&newtres.first_guess_newt_im, r->im_s);
    if (newt_valid)
    {
      currentWorker->sub(&newtres.first_guess_newt_re, newtX.re_s);
      currentWorker->sub(&newtres.first_guess_newt_im, newtX.im_s);
    };

    currentWorker->assign(&newtres.first_guess_lagu_re, r->re_s);
    currentWorker->assign(&newtres.first_guess_lagu_im, r->im_s);
    if (lagu_valid)
    {
      currentWorker->sub(&newtres.first_guess_lagu_re, laguX.re_s);
      currentWorker->sub(&newtres.first_guess_lagu_im, laguX.im_s);
    };
  };
#endif
  if (!newt_valid)
  {
    if (!lagu_valid)
    {
      return false;
    };
    step.assign(laguX);
  }
  else if (!lagu_valid)
  {
    //keep newtX
  }
  else
  {
    if (m>1)//(fastHoming && (newtonCycle<2) && (m>1))
    {
      step.assign(laguX);
    }
    else
    {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
      double N_mag=step.getMag_double();
      double L_mag=laguX.getMag_double();
      if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
      {
        step.assign(laguX);
      };
    }
  }
#if 0
  if ((g_r_mag>bestfm) && (newtonCycle>30))
  {
    currentWorker->lshift(newtX.re_s, -2);
    currentWorker->lshift(newtX.im_s, -2);
  };
#endif
  //currentWorker->sub(r->re_s, newtX.re_s);
  //currentWorker->sub(r->im_s, newtX.im_s);
  return true;
}

template<typename BASE>
MandelLoopEvaluator<BASE>::MandelLoopEvaluator(MandelMath::number<BASE>::Scratchpad *spad):
  f(spad), f_z(spad), f_c(spad),
  f_zz(spad), f_zc(spad), f_cc(spad), f_zzc(spad),
  multi(0), first_multi(spad), near0iter_1(0), sumA(spad), f_z_mag(spad), f_z_minmag(spad),
  s1(spad), s2(spad)
{
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::evalg(int period, MandelMath::complex<BASE> const &c)
{
  f.assign(c);
  f_c.zero(1, 0);
  f_cc.zero(0, 0);
  //also using s2
  for (int i=0; i<period; i++)
  {
    //g_cc=2*(g_cc*g + g_c*g_c)
    f_cc.mul(f);
    s2.assign(f_c);
    s2.sqr();
    f_cc.add(s2);
    f_cc.lshift(1);
    //g_c=2*g_c*g+1
    f_c.mul(f);
    f_c.lshift(1);
    //g=g^2+c
    f.sqr();
    if (i+1!=period)
    {
      f_c.re.add_double(1);
      f.add(c);
    };
    double f_mag=f.getMag_double();
    double allmag=f_mag+
                  f_c.getMag_double()+
                  f_cc.getMag_double();
    if (allmag>1e60)
      return false;
  }
  return true;
}

/**
 * @brief find f, f_z, f_zz, f_c, f_zc, f_cc
 * @param period
 * @param c
 * @param z
 * @return false on overflow
 */
template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval2(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z, bool minusZ, bool doInit)
{
  if (doInit)
  {
    f.assign(z);
    f_z.zero(1, 0);
    f_zz.zero(0, 0);
    f_c.zero(0, 0);
    f_zc.zero(0, 0);
    f_cc.zero(0, 0);
  };
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double()+
                 f_c.getMag_double()+
                 f_zc.getMag_double()+
                 f_cc.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;

    //f^2+c -> 2f fc+1 -> 2f fcc+2fc fc
    //f_cc = 2 * (f_c^2 + f * f_cc)
    s2.assign(f_c);
    s2.sqr();
    f_cc.mul(f);
    f_cc.add(s2);
    f_cc.lshift(1);

    //f^2+c -> 2f fz -> 2fc fz+2f fzc
    //f_zc = 2 * (f_z * f_c + f * f_zc);
    s2.assign(f_c);
    s2.mul(f_z);
    f_zc.mul(f);
    f_zc.add(s2);
    f_zc.lshift(1);

    //f^2+c -> 2f fz -> 2fz fz+2f fzz
    //f_zz = 2 * (f_z^2 + f * f_zz)
    s2.assign(f_z);
    s2.sqr();
    f_zz.mul(f);
    f_zz.add(s2);
    f_zz.lshift(1);

    //f^2+c -> 2f fc+1
    //f_c = 2 * f * f_c + 1
    f_c.mul(f);
    f_c.lshift(1);
    f_c.re.add_double(1);

    //f^2+c -> 2f fz
    //f_z = 2 * f * f_z
    f_z.mul(f);
    f_z.lshift(1);
    f_z_mag.mul(*f.getMag_tmp());
    f_z_mag.lshift(2); //f_z_mag*=4*mag(f)

    //f = f^2 + c
    f.sqr();
    if (minusZ && i+1==period)
    {
      s2.assign(c);
      s2.sub(z);
      f.add(s2);
      f_z.re.add_double(-1);
    }
    else
      f.add(c);
  }
  return true;
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval2_mag(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_c.zero(0, 0);
  f_zz.zero(0, 0);
  f_zc.zero(0, 0);
  f_z_mag.zero(1);
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 //f_z.getMag_double()+
                 f_z_mag.toDouble()+
                 f_zz.getMag_double()+
                 f_c.getMag_double()+
                 f_zc.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;

    //f^2+c -> 2f fz -> 2fc fz+2f fzc
    //f_zc = 2 * (f_z * f_c + f * f_zc);
    s2.assign(f_c);
    s2.mul(f_z);
    f_zc.mul(f);
    f_zc.add(s2);
    f_zc.lshift(1);

    //f^2+c -> 2f fz -> 2fz fz+2f fzz
    //f_zz = 2 * (f_z^2 + f * f_zz)
    s2.assign(f_z);
    s2.sqr();
    f_zz.mul(f);
    f_zz.add(s2);
    f_zz.lshift(1);

    //f^2+c -> 2f fc+1
    //f_c = 2 * f * f_c + 1
    f_c.mul(f);
    f_c.lshift(1);
    f_c.re.add_double(1);

    //f^2+c -> 2f fz
    //f_z = 2 * f * f_z
    f_z.mul(f);
    f_z.lshift(1);
    f_z_mag.mul(*f.getMag_tmp());
    f_z_mag.lshift(2); //f_z_mag*=4*mag(f)

    //f = f^2 + c
    f.sqr();
    f.add(c);
  }
  return true;
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval_zz(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z, bool minusZ, bool doInit)
{
  if (doInit)
  {
    f.assign(z);
    f_z.zero(1, 0);
    f_zz.zero(0, 0);
    f_z_minmag.zero(10000);
  };
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(f);
    s2.assign(f_z);
    s2.sqr();
    f_zz.add(s2);
    f_zz.lshift(1);
    //f_z:=2*f*f_z
    f_z.mul(f);
    f_z.lshift(1);
    f_z_minmag.min(*f_z.getMag_tmp());
    //bigger than 3 is common
    //f:=f^2+c
    f.sqr();
    if (minusZ && i+1==period)
    { //f+=c-r instead of f:=f+c-r for better precision
      s2.assign(c);
      s2.sub(z);
      f.add(s2);
      f_z.re.add_double(-1);
    }
    else
      f.add(c);
  }
  return true;
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval2zzc(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_c.zero(0, 0);
  f_zz.zero(0, 0);
  f_zc.zero(0, 0);
  f_cc.zero(0, 0);
  f_zzc.zero(0, 0);
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double()+
                 f_c.getMag_double()+
                 f_zc.getMag_double()+
                 f_cc.getMag_double()+
                 f_zzc.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;
    //f_zzc=d/dc 2*(f_z*f_z + f*f_zz)=2*(2*f_z*f_zc + f_c*f_zz+f*f_zzc)
    f_zzc.mul(f);
    s2.assign(f_c);
    s2.mul(f_zz);
    f_zzc.add(s2);
    s2.assign(f_z);
    s2.mul(f_zc);
    s2.lshift(1);
    f_zzc.add(s2);
    f_zzc.lshift(1);
    //f_cc=2 * (f_c^2 + f * f_cc)
    f_cc.mul(f);
    s2.assign(f_c);
    s2.sqr();
    f_cc.add(s2);
    f_cc.lshift(1);
    // f_zc = 2 * (f_z * f_c + f * f_zc);
    f_zc.mul(f);
    s2.assign(f_c);
    s2.mul(f_z);
    f_zc.add(s2);
    f_zc.lshift(1);
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(f);
    s2.assign(f_z);
    s2.sqr();
    f_zz.add(s2);
    f_zz.lshift(1);
    // f_c = 2 * f * f_c + 1;
    f_c.mul(f);
    f_c.lshift(1);
    f_c.re.add_double(1);
    //f_z:=2*f*f_z
    f_z.mul(f);
    f_z.lshift(1);
    //f:=f^2+c
    f.sqr();
    f.add(c);
  }
  return true;
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval_multi(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z, MandelMath::complex<BASE> const &f_z_target, double dist_tolerance)
{
  f.assign(z);
  f_z.zero(1, 0);
  f_zz.zero(0, 0);
  f_z_minmag.zero(10000);
  //int near1=0;
  //int sumnear1=0;
  //bool dangerzone=false;
  int reducedbymag=period;
  double closest_accepted=10000, closest_rejected=10000;
  int sumA_cnt=0;
  for (int i=0; i<period; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(f);
    s2.assign(f_z);
    s2.sqr();
    f_zz.add(s2);
    f_zz.lshift(1);
    //f_z:=2*f*f_z
    f_z.mul(f);
    f_z.lshift(1);
    f_z_minmag.min(*f_z.getMag_tmp());
    //f:=f^2+c
    f.sqr();
    f.add(c);

    /*
    assume all nearby roots are on a circle and we visit them all during period (plus other points that we want to filter out)
    z is on the circle at angle f_z_target
    f is on the circle at angle f_z (normalized to |f_z|=1 ?)
    center of circle=A, circle at angle 0 =D
    (z-A)/(D-A)=f_z_target
    (f-A)/(D-A)=f_z
    A=? (D=?)
    z-A=D*f_z_target-A*f_z_target
    f-A=D*f_z-A*f_z
    z*f_z=D*f_z_target*f_z+A*f_z*(1-f_z_target)
    f*f_z_target=D*f_z_target*f_z+A*f_z_target*(1-f_z)
    (z*f_z-f*f_z_target)/(f_z-f_z_target)=A
    */
    if (f_z.dist2_double(f_z_target)*period>0.5) //assuming all points are spaced regularly at 1/period around the circle,
    {    //TODO: should be *multi not *period                //skip those close to target to avoid div by 0
      s1.assign(z);
      s1.mul(f_z);
      s2.assign(f);
      s2.mul(f_z_target);
      s2.rsub(s1); //z*f_z-f*f_z_target
      s1.assign(f_z);
      s1.sub(f_z_target);
      s1.recip();
      s1.mul(s2); //s1=A

      double dist_to_A=s1.dist2_double(z);
      double f_z_err=f_z.dist2_double(f_z_target);
      double f_zz_mag=f_zz.getMag_double();
      double expected=f_z_err/f_zz_mag;
      (void)expected;
      if (dist_to_A*3+dist_tolerance<closest_accepted)
      {
        closest_rejected=closest_accepted;
        closest_accepted=dist_to_A;
        sumA_cnt=1;
        sumA.assign(s1);
        reducedbymag=MandelMath::gcd(period, i+1);
        first_multi.assign(f_z);
      }
      else if (dist_to_A<3*closest_accepted+dist_tolerance)
      {
        if (dist_to_A<closest_accepted)
          closest_accepted=dist_to_A;
        reducedbymag=MandelMath::gcd(reducedbymag, i+1);
        sumA.add(s1);
        sumA_cnt++;
      }
      else if (dist_to_A<closest_rejected)
        closest_rejected=dist_to_A;
    }
    /*double f_z_mag=currentWorker->toDouble(f_z.getMagTmp());
    if (f_z_mag<0.99)
    {
      //dangerzone=true;
      reducedbymag=MandelMath::gcd(reducedbymag, i+1);
    }
    else if (f_z_mag<1.01)
    {
      near1++;
      sumnear1+=i;
      if (near1==1)
      {
        currentWorker->assign(&first_multi_re, f_z.re_s);
        currentWorker->assign(&first_multi_im, f_z.im_s);
      };
    }
    else if (f_z_mag<2.59) //period 15/5 iter 0+3k: 2.682..2.687
      dangerzone=true;     //period 9/3 iter 0: 2.5905
                           //
                           //
    //bigger than 3 is common
    */
  }
  s1.assign(f_z_target);
  //maybe sqrt(dist2) but not this  double f_z_err=f_z.dist2_double(&s1);//f_z.getMag_double()-s1.getMag_double(); //should==0 when called from periodCheck?
  double f_z_err=f_z.getMag_double()-1;//s1.getMag_double(); //should==0 when called from periodCheck?
  double f_zz_mag=f_zz.getMag_double();
  //if f_z_err==0, it can be anything between 0 and eps/2
  //TODO: if f_z_err=0 & f_zz==0 then expected ~ eps^(1/3) but we just use infinity
  //f_z ~ 1+f_zz*deltaz -> deltaz=(f_z-1)/f_zz
  double expected=(std::abs(f_z_err)+c.re.eps2()/4)/f_zz_mag;
  double expected_more=c.re.eps2()/f_z.getDist1_tmp()->toDouble();
  if (expected_more*expected_more>c.re.eps2())
    expected_more=std::sqrt(c.re.eps2());
  if (expected_more>expected/10+5*c.re.eps2())
    nop();
  expected+=expected_more;
  /* try to make it work with dist_tolerance instead
  double min_reject_ratio;
  {
    double min_fz=f_z_minmag.toDouble()/f_z_target->getMag_double();
    if (min_fz>0.99)
      min_reject_ratio=13;
    else if (min_fz>0.9)
      min_reject_ratio=10;
    else  if (min_fz>0.7)
      min_reject_ratio=3;
    else
      min_reject_ratio=13;
    //
    //
    //
    //
    //
    //
  }*/
  if (reducedbymag>=period)
    multi=1;
  else if (closest_accepted>expected*1000)
    multi=1;
  else if (closest_accepted>expected*3)
  {
    //dbgPoint();
    multi=1;
  }
  else if (closest_rejected<13*closest_accepted) //12.27 at second 76/38/19
    multi=1;
  else if (closest_rejected>100000*closest_accepted)
  {
    multi=period/reducedbymag;
    s1.re.zero(sumA_cnt);
    s1.re.recip();
    sumA.mul(s1.re);
  }
  else if (closest_rejected<100*closest_accepted)
    multi=1;
  else if (closest_rejected>1000*closest_accepted)
  {
    multi=period/reducedbymag;
    s1.re.zero(sumA_cnt);
    s1.re.recip();
    sumA.mul(s1.re);
  }
  /*else if (near1<1)// || dangerzone)
  {
    multi=0;
  }
  else if (near1==1)
  {
    multi=1;
  }
  else
  {
    //expect near1 numbers, each a multiple of period/near1, and that's 1,2,3..period/near1 multiple
    //(period/near1)*near1*(near1+1)/2=period*(near1+1)/2
    //ex: p=15 near1=3 5+10+15=30=15*4/2
    //except we add 2,5,8... not 3,6,9 so period*(near1+1)/2-near1
    if (period%near1!=0)
    { dbgPoint(); multi=0; }
    else if (sumnear1==period*(near1+1)/2-near1)
    {
      multi=near1;
    }
    else
    { dbgPoint(); multi=0; }
  }*/
  else
    multi=0;
  return true;
}

template <typename BASE>
bool MandelLoopEvaluator<BASE>::eval_near0(int period, MandelMath::complex<BASE> const &c)
{
  f.assign(c);
  f_z.zero(1, 0); //are they even needed?
  f_zz.zero(0, 0);
  near0iter_1=1;
  f_z_minmag.assign(*c.getMag_tmp());
  for (int i=0; i<period-1; i++)
  {
    double m_sum=f.getMag_double()+
                 f_z.getMag_double()+
                 f_zz.getMag_double();
    if (m_sum>MandelEvaluator<BASE>::LARGE_FLOAT2)
      return false;
    //f_zz:=2*(f_z*f_z + f*f_zz)
    f_zz.mul(f);
    s2.assign(f_z);
    s2.sqr();
    f_zz.add(s2);
    f_zz.lshift(1);
    //f_z:=2*f*f_z
    f_z.mul(f);
    f_z.lshift(1);
    //f:=f^2+c
    f.sqr();
    f.add(c);

    //if (i+2<=period)
    {
      const MandelMath::number<BASE> *fmag=f.getMag_tmp();
      if (!f_z_minmag.isle(*fmag)) //fmag<near0mag
      {
        near0iter_1=i+2;
        f_z_minmag.assign(*fmag);
      };
    };
  }
  return true;
}

template<typename BASE>
bool MandelLoopEvaluator<BASE>::eval_ext_mandMJ(double mandel, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &z, int iter)
{ //mandel=1.0 -> d/dc f_c(c) |c, mandel=0.0 -> d/dz f_c(z) |z
  f.assign(z);
  f_c.zero(1, 0);
  f_cc.zero(0, 0);
  //also using s2
  for (int i=0; i<iter; i++)
  {
    //g_cc=2*(g_cc*g + g_c*g_c)
    f_cc.mul(f);
    s2.assign(f_c);
    s2.sqr();
    f_cc.add(s2);
    f_cc.lshift(1);
    //g_c=2*g_c*g for Julia, g_c=2*g_c*g+1 for Mandelbrot
    f_c.mul(f);
    f_c.lshift(1);
    f_c.re.add_double(mandel);
    //g=g^2+c
    f.sqr();
    f.add(c);
    double f_mag=f.getMag_double();
    double allmag=f_mag+
                  f_c.getMag_double()+
                  f_cc.getMag_double();
    if (allmag>1e130) //goes pretty high
      return false;
  }
  return true;
}



MandelEvaluatorThread::MandelEvaluatorThread(MandelEvaluator<MandelMath::number_any> *owner):
  QThread(nullptr), owner(owner)
{

}

void MandelEvaluatorThread::syncMandel()
{
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      //(specific_cast<MandelEvaluator<MandelMath::number<MandelMath::number_a *>> *, MandelEvaluator<MandelMath::worker_multi> *>(owner))->syncComputeSplit();
      owner->syncMandelSplit();
      return;
    case MandelMath::NumberType::typeDouble:
      (specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncMandelSplit();
      return;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
      (specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncMandelSplit();
      return;
    case MandelMath::NumberType::typeDDouble:
      (specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncMandelSplit();
      return;
    case MandelMath::NumberType::typeQDouble:
      (specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncMandelSplit();
      return;
    case MandelMath::NumberType::typeReal642:
      (specific_cast<MandelEvaluator<MandelMath::real642> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncMandelSplit();
      return;
#endif
  }
}

void MandelEvaluatorThread::syncJulia()
{
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      //(specific_cast<MandelEvaluator<MandelMath::number<MandelMath::number_a *>> *, MandelEvaluator<MandelMath::worker_multi> *>(owner))->syncComputeSplit();
      owner->syncJuliaSplit();
      return;
    case MandelMath::NumberType::typeDouble:
      (specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncJuliaSplit();
      return;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
      (specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncJuliaSplit();
      return;
    case MandelMath::NumberType::typeDDouble:
      (specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncJuliaSplit();
      return;
    case MandelMath::NumberType::typeQDouble:
      (specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncJuliaSplit();
      return;
    case MandelMath::NumberType::typeReal642:
      (specific_cast<MandelEvaluator<MandelMath::real642> *, MandelEvaluator<MandelMath::number_any> *>(owner))->syncJuliaSplit();
      return;
#endif
  }
}

int MandelEvaluatorThread::syncLaguerre(int period, MandelMath::complex<MandelMath::number_any> *r, const bool fastHoming)
{
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      return owner->laguerre(period, owner->currentParams.mandel.c, r, fastHoming);
    case MandelMath::NumberType::typeDouble:
    {
      MandelEvaluator<double> *owner_access=specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->laguerreData.r.assign_across(*r);
      return owner_access->laguerre(period, owner_access->currentParams.mandel.c, &owner_access->laguerreData.r, fastHoming);
    } break;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
    {
      MandelEvaluator<__float128> *owner_access=specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->laguerreData.r.assign_across(*r);
      return owner_access->laguerre(period, owner_access->currentParams.mandel.c, &owner_access->laguerreData.r, fastHoming);
    } break;
    case MandelMath::NumberType::typeDDouble:
    {
      MandelEvaluator<MandelMath::dd_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->laguerreData.r.assign_across(*r);
      return owner_access->laguerre(period, owner_access->currentParams.mandel.c, &owner_access->laguerreData.r, fastHoming);
    } break;
    case MandelMath::NumberType::typeQDouble:
    {
      MandelEvaluator<MandelMath::dq_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->laguerreData.r.assign_across(*r);
      return owner_access->laguerre(period, owner_access->currentParams.mandel.c, &owner_access->laguerreData.r, fastHoming);
    } break;
    case MandelMath::NumberType::typeReal642:
    {
      MandelEvaluator<MandelMath::dq_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->laguerreData.r.assign_across(*r);
      return owner_access->laguerre(period, owner_access->currentParams.mandel.c, &owner_access->laguerreData.r, fastHoming);
    } break;
#endif
  }
  return 0;
}

void MandelEvaluatorThread::doMandelThreaded(int epoch)
{
  //switch (((MandelEvaluator<worker_multi> *)(this))->currentWorker->ntype())
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      owner->doMandelThreadedSplit(epoch);
      return;
    case MandelMath::NumberType::typeDouble:
    {
      MandelEvaluator<double> *owner_access=specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doMandelThreadedSplit(epoch);
    } return;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
    {
      MandelEvaluator<__float128> *owner_access=specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doMandelThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeDDouble:
    {
      MandelEvaluator<MandelMath::dd_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doMandelThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeQDouble:
    {
      MandelEvaluator<MandelMath::dq_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doMandelThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeReal642:
    {
      MandelEvaluator<MandelMath::real642> *owner_access=specific_cast<MandelEvaluator<MandelMath::real642> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doMandelThreadedSplit(epoch);
    } return;
#endif
  }
}

void MandelEvaluatorThread::doJuliaThreaded(int epoch)
{
  //switch (((MandelEvaluator<worker_multi> *)(this))->currentWorker->ntype())
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      owner->doJuliaThreadedSplit(epoch);
      return;
    case MandelMath::NumberType::typeDouble:
    {
      MandelEvaluator<double> *owner_access=specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doJuliaThreadedSplit(epoch);
    } return;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
    {
      MandelEvaluator<__float128> *owner_access=specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doJuliaThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeDDouble:
    {
      MandelEvaluator<MandelMath::dd_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doJuliaThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeQDouble:
    {
      MandelEvaluator<MandelMath::dq_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doJuliaThreadedSplit(epoch);
    } return;
    case MandelMath::NumberType::typeReal642:
    {
      MandelEvaluator<MandelMath::real642> *owner_access=specific_cast<MandelEvaluator<MandelMath::real642> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doJuliaThreadedSplit(epoch);
    } return;
#endif
  }
}

void MandelEvaluatorThread::doLaguerreThreaded(int epoch, int laguerrePeriod)
{
  //switch (((MandelEvaluator<worker_multi> *)(this))->currentWorker->ntype())
  switch (owner->ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      owner->doLaguerreThreadedSplit(epoch, laguerrePeriod);
      return;
    case MandelMath::NumberType::typeDouble:
    {
      MandelEvaluator<double> *owner_access=specific_cast<MandelEvaluator<double> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doLaguerreThreadedSplit(epoch, laguerrePeriod);
    } return;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
    {
      MandelEvaluator<__float128> *owner_access=specific_cast<MandelEvaluator<__float128> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doLaguerreThreadedSplit(epoch, laguerrePeriod);
    } return;
    case MandelMath::NumberType::typeDDouble:
    {
      MandelEvaluator<MandelMath::dd_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dd_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doLaguerreThreadedSplit(epoch, laguerrePeriod);
    } return;
    case MandelMath::NumberType::typeQDouble:
    {
      MandelEvaluator<MandelMath::dq_real> *owner_access=specific_cast<MandelEvaluator<MandelMath::dq_real> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doLaguerreThreadedSplit(epoch, laguerrePeriod);
    } return;
    case MandelMath::NumberType::typeReal642:
    {
      MandelEvaluator<MandelMath::real642> *owner_access=specific_cast<MandelEvaluator<MandelMath::real642> *, MandelEvaluator<MandelMath::number_any> *>(owner);
      owner_access->doLaguerreThreadedSplit(epoch, laguerrePeriod);
    } return;
#endif
  }
}


template <typename BASE>
MandelEvaluator<BASE>::MandelEvaluator(MandelMath::NumberType ntype, bool dontRun):
  busyEpoch(0), thread((MandelEvaluator<MandelMath::number_any> *)this),
  ntype(MandelMath::NumberTypeFromBase<BASE>::ntype),
  totalNewtonIterations(0),
  tmp(ntype),
  currentParams(&tmp),
  laguerreData(&laguerreStore, &tmp),
  mandelData(&mandelDataStore, &tmp),
  juliaData(&juliaStore, &tmp),
  loope(&tmp),
  newtres(&tmp),
  eval(&tmp),
  newt(&tmp),
  interior(&tmp),
  bulb(&tmp, &loope),
  extangle(&tmp, &loope, &bulb.lagu, this)
{
  if (!dontRun)
    thread.start(QThread::Priority::LowestPriority);
  //working_assert(currentWorker->getAllocator()->checkFill());
  workIfEpoch=-1;
  pointsComputed=0;
  timeOuterTotal=0;
  timeInnerTotal=0;
  timeInvokePostTotal=0;
  timeInvokeSwitchTotal=0;
  timeThreadedTotal=0;
  if (!dontRun)
    thread.moveToThread(&thread);
}

template <class WORKER_MULTI>
void MandelEvaluator<WORKER_MULTI>::startRunning()
{
  thread.start(QThread::Priority::LowestPriority);
  thread.moveToThread(&thread);
}

template <class WORKER_MULTI>
MandelEvaluator<WORKER_MULTI>::~MandelEvaluator()
{
}

#if 0
#if NUMBER_DOUBLE_EXISTS
void MandelEvaluator::simple_double(double cr, double ci, MandelPoint &data, int maxiter)
{
  double zr=0;
  double zi=0;
  for (int iter=0; iter<maxiter; iter++)
  {
    if (zr*zr+zi*zi>4)
    {
      data.store->rstate=MandelPointStore::ResultState::stOutside;
      data.store->iter=iter;
      data.f.zero(zr, zi);
      return;
    };
    double tmp=zr*zr-zi*zi+cr;
    zi=2*zr*zi+ci;
    zr=tmp;
  }
  //data.state=MandelPoint::State::stMaxIter;
  data.store->iter=maxiter;
  data.f.zero(zr, zi);
}
#endif //NUMBER_DOUBLE_EXISTS

#if 0
maybe with local worker
void MandelEvaluator::simple_ddouble(MandelMath::dd_real *cr, MandelMath::dd_real *ci, MandelPoint &data, int maxiter)
{
  MandelMath::dd_real zr;
  MandelMath::dd_real zi;
  MandelMath::dd_real r2;
  MandelMath::dd_real i2;
  MandelMath::dd_real t;
  for (int iter=0; iter<maxiter; iter++)
  {
    r2.assign(zr); r2.sqr();
    i2.assign(zi); i2.sqr();
    t.assign(r2); t.add(i2.hi, i2.lo_);
    if (t.hi>4)
    {
      data.state=MandelPoint::State::stOutside;
      data.iter=iter;
      data.f_re.as.ddouble_.dd->assign(zr);
      data.f_im.as.ddouble_.dd->assign(zi);
      return;
    };
    t.assign(r2); t.add(-i2.hi, -i2.lo_); t.add(cr->hi, cr->lo_); //double tmp=zr*zr-zi*zi+cr;
    zi.mul(2*zr.hi, 2*zr.lo_); zi.add(ci->hi, ci->lo_);
    zr.assign(t);
  }
  data.iter=maxiter;
  data.f_re.as.ddouble_.dd->assign(zr);
  data.f_im.as.ddouble_.dd->assign(zi);
}

void MandelEvaluator::simple_multi(MandelMath::multiprec *cr, MandelMath::multiprec *ci, MandelPoint &data, int maxiter)
{
  (void)cr;
  (void)ci;
  data.state=MandelPoint::State::stMaxIter;
  data.iter=maxiter;
}
#endif
#endif

template <class WORKER_MULTI>
void MandelEvaluator<WORKER_MULTI>::syncMandelSplit()
{
  //still called from paintOrbit()
  evaluateMandel();
  pointsComputed++;
}

template <class WORKER_MULTI>
void MandelEvaluator<WORKER_MULTI>::syncJuliaSplit()
{
  //still called from paintOrbit()
  evaluateJulia();
  pointsComputed++;
}

#if 0
void MandelEvaluator::doCompute()
{
  dbgPoint();
  timeInvokeSwitchTotal+=timeInvoke.nsecsElapsed();
  timeInner.start();
  //simple_double(currentParams.cr_n.impl->store->as.doubl, currentParams.ci_n.impl->store->as.doubl, currentData, currentParams.maxiter);
  evaluate();
  pointsComputed++;
  //msleep(10);
  timeInnerTotal+=timeInner.nsecsElapsed();
  emit doneCompute(this);
}
#endif

template <typename BASE>
void MandelEvaluator<BASE>::doMandelThreadedSplit(int epoch)
{
  busyEpoch=epoch;
  timeThreaded.start();
  threaded_errorcode=0;
#if THREADED_DONE_GIVE_WORK
  if (threaded.give(this)) //debug && !wantStop)
#endif
  {
    for (;;)
    {
      timeInvoke.start();
#if THREADED_DONE_GIVE_WORK
#else
      threaded_errorcode=threaded.give(this); //debug && !wantStop)
      if (threaded_errorcode)
        break;
#endif
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      timeInner.start();
      evaluateMandel();
      timeInnerTotal+=timeInner.nsecsElapsed();
      timeOuter.start();
      threaded_errorcode=threaded.doneMandel(this, THREADED_DONE_GIVE_WORK);
      if (threaded_errorcode)
        break;
      timeOuterTotal+=timeOuter.nsecsElapsed();
    }
    /*while (threaded.give(this))//debug && !wantStop)
    {
      evaluate();
      if (!threaded.done(this))
        break;
    }*/
  }
  timeThreadedTotal+=timeThreaded.nsecsElapsed();
  emit thread.doneMandelThreaded(&thread);
}

template <typename BASE>
void MandelEvaluator<BASE>::doJuliaThreadedSplit(int epoch)
{
  busyEpoch=epoch;
  timeThreaded.start();
  threaded_errorcode=0;
#if THREADED_DONE_GIVE_WORK
  if (threaded.give(this)) //debug && !wantStop)
#endif
  {
    for (;;)
    {
      timeInvoke.start();
#if THREADED_DONE_GIVE_WORK
#else
      threaded_errorcode=threaded.give(this); //debug && !wantStop)
      if (threaded_errorcode)
        break;
#endif
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      timeInner.start();
      evaluateJulia();
      timeInnerTotal+=timeInner.nsecsElapsed();
      timeOuter.start();
      threaded_errorcode=threaded.doneJulia(this, THREADED_DONE_GIVE_WORK);
      if (threaded_errorcode)
        break;
      timeOuterTotal+=timeOuter.nsecsElapsed();
    }
    /*while (threaded.give(this))//debug && !wantStop)
    {
      evaluate();
      if (!threaded.done(this))
        break;
    }*/
  }
  timeThreadedTotal+=timeThreaded.nsecsElapsed();
  emit thread.doneJuliaThreaded(&thread);
}

template <typename BASE>
void MandelEvaluator<BASE>::doLaguerreThreadedSplit(int epoch, int laguerrePeriod)
{
  busyEpoch=epoch;
  timeThreaded.start();
  threaded_errorcode=0;
  MandelMath::complex<BASE> &c=currentParams.mandel.c;
  MandelMath::complex<BASE> &root=laguerreData.r;
#if THREADED_DONE_GIVE_WORK
  if (threaded.give(this)) //debug && !wantStop)
#endif
  {
    for (;;)
    {
      timeInvoke.start();
#if THREADED_DONE_GIVE_WORK
#else
      threaded_errorcode=threaded.give(this); //debug && !wantStop)
      if (threaded_errorcode)
        break;
#endif
      timeInvokePostTotal+=timeInvoke.nsecsElapsed();
      timeInner.start();
      int result=laguerre(laguerrePeriod, c, &root, true);
      timeInnerTotal+=timeInner.nsecsElapsed();
      timeOuter.start();
      threaded_errorcode=threaded.doneLaguerre(this, result, THREADED_DONE_GIVE_WORK);
      if (threaded_errorcode)
        break;
      timeOuterTotal+=timeOuter.nsecsElapsed();
    }
    /*while (threaded.give(this))//debug && !wantStop)
    {
      evaluate();
      if (!threaded.done(this))
        break;
    }*/
  }
  timeThreadedTotal+=timeThreaded.nsecsElapsed();
  emit thread.doneLaguerreThreaded(&thread);
}
/*
void MandelEvaluator::startNewton(int period, const MandelMath::complex *c)
{
  currentData.store->period=period;
  currentParams.c.assign(c);
  QMetaObject::invokeMethod(this, &MandelEvaluator::doNewton, Qt::ConnectionType::QueuedConnection);
}

void MandelEvaluator::doNewton()
{
  //MandelMath::complex tmpc(currentWorker, &currentParams.c_re, &currentParams.c_im, true);
  //MandelMath::complex root(currentWorker, &currentData.f_re, &currentData.f_im, true);
  int result=newton(currentData.store->period, &currentParams.c, &currentData.root, true);
  emit doneNewton(this, result);
}
*/

template <typename BASE>
void MandelEvaluator<BASE>::Bulb::findBulbBase(int period2, MandelMath::complex<BASE> const &c/*,
    MandelMath::complex<WORKER_MULTI> *cb, MandelMath::complex<WORKER_MULTI> *rb, MandelMath::complex<WORKER_MULTI> *xc,
    MandelMath::complex<WORKER_MULTI> *baseZC, MandelMath::complex<WORKER_MULTI> *baseCC, bool *is_card, int *foundMult,
    MandelMath::complex<WORKER_MULTI> *baseFz*/)
//on input foundMult=0 -> guess rb here; =1 -> rb already set
//xc: bulb center, both z and c
//cb: bulb base c, rb: bulb base root (final point)
{ //"findBulbBaseOri"
  //note: cb, rb, xc are not this->cb, rb, xc when called from Orbit
  if (period2==1)
  {
    first_cb.zero(0.25, 0);
    dbg_first_rb.zero(0.5, 0);
    res_cb.zero(0.25, 0);
    res_rb.zero(0.5, 0);
    res_xc.zero(0, 0);
    res_baseZC.zero(0, 0);
    res_baseCC.zero(0, 0);
    res_baseFz.zero(0, 0);
    res_card=true;
    res_foundMult=1;//other atoms say 0/1   2;
    res_valid=true;
    return;// true;
  };
  int period=period2;
  res_xc.assign(c);
  res_card=false;
  bool did_reduce_period=false;

  /*MandelMath::complex f(currentWorker, &bulb.bulbe.f_re, &bulb.bulbe.f_im, true);
  MandelMath::complex f_z(currentWorker, &bulb.bulbe.f_z_re, &bulb.bulbe.f_z_im, true);
  MandelMath::complex f_c(currentWorker, &bulb.bulbe.f_c_re, &bulb.bulbe.f_c_im, true);
  MandelMath::complex f_zz(currentWorker, &bulb.bulbe.f_zz_re, &bulb.bulbe.f_zz_im, true);
  MandelMath::complex f_zc(currentWorker, &bulb.bulbe.f_zc_re, &bulb.bulbe.f_zc_im, true);
  MandelMath::complex f_cc_(currentWorker, &bulb.bulbe.f_cc_re, &bulb.bulbe.f_cc_im, true);
  MandelMath::complex f_zzc(currentWorker, &bulb.bulbe.f_zzc_re, &bulb.bulbe.f_zzc_im, true);
  MandelMath::complex s1(currentWorker, &bulb.s1_re, &bulb.s1_im, true);
  MandelMath::complex s2_(currentWorker, &bulb.s2_re_, &bulb.s2_im_, true);
  MandelMath::complex s3(currentWorker, &bulb.s3_re, &bulb.s3_im_, true);
  MandelMath::complex deltac(currentWorker, &bulb.deltac_re, &bulb.deltac_im, true);
  MandelMath::complex deltar(currentWorker, &bulb.deltar_re, &bulb.deltar_im, true);
  MandelMath::complex target_f_z(currentWorker, &bulb.target_f_z_re, &bulb.target_f_z_im, true);
  currentWorker->zero(target_f_z.re_s, 1);
  currentWorker->zero(target_f_z.im_s, 0);*/
  target_f_z.zero(1, 0);
  res_foundMult=1;
  /*MandelMath::complex g(currentWorker, &bulb.bulbe.g_re, &bulb.g_im, true);
  MandelMath::complex g_c(currentWorker, &bulb.g_c_re, &bulb.g_c_im, true);
  MandelMath::complex g_c2(currentWorker, &bulb.g_c2_re, &bulb.g_c2_im, true);
  MandelMath::complex g_cc(currentWorker, &bulb.g_cc_re, &bulb.g_cc_im, true);*/
  //1) find bulb center
  for (int cyc=0; cyc<10; cyc++)
  {
    if (!loope->evalg(period, res_xc))
    {
      res_valid=false;
      return;// false;
    }
    double g_mag=loope->f.getMag_double();
    //we know multiplicity==2: first twice the newton's step, then solve for f'=0 using x:=x-f'/f''
    double g_c_mag=loope->f_c.getMag_tmp()->toDouble();
    double g_cc_mag;
    if (g_mag>c.re.eps2()*1000)
    { //xc=xc-2*g/g_c
      g_cc_mag=1;
      loope->f_c.recip_prepared();
      loope->f.mul(loope->f_c);
      loope->f.lshift(1);
      res_xc.sub(loope->f);
    }
    else
    { //xc=xc-g_c/g_cc
      g_cc_mag=loope->f_cc.getMag_tmp()->toDouble();
      loope->f_cc.recip_prepared();
      loope->f_c.mul(loope->f_cc);
      res_xc.sub(loope->f_c);
    }
    if (g_c_mag<g_cc_mag*c.re.eps2()*2*(res_xc.re.radixfloor()+
                                        res_xc.im.radixfloor()))
      break;
    if (cyc==9)
      nop();
  }
  //f=0 f_c~0 f_cc=-0.00860858-i0.03615690
  //deus ex machina   bulb 1/4: cent~0.2822713907669139+i0.5300606175785253 base c~0.249+i0.500 r~-0.01+i0.499
  /*currentWorker->assign(rb->re_s, xc->re_s);
  currentWorker->assign(rb->im_s, xc->im_s);
  currentWorker->assign(cb->re_s, xc->re_s);
  currentWorker->assign(cb->im_s, xc->im_s);
  bulb.bulbe.eval2(period, cb, rb);*/ //currentWorker->sub(f.re_s, rb->re_s); currentWorker->sub(f.im_s, rb->im_s); currentWorker->add_double(f_z.re_s, -1);
  /*currentWorker->assign(s1.re_s, &bulb.bulbe.f_zc_re);
  currentWorker->assign(s1.im_s, &bulb.bulbe.f_zc_im);
  currentWorker->add(s1.re_s, &bulb.bulbe.f_zz_re);
  currentWorker->add(s1.im_s, &bulb.bulbe.f_zz_im);
  s1.recip();
  cb->add(&s1);
  currentWorker->lshift(s1.re_s, 1);
  currentWorker->lshift(s1.im_s, 1);
  rb->add(&s1);*/

  /*bulb.bulbe.eval2(period, cb, rb);
  currentWorker->assign(s1.re_s, cb->re_s);
  currentWorker->assign(s1.im_s, cb->im_s);
  currentWorker->sub(s1.re_s, xc->re_s); //s1=cb-xc
  currentWorker->sub(s1.im_s, xc->im_s);
  currentWorker->assign(s3.re_s, &bulb.bulbe.f_zz_re);
  currentWorker->assign(s3.im_s, &bulb.bulbe.f_zz_im);
  s3.mul(&f_z);
  currentWorker->rsub(s3.re_s, &bulb.bulbe.f_zc_re); //s3=f_zc-f_z*f_zz
  currentWorker->rsub(s3.im_s, &bulb.bulbe.f_zc_im);
  currentWorker->assign(s2_.re_s, f_z.re_s); //s2=f_z+1
  currentWorker->assign(s2_.im_s, f_z.im_s);
  s2_.recip();
  s2_.mul(&s1);
  s2_.mul(&s3);
  bulb.dbg_guessmult=currentWorker->toDouble(s2_.re_s); //1+1/x*/


  /*
  find r, c where f=0 fz=1
  we already have f=0 so just keep that: fz*(r-xc)+fc*(c-xc)=0
      //but fz(xc,xc)==0 fc(xc,xc)==1 so we need 2nd derivatives
      //fz*(r-xc)+fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+fc*(c-xc)+fcc/2*(c-xc)^2=0
      //fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+(c-xc)+fcc/2*(c-xc)^2=0
      actually fz(xc,xc)=-1, we need fz=0
      fzz*(r-xc)+fzc*(c-xc)=1 //move from fz=-1 to fz=0
      fz*(r-xc)+fc*(c-xc)=0   //keep f=0
      solve
      fzz*(r-xc)+fzc*(c-xc)=1
      (r-xc)=(c-xc)
      (fzz+fzc)*(c-xc)=1 correct but uhh delta r=delta c ?
      2nd derivatives
      fzz*(r-xc)+fzc*(c-xc)=1 //move from fz=-1 to fz=0
      fz*(r-xc)+fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+fc*(c-xc)+fcc/2*(c-xc)^2=0   //keep f=0
      solve, S=c-xc
      fzz*(r-xc)=(1-fzc*S)    8*0.5=1+4*0.25  4=2
      -fzz*(r-xc)/fzz+fzz/2*fzz*(r-xc)^2/fzz+fzc*fzz*(r-xc)/fzz*S+S+fcc/2*S^2=0
      -(1-fzc*S)/fzz+1/2*(1-fzc*S)^2/fzz+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0   |*2*fzz
      -2*(1-fzc*S)+(1-fzc*S)^2+2*fzc*(1-fzc*S)*S+2*fzz*S+fzz*fcc*S^2=0
      2*(fzz+fzc)*S+(fzz*fcc-fzc*fzc)*S^2-1=0
      fzz*fcc-fzc*fzc=0 at r=c=xc
      2*(fzz+fzc)*S=1

  being at xr, xc assume fzz=fz=0 at target r0, c0: replace fz*(r-xr)+fzz/2*(r-xr)^2 with fzzz*(r-r0)^3 =fzzz*(r-r0)*(r-r0)^2=fzz*(r-r0)^2
  fz~3*fzzz*(xr-r0)^2  fzz~6*fzzz*(xr-r0)

  xc=-0.76 -> xr=1/10 (5 - sqrt(101))
  fz=0.02004975155164389191229403490 fc=-0.00997512422417805404385298255 fzz=0.02014925465493167573688210469  fzc=-2.0199502484483561080877059651038304 fcc=2
  at xc=-0.76 r0=-0.5: f=0.0001  r0-xr=0.004987562112089027021926491275957618694502347002637729057282829...

  S=c0-xc
  fzz/2*(r0-xr)+fzc*S=-fz   //move from fz to fz=0       4*0.5-4*0.25=1  2-1=1
  fzz*(r0-xr)=2*(-fz-fzc*S)
  fz*(r0-xr)+fzz/2*(r0-xr)^2+fzc*(r0-xr)*S+fc*S+fcc/2*S^2=0
  -fzz*(r0-xr)^2+fzc*(r0-xr)*S+fc*S+fcc/2*S^2=0
  -4*(-fz-fzc*S)*(-fz-fzc*S)+fzc*2*(-fz-fzc*S)*S+fc*S*fzz+fzz*fcc/2*S^2=0
  -4*fz^2-4*2*fz*fzc*S-4*fzc*fzc*S*S-2*fzc*fz*S-2*fzc*fzc*S*S+fc*S*fzz+fzz*fcc/2*S^2=0
  -4*fz^2+(fc*fzz-10*fz*fzc)*S+(fzz*fcc/2-6*fzc*fzc)*S^2=0
  4+(8-10*4)*S+(8*2-6*16)*S^2 = 4-32S-80S^2 = 1-8S-10S^2



                        1-|fz+1|^2
  // c0-c=  -----------------------
  //            | fzc + fzz fc/fz |
  if?
  fzz/2*(r-r0)+fzc*(c-c0)=-fz


  change in fz = 1 = fzz*(r-xc)+fzc*(c-xc)
  solve
    fzz*(r-xc)+fzc*(c-xc)=1
    fzz/2*(r-xc)^2+fzc*(r-xc)*(c-xc)+(c-xc)+fcc/2*(c-xc)^2=0
    |  S=c-xc
    fzz*(r-xc)=(1-fzc*S)
    1/fzz/2*fzz^2*(r-xc)^2+fzc*fzz*(r-xc)/fzz*S+S+fcc/2*S^2=0
    1/fzz/2*(1-fzc*S)^2+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0    |*2*fzz
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0
    bulb 1/2: fz=-1  fc=1  fzz=8  fzc=-4  fcc=2  fzcc=0  fzzc=4  fzzz=-24
    correct guess 1/(fzz+fzc)=1/(8-4)=1/4   -1+1/4=-0.75 correct
            8*0.5=(1+4*0.25)   4=1+1

    base of bulb 1/2: c=-0.75 r=-0.5 fz=0 fc=0 fzz=0 fzc=-2 fcc=2  fzcc=0  fzzc=4  fzzz=-12   evaluate d^2/dz^2 (z^2+c)^2+c-z at c=-0.75 and z=-0.5



    (r-xc)=(1-fzc*S)/fzz
    fzz/2*(1-fzc*S)/fzz*(1-fzc*S)/fzz+fzc*(1-fzc*S)/fzz*S+S+fcc/2*S^2=0  |*2*fzz
    (1-fzc*S)*(1-fzc*S)+2*fzc*(1-fzc*S)*S+2*fzz*S+fzz*fcc*S^2=0
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0
    1+2*fzz*S+(fzz*fcc-fzc*fzc)*S^2=0

    fzz^2*(1-fzc*S)/fzz*(1-fzc*S)/fzz+2*fzz*fzc*(1-fzc*S)/fzz*S+S*fzz*2+fcc/2*S^2=0
    1/2*(1-fzc*S)*(1-fzc*S)/fzz+2*fzc*(1-fzc*S)*S+S+fcc/2*S^2=0
    (fcc*fzz-fzc*fzc)*S^2+2*(fzc+fzz)*S+1=0
      proof? that fcc*fzz-fzc^2=0
      w.l.o.g f=a*x^2+b*x*y+c*y^2
      fxx=2a fyy=2c fxy=b
      not really



    fzz/2*R+fzc*C+fz=0    f=fzzz/6*R^3  fz=3*fzzz/6*R^2  fzz=6*fzzz/6*R  fzzz=fzzz  2*fz/fzz=R  2*-1/8=-1/4!=R
    fc=1+fcc*C+fzc*R  0=1+2*0.25-(4+2)/2*0.5=1+0.5-1.5
    lower case at xr, xc; upper case at r0, c0  C=c-c0 R=r-r0
    f=0+FC*C+FZC*R*C+FCC*C^2/2+FZZZ*R^3/6+FZZC*R^2*C/2  (FZ=FZZ=0  (FZZZ),FZZC,FZCC,FCC don't change much)
    fz=FZC*C+FZZZ*R^2/2+FZZC*R*C    -1=-2*-0.25+(-12..-24)*0.5^2/2+4*-0.5*-0.25=1+(-1.5..-3)=-0.5..-2~-1
    fzz=FZZZ*R+FZZC*C               8=(-12..-24)*-0.5+4*-0.25=6..12-1=5..11~8
    fc=FC+FZC*R+FCC*C+FZZC*R^2/2    1=0+-2*-0.5+2*-0.25+4*0.5^2/2=1-0.5+1/2=1
    fzc=FZC+FZZC*R                  -4=-2+4*-0.5=-2-2=-4
    fcc=FCC            -> FC FZC FCC FZZZ FZZC R C
    fzzc=FZZC
    now find R, C from fz..fzzc and f=0
    0=FC*C+(fzc-fzzc*R)*R*C+fcc*C^2/2+(fzz-fzzc*C)*R^2/6+fzzc*R^2*C/2
    fz=(fzc-fzzc*R)*C+(fzz-fzzc*C)*R/2+fzzc*R*C
    //(fzz-fzzc*C)=FZZZ*R
    fc=FC+(fzc-fzzc*R)*R+fcc*C+fzzc*R^2/2  |*C-(1)
    //(fzc-fzzc*R)=FZC
    -
    //fc*C= FC*C+(fzc-fzzc*R)*R*C+fcc*C*C+fzzc*R^2*C/2
    //0   =-FC*C-(fzc-fzzc*R)*R*C-fcc*C*C/2-(fzz-fzzc*C)*R^2/6-fzzc*R^2*C/2
    fc*C=fcc*C*C/2-fzz*R^2/6+fzzc*R^2*C/6
    fz=fzc*C+fzz*R/2-fzzc*C*R/2
    -
    2*fzc*C-2*fz=(-fzz+fzzc*C)*R   (2*fzc*C-2*fz)/(-fzz+fzzc*C)=R
    fc*C=fcc*C*C/2+(-fzz+fzzc*C)*R*R/6
    -
    fc*C=fcc*C*C/2+(2*fzc*C-2*fz)*(2*fzc*C-2*fz)/(-fzz+fzzc*C)/6
    -
    0=4*fz*fz+(6*fc*fzz-8*fz*fzc)*C+(4*fzc*fzc-3*fcc*fzz-6*fc*fzzc)*C*C+3*fzzc*fcc*C*C*C
    at 1/2 bulb center: ->C=-0.212873 R=-0.418346
    first newton step: 4*fz*fz/(6*fc*fzz-8*fz*fzc), at bulb center 2/(3*fzz+4*fzc)

    FC=d FZC=e FCC=f FZZZ=g FZZC=h
    f=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2
    fz=e*C+g*R^2/2+h*R*C
    fzz=g*R+h*C
    fc=d+e*R+f*C+h*R^2/2
    fzc=e+h*R
    fcc=f
    fzzc=h
    -> C≈-0.212873 R≈-0.418346 fc≈0.102389 fzc≈-2.32662 fcc=2 fzzz≈-21.1583 fzzc=4
    correct C=-0.25 R=-0.5     fc=0        fzc=-2       fcc=2 fzzz=-12      fzzc=4
  */

  //so let's find r0, c0 aka rb, cb
  res_rb.assign(res_xc);
  res_cb.assign(res_xc);
  firstStep=0;
  //double test_x0_mag=0;
  for (int cycle=0; cycle<6; cycle++)
  {
#if 0 //attempt at 2-nd order approximation... quite a fail, and not needed any more since we reduce period later
    bulb.bulbe.eval2zzc(period, cb, rb);
    currentWorker->sub(f.re_s, rb->re_s);
    currentWorker->sub(f.im_s, rb->im_s);
    currentWorker->add_double(f_z.re_s, -1);
    //0=4*fz*fz+(6*fc*fzz-8*fz*fzc)*C+(4*fzc*fzc-3*fcc*fzz-6*fc*fzzc)*C*C+3*fzzc*fcc*C*C*C
    currentWorker->assign(s2_.re_s, f_zzc.re_s);
    currentWorker->assign(s2_.im_s, f_zzc.im_s);
    s2_.mul(&f_cc_);
    double A3_re=3*currentWorker->toDouble(s2_.re_s);
    double A3_im=3*currentWorker->toDouble(s2_.im_s);

    currentWorker->assign(s2_.re_s, f_zzc.re_s);
    currentWorker->assign(s2_.im_s, f_zzc.im_s);
    s2_.mul(&f_c);
    currentWorker->lshift(s2_.re_s, 1);
    currentWorker->lshift(s2_.im_s, 1);
    currentWorker->assign(s1.re_s, f_zz.re_s);
    currentWorker->assign(s1.im_s, f_zz.im_s);
    s1.mul(&f_cc_);
    s2_.add(&s1);
    currentWorker->mul_double(s2_.re_s, -3);
    currentWorker->mul_double(s2_.im_s, -3);
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    currentWorker->lshift(s1.re_s, 1);
    currentWorker->lshift(s1.im_s, 1);
    s1.sqr();
    s2_.add(&s1);
    double B3_re=currentWorker->toDouble(s2_.re_s);
    double B3_im=currentWorker->toDouble(s2_.im_s);

    currentWorker->assign(s2_.re_s, f_zc.re_s);
    currentWorker->assign(s2_.im_s, f_zc.im_s);
    s2_.mul(&f_z);
    currentWorker->lshift(s2_.re_s, 3);
    currentWorker->lshift(s2_.im_s, 3);
    currentWorker->assign(s1.re_s, f_zz.re_s);
    currentWorker->assign(s1.im_s, f_zz.im_s);
    s1.mul(&f_c);
    currentWorker->mul_double(s1.re_s, 6);
    currentWorker->mul_double(s1.im_s, 6);
    currentWorker->sub(s1.re_s, s2_.re_s);
    currentWorker->sub(s1.im_s, s2_.im_s);
    double C3_re=currentWorker->toDouble(s1.re_s);
    double C3_im=currentWorker->toDouble(s1.im_s);

    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    currentWorker->lshift(s2_.re_s, 1);
    currentWorker->lshift(s2_.im_s, 1);
    s2_.sqr();
    double D3_re=currentWorker->toDouble(s2_.re_s);
    double D3_im=currentWorker->toDouble(s2_.im_s);

    //initial guess between (1-(fz+1)^2)/(fzc-fzz*fc/fz) and 1/2x that -> mul by (1+|fz|)/2
    //-4*fz-4*fz*fz-fz^3    (2-1)*(1-(-1+1)^2)/(-4+8*1/-1)=1/(-4+8*1/-1)
    //(1+fzmag)/2*fz^2*(2+fz))/(fzc*fz-fzz*fc)  (1+1)/2*1*(2-1)/(-4*-1-8*1)=1/-4
    //deltac:=fz^2*(4+4*fz+fz^2)/(fz*fzc-fc*fzz)   -(4-4+1)/(-1*-4-1*8)=-1/(-4)=0.25
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    s1.mul(&f_z);
    currentWorker->assign(s2_.re_s, f_zz.re_s);
    currentWorker->assign(s2_.im_s, f_zz.im_s);
    s2_.mul(&f_c);
    currentWorker->sub(s1.re_s, s2_.re_s);
    currentWorker->sub(s1.im_s, s2_.im_s);
    s1.recip();
    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    s2_.sqr();
    s1.mul(&s2_); //fz^2/(fz*fzc-fc*fzz)
    /*currentWorker->assign(s3.re_s, f_z.re_s);
    currentWorker->assign(s3.im_s, f_z.im_s);
    currentWorker->add_double(s3.re_s, 1);
    currentWorker->lshift(s3.re_s, 2);
    currentWorker->lshift(s3.im_s, 2);
    s2_.add(&s3);*/
    currentWorker->assign(s2_.re_s, f_z.re_s);
    currentWorker->assign(s2_.im_s, f_z.im_s);
    currentWorker->add_double(s2_.re_s, 2);
    s1.mul(&s2_); // *(2+fz)
    currentWorker->assign(s2_.re_s, f_z.getMagTmp());
    currentWorker->add_double(s2_.re_s, 1);
    currentWorker->lshift(s2_.re_s, -1);
    currentWorker->zero(s2_.im_s);
    s1.mul(&s2_); // *(1+fz_mag)/2   not sure what function is best but we need fzmag=1 -> *1 fzmag=0 -> *0.5
    double Z_re=currentWorker->toDouble(s1.re_s); //-0.158439 -> -0.03827124; 0.00449088903772943->/4; -4.8888e-8->-1.2222e-8
    double Z_im=currentWorker->toDouble(s1.im_s);

    //now let's solve A3*Z^3+B3*Z^2+C3*Z+D3=0
    for (int cycle3=0; cycle3<9; cycle3++)
    {
      double F_re=A3_re, F_im=A3_im;
      double t;
      t=F_re*Z_re-F_im*Z_im+B3_re;
      F_im=F_im*Z_re+F_re*Z_im+B3_im;
      F_re=t; //F=A3*F+B3
      t=F_re*Z_re-F_im*Z_im+C3_re;
      F_im=F_im*Z_re+F_re*Z_im+C3_im;
      F_re=t; //F=(A3*Z+B3)*Z+C3
      t=F_re*Z_re-F_im*Z_im+D3_re;
      F_im=F_im*Z_re+F_re*Z_im+D3_im;
      F_re=t; //F=((A3*Z+B3)*Z+C3)*Z+D3
      double Fmag=(F_re*F_re+F_im*F_im);
      if (Fmag==0)
        break;
      t=1.0/Fmag;
      double f1_re=F_re*t;
      double f1_im=-F_im*t; //f1=1/F

      double laguG_re=3*A3_re;
      double laguG_im=3*A3_im; //laguG=3*A3
      t=laguG_re*Z_re-laguG_im*Z_im+2*B3_re;
      laguG_im=laguG_im*Z_re+laguG_re*Z_im+2*B3_im;
      laguG_re=t; //laguG=3*A3*Z+2*B3
      t=laguG_re*Z_re-laguG_im*Z_im+C3_re;
      laguG_im=laguG_im*Z_re+laguG_re*Z_im+C3_im;
      laguG_re=t; //laguG=(3*A3*Z+2*B3)*Z+C3=f' = 3*A3*Z^2+2*B3*Z+C3
      t=laguG_re*f1_re-laguG_im*f1_im;
      laguG_im=laguG_im*f1_re+laguG_re*f1_im;
      laguG_re=t; //laguG=f'/f

      double fzzf_re=6*A3_re;
      double fzzf_im=6*A3_im;
      t=fzzf_re*Z_re-fzzf_im*Z_im+2*B3_re;
      fzzf_im=fzzf_im*Z_re+fzzf_re*Z_im+2*B3_im;
      fzzf_re=t; //fzzf=6*A3*Z+2*B3=f'' = 6*A3*Z+2*B3
      t=fzzf_re*f1_re-fzzf_im*f1_im;
      fzzf_im=fzzf_im*f1_re+fzzf_re*f1_im;
      fzzf_re=t; //fzzf=f''/f

      t=laguG_re*laguG_re-laguG_im*laguG_im;
      double laguG2_im=2*laguG_re*laguG_im;
      double laguG2_re=t; //laguG2=laguG^2
      // laguH=fzf^2-fzzf
      // m=Round( Re(G^2*H^T)/mag(H) )
      // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
      double laguH_re=laguG2_re-fzzf_re;
      double laguH_im=laguG2_im-fzzf_im;

      // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      double laguX_re=(laguH_re-laguG2_re/3.0)*(2/3.0);
      double laguX_im=(laguH_im-laguG2_im/3.0)*(2/3.0); //laguX=(1/m-1/n)*(H-G^2/n)  m=1 n=3
      MandelMath::complex_double_sqrt(&laguX_re, &laguX_im, laguX_re, laguX_im);
      //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
      t=laguX_re*laguG_re+laguX_im*laguG_im; //laguX.mulreT(&laguG)
      if (t<0)
      {
        laguX_re=-laguX_re;
        laguX_im=-laguX_im;
      };
      laguX_re+=laguG_re/3.0;
      laguX_im+=laguG_im/3.0; //(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      t=1/(laguX_re*laguX_re+laguX_im*laguX_im);
      laguX_re=laguX_re*t;
      laguX_im=-laguX_im*t; //1/(...)

      double oldzre=Z_re; Z_re-=laguX_re;
      double oldzim=Z_im; Z_im-=laguX_im;
      if ((Fmag<1e-37) || (Z_re==oldzre && Z_im==oldzim))
        break;
    }
    //now Z is the step in C (in cb)
    //(2*fzc*C-2*fz)/(fzzc*C-fzz)=R
    currentWorker->zero(deltac.re_s, Z_re);
    currentWorker->zero(deltac.im_s, Z_im);
    currentWorker->assign(deltar.re_s, deltac.re_s);
    currentWorker->assign(deltar.im_s, deltac.im_s);
    deltar.mul(&f_zc);
    currentWorker->sub(deltar.re_s, f_z.re_s);
    currentWorker->sub(deltar.im_s, f_z.im_s);
    currentWorker->lshift(deltar.re_s, 1);
    currentWorker->lshift(deltar.im_s, 1);
    currentWorker->assign(s1.re_s, deltac.re_s);
    currentWorker->assign(s1.im_s, deltac.im_s);
    s1.mul(&f_zzc);
    currentWorker->sub(s1.re_s, f_zz.re_s);
    currentWorker->sub(s1.im_s, f_zz.im_s);
    s1.recip();
    deltar.mul(&s1); //R
#else //4-card c~-0.154723+I*1.031046 r~-0.153526+I*1.029663
    loope->eval2(period, res_cb, res_rb, true); //TODO: don't need f_cc here
    if (cycle==0)
    { //test if period needs a cleanup
      //https://www.researchgate.net/publication/337838756_The_size_of_Mandelbrot_bulbs
      //The size of Mandelbrot bulbs - A.C. Fowler, M.J. McGuinness
      //radius of (p/q) at main card ~ 1/q^2 sin(pi p/q)
      //  =>  radius of bulb <= ~ 1/q^2 = 1/period^2
      //  largest I found is 1.289/period^2 (lim=0.602 per^4) for p=1024, bulb=(3/8)^3*1/2
      //                     0.474 at xc=-0.7496491297699771+i0.02178217774381104 p=1442
      //
      t1.assign(loope->f_zz);
      t1.add(loope->f_zc); //usually f_zz+f_zc is large (small bulb), but close to 0 (very large bulb) at multiple iterations
      double radius_2=t1.getMag_double(); //radius^-2
      double limit=(period*0.474*period)*period*period; // 1/radius^2=1/(1.4/period^2)^2=0.5 period^4
      //
      if (radius_2<limit && 4<radius_2)
        nop();
      if (radius_2<4) //smallest for main cardioid, everything else should be larger
                         //
      {
        nop();
        t1.assign(res_rb);
        double near_limit=loope->f.getMag_double()+25*25*c.re.eps2();
        for (int i=1; i<period; i++) //don't do last period, we know it works
        {
          t1.sqr();
          t1.add(res_cb);
          double dist=t1.dist2_double(res_rb);
          if (dist<near_limit) //543? at period=1342154 i=61007
          {
            nop();
            period=i;
            loope->eval2(period, res_cb, res_rb, true);
            break;
          }
        }
        if (period==period2)
          nop(); //shortening is pretty much guaranteed here
      };
    };
#endif
    //for cardioid, we need to use f_zz as well in the equation for f because it does not go to 0
    //0=f=0+C*fc+R*fz+R*R*fzz/2    C=-R*(fz+R*fzz/2)/fc
    //target_fz-1=fz+C*fzc+R*fzz   target_fz-1=fz-R*(fz+R*fzz/2)/fc*fzc+R*fzz    target_fz-1-fz=-R*R*fzz/fc*fzc/2+R*(fzz-fz/fc*fzc)
    //0=R*R*fzz/2/fc*fzc+R*(fz/fc*fzc-fzz)+target_fz-fz-1  : a=fzz/2/fc*fzc b=(fz/fc*fzc-fzz) c=target_fz-fz-1
    //x1,2= -(2*c)/(b+-sqrt(b^2-2*a*2*c))         good except both b,c small e.g. 0
    t1.assign(loope->f_c);
    t1.recip();
    t1.mul(loope->f_zc);
    t2.assign(t1); //f_zc/f_c
    t1.mul(loope->f_z);
    t1.sub(loope->f_zz); //s1=b
    t2.mul(loope->f_zz); //s2=2*a
    t3.assign(loope->f_z);
    t3.sub(target_f_z);
    t3.re.add_double(1); //(fz-target_fz+1)
    t3.lshift(1); //s3=-2*c
    t2.mul(t3);
    deltar.assign(t1);
    deltar.sqr();
    deltar.add(t2); //b^2-4ac
    deltar.sqrt();
    if (deltar.mulreT_tmp(t1)->toDouble()<0)
    {
      deltar.chs();
    };
    deltar.add(t1);
    deltar.recip();
    deltar.mul(t3);
    //C=-R*(fz+R*fzz/2)/fc
    deltac.assign(loope->f_zz);
    deltac.mul(deltar);
    deltac.lshift(-1);
    deltac.add(loope->f_z);
    deltac.mul(deltar);
    t1.assign(loope->f_c);
    t1.recip();
    deltac.mul(t1);
    deltar.chs();
    t1.assign(deltac);
    t2.assign(deltar);

    //check the easy way:
    //        0=f=0+C*fc+R*fz -> R=-C*fc/fz
    //target_fz-1=fz+C*fzc+R*fzz   fz/(fc/fz*fzz-fzc)=C=-deltaC    (fz-target_fz+1)/(fzc-fc/fz*fzz)=deltaC
    deltac.assign(loope->f_z);
    deltac.recip();
    deltac.mul(loope->f_c);
    deltar.assign(deltac); // fc/fz
    deltac.mul(loope->f_zz);
    deltac.rsub(loope->f_zc);
    { //check fzc-fc/fz*fzz before .recip()
      double mag=deltac.getMag_tmp()->toDouble();
      if (mag<1e-8)
      {
        nop(); //if period is a multiple of real period, in first cycle we get f_zz~0, f_zc~0 -> deltac~infinity
        res_valid=false;
        return;// false; //also happens in later cycles, not sure why the first cycle doesn't catch it
      }
      else if (mag<0.1)
        nop();
    }
    deltac.recip_prepared(); // 1/(fzc-fc/fz*fzz)
    t3.assign(loope->f_z);
    t3.sub(target_f_z);
    t3.re.add_double(1); //(fz-target_fz+1)
    deltac.mul(t3);  //deltac
    deltar.mul(deltac);
    deltar.chs(); // deltar

#if 0
    //at card, we are not at result yet so f_zz/f_z does not blow clearly enough
    if (!did_reduce_period && !*is_card) //after reduction, we can arrive at cardioid but must not report it any more
    { //preferably we need to decide at cycle 1 because otherwise it never converges
      currentWorker->sub(&bulb.test_xn_re, &bulb.test_x0_re);
      currentWorker->sub(&bulb.test_xn_im, &bulb.test_x0_im);
      double test_xn=MandelMath::sqr_double(currentWorker->toDouble(&bulb.test_xn_re))+MandelMath::sqr_double(currentWorker->toDouble(&bulb.test_xn_im));
      if (test_xn>test_x0_mag*100)
        *is_card=true;
      else if (test_xn>test_x0_mag*4.1)
        nop(); //?
      else if (test_xn*4.1>test_x0_mag)
        nop(); //xn~x0
      else
        nop(); //xn<x0 should not be
      /*double test_x1=currentWorker->toDouble(f_zz.getMagTmp());
      double test_x2=currentWorker->toDouble(f_zz.getMagTmp())/currentWorker->toDouble(f_z.getMagTmp())*currentWorker->toDouble(f_c.getMagTmp());
        //does not work, maybe f_zz^2
      double test_zc=currentWorker->toDouble(f_zc.getMagTmp());
      if (test_x1*1000<test_zc)
        *is_card=true;
      else if (test_x1<test_zc) // *0.25 is bulb
        nop();
      else
        nop();
      (void)test_x1;
      (void)test_x2;*/
    }
#endif
    if (res_card)
    {
      //no idea why but the quadratic approximation gives deltac=easy_deltac/2, deltar=easy_deltar
      //easy/2 actually seems better than quadratic approx
      /*currentWorker->assign(deltac.re_s, s1.re_s);
      currentWorker->assign(deltac.im_s, s1.im_s);
      currentWorker->assign(deltar.re_s, s2_.re_s);
      currentWorker->assign(deltar.im_s, s2_.im_s);*/
      deltac.lshift(-1);
    }
    else if (cycle>0)
    {
      t3.assign(res_cb);
      t3.sub(deltac);
      if (t3.dist2_double(first_cb)>firstStep/16)
      { //either wrong reduction of period, or reduced all the way to previous bulb - so we're done
        dbgPoint(); //should be caught at the end of last cycle
        res_rb.assign(prev_rb);
        //TODO: baseZC, baseCC are wrong as well
        break;
      };
    }
    else
    {
      firstStep=deltac.getMag_double();
    }



    //s2:=f_zc^2-f_zz*f_cc
    /* currentWorker->assign(s2_.re_s, f_zz.re_s);
    currentWorker->assign(s2_.im_s, f_zz.im_s);
    s2_.mul(&f_cc_);
    currentWorker->assign(s1.re_s, f_zc.re_s);
    currentWorker->assign(s1.im_s, f_zc.im_s);
    s1.sqr();
    currentWorker->sub(s2_.re_s, s1.re_s);
    currentWorker->sub(s2_.im_s, s1.im_s);
    //g always=0+0i at the bulb center (or it should) and 1+0i at cardioid
    double cardioid_discrim=currentWorker->toDouble(s2_.getMagTmp());
    double precis=currentWorker->toDouble(f_zc.getMagTmp());
    if (cardioid_discrim<precis*1e-10)
      *is_card=false;
    else if (cardioid_discrim<1-precis*1e-10)
      *is_card=false; //?
    else if (cardioid_discrim<1+precis*1e-10)
      *is_card=true;
    else
      *is_card=true; //? */


    res_cb.sub(deltac);
    res_rb.sub(deltar);
    if (cycle==0)
    {
      //try to snap to the central root, makes everything work
      loope->eval_near0(period, res_cb);
      int step=MandelMath::gcd(period, loope->near0iter_1);
      //loope->first_multi.zero(0, 0);
      loope->eval_zz(0, res_cb, res_cb, false, true);
      loope->sumA.assign(res_cb);
      for (int i=step; i<period; i+=step)
      {
        loope->eval_zz(step, res_cb, loope->first_multi, false, false);
        loope->sumA.add(loope->f);
      }
      loope->sumA.mul_double(step/double(period));
      res_rb.assign(loope->sumA);


      first_cb.assign(res_cb);
      dbg_first_rb.assign(res_rb);
    };


/*
    //verification
    would have to compute all uppercase, which I don't really need for anything else
    f=0+FC*C+FZC*R*C+FCC*C^2/2+FZZZ*R^3/6+FZZC*R^2*C/2  (FZ=FZZ=0  (FZZZ),FZZC,FZCC,FCC don't change much)
    fz=FZC*C+FZZZ*R^2/2+FZZC*R*C    -1=-2*-0.25+(-12..-24)*0.5^2/2+4*-0.5*-0.25=1+(-1.5..-3)=-0.5..-2~-1
    fzz=FZZZ*R+FZZC*C               8=(-12..-24)*-0.5+4*-0.25=6..12-1=5..11~8
    fc=FC+FZC*R+FCC*C+FZZC*R^2/2    1=0+-2*-0.5+2*-0.25+4*0.5^2/2=1-0.5+1/2=1
    fzc=FZC+FZZC*R                  -4=-2+4*-0.5=-2-2=-4

           f                                        fz                    fzz        fc                   fzc       fcc  fzzc
    solve [0=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2, -1=e*C+g*R^2/2+h*R*C, 8=g*R+h*C, 1=d+e*R+f*C+h*R^2/2, -4=e+h*R, 2=f, 4=h]
solve [0=0+d*C+e*R*C+f*C^2/2+g*R^3/6+h*R^2*C/2, -1=e*C+g*R^2/2+h*R*C, -124.796127+I*6.194576=g*R+h*C, 1=d+e*R+f*C+h*R^2/2, 48.31925+I*48.506833=e+h*R, 2.0051-I*37.4627558=f, -1020.0523+I*3920.8781217=h]
*/

    t1.assign(res_rb);

    //now we guessed step in cb and rb, need to tune rb so that f(cb, rb)=0 again
    double prev_f_mag=1e10;
    for (int lcycle=0; lcycle<10; lcycle++)
    {
      //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
      //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
      //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
      if (!loope->eval_zz(period, res_cb, res_rb, true))
      { //happens way too often... deltar guess suks
        res_foundMult=1;
        res_valid=false;
        return;// false;
      };
      //bulbe.f.sub(rb);
      //currentWorker->add_double(bulbe.f_z.re, -1);
#if 0 //worse than without; leave fixing of f_z to the 2D newton above
      if (currentWorker->toDouble(f.getMagTmp())<3*currentWorker->eps2())
      {
        //half-step to improve f'==target: rb:=rb-(f'-target)/(f'')
        currentWorker->assign(s2_.re_s, f_z.re_s);
        currentWorker->assign(s2_.im_s, f_z.im_s);
        currentWorker->rsub(s2_.re_s, target_f_z.re_s);
        currentWorker->rsub(s2_.im_s, target_f_z.im_s);
        currentWorker->add_double(s2_.re_s, -1);
        currentWorker->assign(s3.re_s, f_zz.re_s);
        currentWorker->assign(s3.im_s, f_zz.im_s);
        s3.recip();
        s2_.mul(&s3);
        rb->add(&s2_);
        break;
      }
#endif
      if (loope->f.is0())
        break;
      double f_mag=loope->f.getMag_double();
      if (f_mag>=prev_f_mag)
      {
        lagu.step.lshift(-1);
        res_rb.add(lagu.step);
      }
      else
      {
        prev_f_mag=f_mag;
        if (!lagu.eval(period, loope->f, loope->f_z, loope->f_zz))
        {
          res_foundMult=1;
          res_valid=false;
          return;// false;
        };
        res_rb.sub(lagu.step);
        if (f_mag<3*c.re.eps2())
          break;
      }
    }

    if (cycle==0)
    { //no idea why but first guess, after fixing rb, has f_z either 0+0i at bulb, or 0+-i at cardioid
      t2.assign(loope->f_z);
      double dist_to_0=t2.getMag_double();
      t2.im.add_double(1);
      double dist_to_ni=t2.getMag_double();
      t2.im.assign(loope->f_z.im);
      t2.im.add_double(-1);
      double dist_to_pi=t2.getMag_double();
      if (dist_to_0<0.01)
        nop(); //good, bulb
      else if (dist_to_ni<0.186 || dist_to_pi<0.186)
        //0.185 per=27720 at r=-1.2623031085342689 im=0.38359354974812726
        //0.121 [actually 0.00061 in float128] per=2310 at r=-1.2623031085359033 i=0.3835935497471468
        //0.0217 per=100 at 0.249571+i0.562510
        //0.0332 per=280 at -0.5921925+i0.6198388
        res_card=true;
      else if (dist_to_0<0.25)
        nop(); //bulb but not so certain
      else if (dist_to_ni<0.25 || dist_to_pi<0.25)
      {
        nop();
        res_card=true;
      }
      else
        nop();
    };

    prev_rb.assign(res_rb);

    double f_error, fz_error, fz_error_fromr;
    t2.assign(loope->f_z);
    t2.re.add_double(1);
    fz_error_fromr=c.re.eps2()/t2.getMag_double()/4*loope->f_zz.getMag_double();//error in (prev)r~eps2/fz_error; error in fz~err(r)^2/2*f_zz~eps2/fz_error/2*f_zz
    t2.sub(target_f_z);
    f_error=loope->f.getMag_double(); //supposed to be 0 but just in case
    fz_error=t2.getMag_double();
    if (!res_card && f_error<1e-15 && fz_error<1e-4) //close enough to bulb base to start trying
    {
#if 0
      if (!loope->eval_multi(tmp, period, &res_cb, &res_rb, &target_f_z, 0)) //toler=0, here we get close to base so the iterations should eventually be in a circle
      {
        res_foundMult=1;
        res_valid=false;
        return;// false;
      };
#else
      if (!loope->eval_near0(period-1, res_cb)) //just find near0iter, best(?) with f_c(c) not f_c(r)
                                                           //period-1 so we don't incidentally overwrite a good guess with last f... I think :-/
      {
        res_foundMult=1;
        res_valid=false;
        return;// false;
      };
      if (loope->near0iter_1<1)
        loope->multi=1;
      else
      {
        loope->multi=period/MandelMath::gcd(period, loope->near0iter_1);
        loope->sumA.assign(res_rb);
      }
#endif
      if (loope->multi>1)
      {
        res_rb.assign(loope->sumA);
        //reduce period
        period /= loope->multi;

        //fix rb after moving to guessed root position, have to maintain f(cb, rb)==0
        for (int lcycle=0; lcycle<7; lcycle++)
        {
          //x^5=1 -> x:=x-(x^5-1)/(5*x^4)   NestList[#-(#^5-1)/(5*#^4) &, 0.6+0.8I, 10]
          //x=a^(1/5) |y=x/a| a*y=a^(1/5) y=a^(-4/5)  y^5=a^-4  y:=y-(y^5-a^-4)/(5y^4)  y:=(4/5)*y+a^-4/(5y^4)
          //                                          y^-5=a^4  y:=y-(y^-5-a^4)/(-5*y^-6)=y*(6-a^4*y^5)/5   x=a*y
          if (!loope->eval_zz(period, res_cb, res_rb, true))
          {
            res_foundMult=1;
            res_valid=false;
            return;// false;
          };
          //bulbe.f.sub(rb);
          //currentWorker->add_double(bulbe.f_z.re, -1);
          if (loope->f.is0())
            break;
          if (!lagu.eval(period, loope->f, loope->f_z, loope->f_zz))
          {
            res_foundMult=1;
            res_valid=false;
            return;// false;
          };
          res_rb.sub(lagu.step); //no dist from eval_near0 //TODO: fail if (first) step >smallest..average dist from eval_multi
          if (loope->f.getMag_double()<3*c.re.eps2())
            break;
        }

        //change target f_z to multi-th root of 1, using root near first_multi
        //could do a lot of tricks but let's just use newton to find the root of newf_z^multi=1
        //or rather new_target^multi=old_target
        //starting from bulb.bulbe.first_multi
        //"next_target" = "bulbe.first_multi"

        //current loope->f_z should be a much better guess for first_multi
        loope->first_multi.assign(loope->f_z);
        loope->first_multi.re.add_double(1);

        //there has to be a better way to tell if this was a refinement, or a jump to previous bulb/atom
        //easy: 2nd time always quit
        t3.assign(loope->first_multi);
        if (abs(t3.getMag1_tmp()->toDouble())>1e-2) //more like sqrt(fz_error)/multi
        {
          res_rb.assign(prev_rb);
          //TODO: baseZC, baseCC are wrong as well
          break;
        };
        t3.pow_uint(loope->multi);
        double fz1_error=t3.getDist1_tmp()->toDouble();
        bool primary=did_reduce_period;
        if ((fz1_error>2.42*fz_error) != primary) //bulb p=22032 xc=-0.1744832996+i0.6619272216 mult=2 res_c=-0.1744832974710202+i0.6619272205853156 -> fz_error=2.48e-11,1.78e-11 fz1_error=2.5840570306057777e-11=1.0418,1.4504 fz_error
          nop(); //happens rarely, still tuning the constant   2.42 at c=-0.56221793769583428+I*0.6768911909908297 p=4680->1560
        if ((fz1_error>2.42*fz_error+fz_error_fromr) != primary)
          nop(); //fz_error_fromr doesn't help
        if (primary)
        { //how close do we get on legit exits?
          if (fz1_error<5*fz_error)
            nop(); //not hit yet
          if (fz1_error<50*fz_error)
            nop(); //not hit yet
          if (fz1_error<500*fz_error)
            nop(); //not hit yet
          if (fz1_error<5000*fz_error)
            nop(); //not hit yet
          if (fz1_error<50000*fz_error)
            nop(); //not hit yet
          /* testing against a constant will not work, ever
          if (fz1_error<1e-6)
            nop(); //3.6e-9 at (1/2)*(1/52360)
          if (fz1_error<1e-4)
            nop(); //gets hit at (1/231)*(1/2)
          if (fz1_error<1e-3)
            nop(); //gets hit at (1/2)*(1/111), ...   */
          double k=fz1_error*res_foundMult*loope->multi;
          (void)k;
          /* seems k gets arbitrarily small
          if (k<0.1)
            nop(); //gets hit a lot
          if (k<0.03)
            nop(); //not any more //gets hit at (1/2)*(1/111), ...
          if (k<0.01)
            nop(); //not any more //gets hit at (1/2)*(1/111), ...
          if (k<0.003)
            nop(); //3.8e-4 at (1/2)*(1/52360)*/
          double kk=fz1_error*res_foundMult*loope->multi*res_foundMult*loope->multi;
          if (kk<4.0)
            nop(); //not hit yet
          if (kk<1.0)
            nop(); //not hit yet
          if (kk<0.1)
            nop(); //not hit yet
          if (kk<0.01)
            nop(); //not hit yet
          if (kk<0.001)
            nop(); //not hit yet
          res_rb.assign(prev_rb);
          //TODO: baseZC, baseCC are wrong as well
          break;
        }
        else
        {
          if (fz1_error>1.01*fz_error)
            nop(); //1.2595 at c=-0.5622179376996712+I*0.67689119098929496 p=8970->390
          if (fz1_error>1.003*fz_error)
            nop(); //1.003369 at c=-0.56221793788427021+I*0.6768911909821327 p=260->130
          if (fz1_error>1.001*fz_error)
            nop(); //1.00101 at preset "11016"*(1/2)
          /* testing against a constant will not work, ever
          if (fz1_error>1e-4)
            nop(); //gets as close to 1e-4 as fz_error
          if (fz1_error>1.01e-4)
            nop(); //
          if (fz1_error>1.1e-4)
            nop(); // */
          double k=fz1_error*res_foundMult*loope->multi;
          if (k>0.1)
            nop(); //not hit yet
          if (k>0.03)
            nop(); //not hit yet
          if (k>0.01)
            nop(); //not hit yet
          if (k>0.003)
            nop(); //not hit yet
          double kk=fz1_error*res_foundMult*loope->multi*res_foundMult*loope->multi;
          if (kk>1.0)
            nop(); //not hit yet
          if (kk>0.1)
            nop(); //not hit yet
          if (kk>0.01)
            nop(); //not hit yet
          if (kk>0.003)
            nop(); //0.0032 at 7/18, 0.0039 at 1/10
        }

        res_foundMult *= loope->multi;
        if (did_reduce_period)
          nop(); //did refine twice
        else
          did_reduce_period=true;
        loope->first_multi.sign(); //help Newton a little
        for (int multicyc=0; multicyc<10; multicyc++)
        {
          //a=target_f_z  y=next_target  y:=y+y*(1-a^4*y^5)/5   x=a*y
          t3.assign(loope->first_multi);
          t3.mul(target_f_z);
          //t2.assign(&t3);
          t3.pow_uint(loope->multi-1);
          t3.mul(loope->first_multi); //a^4*y^5
          t3.chs();
          t3.re.add_double(1); //(1-a^4*y^5)
          t3.mul_double(1.0/loope->multi); //(1-a^4*y^5)/5
          double dist1=t3.getMag_double();
          t3.mul(loope->first_multi); //y*(1-a^4*y^5)/5
          loope->first_multi.add(t3); //y+=y*(1-a^4*y^5)/5
          if (dist1<3*c.eps2())
            break;
        }
        target_f_z.mul(loope->first_multi);
      }
    }
    /*s1 does not hold old deltar any more
    currentWorker->sub(s1.re_s, rb->re_s);
    currentWorker->sub(s1.im_s, rb->im_s);
    currentWorker->sub(s1.re_s, deltar.re_s);
    currentWorker->sub(s1.im_s, deltar.im_s); //should be around 0*/
    if (f_error<3*c.eps2() &&
        fz_error<3*(1+loope->f_zz.getMag_double())*c.eps2() &&
        (did_reduce_period || res_card))
    {
      nop(); //ok
      break;
    };
    nop();
  }
  res_baseZC.assign(loope->f_zc);
  res_baseCC.assign(loope->f_cc);
  res_baseFz.assign(target_f_z);
  res_valid=(res_foundMult > 1) || res_card;
  return;







/*

  //suppose the center is exact enough
  //2) find bulb base guess c = xc+1/(f_zc+f_zz)   , derivatives at (z=xc,c=xc) (from estimateInterior)
  bulb.bulbe.eval2(period, xc, xc);
  currentWorker->sub(&bulb.bulbe.f_re, xc->re_s);//is 0
  currentWorker->sub(&bulb.bulbe.f_im, xc->im_s);
  currentWorker->add_double(&bulb.bulbe.f_z_re, -1); //is -1

  //s2:=f_zc^2-f_zz*f_cc   always 0 at bulb and card center
  currentWorker->assign(s2_.re_s, f_zz.re_s);
  currentWorker->assign(s2_.im_s, f_zz.im_s);
  s2_.mul(&f_cc_);
  currentWorker->assign(s1.re_s, f_zc.re_s);
  currentWorker->assign(s1.im_s, f_zc.im_s);
  s1.sqr();
  currentWorker->sub(s2_.re_s, s1.re_s);
  currentWorker->sub(s2_.im_s, s1.im_s);
  //g always=0+0i at the bulb center (or it should) and 1+0i at cardioid
  double cardioid_discrim=currentWorker->toDouble(s2_.getMagTmp());
  double precis=currentWorker->toDouble(f_zc.getMagTmp());
  if (cardioid_discrim<precis*1e-10)
    *is_card=false;
  else if (cardioid_discrim<1-precis*1e-10)
    *is_card=false; //?
  else if (cardioid_discrim<1+precis*1e-10)
    *is_card=true;
  else
    *is_card=true; //?



  //taken from estimateInterior: center-base=1/(f_zc+f_zz) ... unless it's a cardioid...
  currentWorker->assign(s1.re_s, f_zz.re_s);
  currentWorker->assign(s1.im_s, f_zz.im_s);
  s1.add(&f_zc);
  double recip_mag=currentWorker->toDouble(s1.getMagTmp());
  //if (recip_mag<1) and (recip_mag>=1e-30) then
  //  recip_mag:=recip_mag; //sure//really?
  if (recip_mag<1e-30) //test for about 1.0 would be enough; it gets larger for small mandels
  {
    *foundMult=1;
    return false;
  };
  s1.recip_prepared();
  currentWorker->assign(cb->re_s, xc->re_s);
  currentWorker->assign(cb->im_s, xc->im_s);
  cb->add(&s1); //cb=xc-1/(f_zc+f_zz), sign is a bit unclear
  if (*foundMult==0)
  {
    currentWorker->assign(rb->re_s, cb->re_s);
    currentWorker->assign(rb->im_s, cb->im_s);
  }
  *foundMult=1;

  //fix r a little
  //move c to estimate of bulb base
  //  and adjust r accordingly
  //repeat
  MandelMath::complex B(currentWorker, &bulb.B_re, &bulb.B_im, true);
  MandelMath::complex C(currentWorker, &bulb.C_re, &bulb.C_im, true);
  MandelMath::complex inte(currentWorker, &interior.inte_re, &interior.inte_im, true);
  for (int cycle=0; cycle<10; cycle++)
  {
    for (int lcycle=0; lcycle<5; lcycle++) //TODO: stop at convergence or just improve everything
    {
      if (!bulb.bulbe.eval2(period, cb, rb))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(f.re_s, rb->re_s);
      currentWorker->sub(f.im_s, rb->im_s);
      currentWorker->add_double(f_z.re_s, -1);
      if (!bulb.lagu.eval(period, &f, &f_z, &f_zz))
      {
        *foundMult=1;
        return false;
      };
      currentWorker->sub(rb->re_s, &bulb.lagu.step_re);
      currentWorker->sub(rb->im_s, &bulb.lagu.step_im);
    }
    int ei=estimateInterior(period, cb, rb);
    if (ei==0)
    {
      *foundMult=1;
      return true;
    }
    else if (ei<=0)
    {
      *foundMult=1;
      return false;
    };
    if (currentWorker->toDouble(&interior.inte_abs)>=4.5)
    {
      *foundMult=1;
      return false;
    };
    //correct step seems to be inte/2 (correct step is guaranteed between inte and inte/4)
    currentWorker->lshift(&interior.inte_re, -1);
    currentWorker->lshift(&interior.inte_im, -1);
    currentWorker->add(cb->re_s, &interior.inte_re);
    currentWorker->add(cb->im_s, &interior.inte_im);

    //f-xc=(cb-xc)*f_c+(rb-xc)*f_z
    //rb-xc=(cb-xc)*f_c+(rb-xc)*f_z
    //(rb-xc)=(cb-xc)*f_c/(1-f_z)   will crash when we reach bulb base because f_z==0

    //rb-xc=(cb-xc)*f_c+(rb-xc)*f_z+(cb-xc)^2*f_cc/2+(rb-xc)^2*f_zz/2+f_zc*(rb-xc)*(cb-xc)
    //-(rb-xc)^2*f_zz/2+(rb-xc)*(1-f_z-f_zc*(cb-xc))=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    //-(rb-xc)^2*f_zz/2+(rb-xc)*(1-f_z-f_zc*(cb-xc))=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    //A=f_zz/2  B=f_z+f_zc*(cb-xc)-1  C=(cb-xc)*f_c+(cb-xc)^2*f_cc/2
    currentWorker->assign(B.re_s, f_zc.re_s);
    currentWorker->assign(B.im_s, f_zc.im_s);
    B.mul(&inte);
    B.add(&f_z);
    currentWorker->add_double(B.re_s, -1);

    currentWorker->assign(C.re_s, f_cc_.re_s);
    currentWorker->assign(C.im_s, f_cc_.im_s);
    currentWorker->lshift(C.re_s, -1);
    currentWorker->lshift(C.im_s, -1);
    C.mul(&inte);
    currentWorker->add(C.re_s, &bulb.bulbe.f_c_re);
    currentWorker->add(C.im_s, &bulb.bulbe.f_c_im);
    C.mul(&inte);

    //let's try in doubles first
    //TODO: pretty much completely wrong
    double r_step_re, r_step_im;
    MandelMath::complex_double_quadratic(&r_step_re, &r_step_im,
                                         currentWorker->toDouble(f_zz.re_s)/2, currentWorker->toDouble(f_zz.im_s)/2,
                                         currentWorker->toDouble(B.re_s)/2, currentWorker->toDouble(B.im_s)/2,
                                         currentWorker->toDouble(C.re_s), currentWorker->toDouble(C.im_s));
    currentWorker->add_double(rb->re_s, -r_step_re);
    currentWorker->add_double(rb->im_s, -r_step_im);
  }
  currentWorker->assign(baseZC->re_s, &bulb.bulbe.f_zc_re);
  currentWorker->assign(baseZC->im_s, &bulb.bulbe.f_zc_im);
  currentWorker->assign(baseCC->re_s, &bulb.bulbe.f_cc_re);
  currentWorker->assign(baseCC->im_s, &bulb.bulbe.f_cc_im);
  return true;
  */
}

template<typename BASE>
MandelEvaluator<BASE>::Bulb::Bulb(MandelMath::number<BASE>::Scratchpad *spad, MandelLoopEvaluator<BASE> *loope):
  loope(loope),
  baseZC(spad), baseCC(spad),
  t1(spad), t2(spad), t3(spad), deltac(spad), deltar(spad),
  first_cb(spad), firstStep(0), dbg_first_rb(spad), prev_rb(spad),
  target_f_z(spad),
  res_cb(spad), res_rb(spad), res_xc(spad),
  res_baseZC(spad), res_baseCC(spad),
  res_card(false), res_foundMult(0), res_valid(false), res_baseFz(spad),
  lagu(spad)
{
  //working_assert(self_allocator.checkFill());
}

#if 0
void MandelEvaluator::Bulb::fixRnearBase(MandelMath::complex *r, const MandelMath::complex *c, int period, int *mult)
{ //"cleverFix" in old code
  //TODO: cycles unused?
  /*MandelMath::complex rb(currentWorker, &bulb.rb_re_, &bulb.rb_im, false);
  MandelMath::complex cb(currentWorker, &bulb.cb_re, &bulb.cb_im, false);
  MandelMath::complex xc(currentWorker, &bulb.xc_re, &bulb.xc_im, false);
  MandelMath::complex baseZC_(currentWorker, &bulb.baseZC_re, &bulb.baseZC_im, false);
  MandelMath::complex baseCC_(currentWorker, &bulb.baseCC_re, &bulb.baseCC_im, false);
  MandelMath::complex s1(currentWorker, &bulb.s1_re, &bulb.s1_im, false);
  MandelMath::complex s2(currentWorker, &bulb.s2_re_, &bulb.s2_im_, false);
  MandelMath::complex cbx_(currentWorker, &bulb.cbx_re, &bulb.cbx_im, false);
  MandelMath::complex rbx_(currentWorker, &bulb.rbx_re, &bulb.rbx_im, false);*/
  if (*mult<=1)
    return;
  rb.assign(r);
  bool is_card=false;
  int foundMult=1;
  bool baseFound=findBulbBase(period, c, &cb, &rb, &xc, &baseZC, &baseCC, &is_card, &foundMult);
  if (!baseFound)
    return;
  if (foundMult<=1)
    return;
  if (*mult!=foundMult)
    *mult=foundMult;
  //s1=c-cb
  s1.assign(c);
  s1.sub(&cb);
  /*currentWorker->assign(rbx.re_s, s1.re_s);
  currentWorker->assign(rbx.im_s, s1.im_s);
  rbx.add(&xc);*/
  //double ratio=currentWorker->toDouble(s1.getMagTmp())/currentWorker->toDouble(cbx.getMagTmp());
  //if (ratio<0.05)
  //s2=s1*baseZC ~= root_z
  s2.assign(&s1);
  s2.mul(&baseZC);
  r->assign(&rb);
  if (!currentWorker->isl0(s2.re) || (*mult==2)) //=0 for baseZC=0 at period=1
  { //above base
    cbx.assign(&cb);
    cbx.sub(&xc);
    rbx.assign(&rb);
    rbx.sub(&xc);
    s2.assign(&cbx);
    s2.recip();
    s2.mul(&rbx);
    double tmpre=s2.getMag_double();
    if (tmpre==0)
    {
      //r:=rb
    }
    else
    {
      double tmpim=std::atan2(currentWorker->toDouble(s2.im), currentWorker->toDouble(s2.re));
      if (tmpim<-M_PI)
        tmpim+=2*M_PI;
      else if (tmpim>=M_PI)
        tmpim-=2*M_PI;
      if (*mult==2)
      {
        tmpre=exp(log(tmpre)/4);
        tmpim/=2;
      }
      else
      {
        tmpre=exp(log(tmpre)/(2*(*mult-1)));
        tmpim/=(*mult-1);
      }
      s2.zero(-tmpre*cos(tmpim), -tmpre*sin(tmpim));
      s2.mul(&rbx);
      //r:=rb+s2
      r->add(&s2);
    }
  }
  else
  { //we're below the bulb base, move close to the base's root
    //dz f_zc+dc/2 f_cc=0
    //dz=-dc/2 f_cc/f_zc
    s1.mul(&baseCC);
    s2.assign(&baseZC);
    s2.recip();
    s2.mul(&s1);
    s2.lshift(-1);
    r->sub(&s2);
  }
}
#endif

//result 0..derivatives or value too large, or other fail (divide by 0)
//result>0 .. tried to return multiplicity but really returns just 1 (1 or >=3) or 2 (mult==2)
template <typename BASE>
int MandelEvaluator<BASE>::laguerre(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> *r, const bool fastHoming) //returns multiplicity...under ideal circumstances
{ //TODO: suggestedMulti = maximumMultip ?
  MandelMath::number<BASE> &newt_tmp2=newt.tmp2;
  double bestfm=1e10; //TODO: actually bestgm? g(z)=f(z)-z
  double prev_grmag=1e10;
  double prev_accy_multiplier=1e10;
  double prev_accy_noise=1e10;
  newt.bestr.assign(*r);
  bool movedOff=false;
  //double accyBound=3e-28/(period*period);
  //was for 80b floats double accyBound2=3e-39*period/log(1+period)*1.5; //1.5=magic
  //double accyBound2=1.23e-32*period/log(1+period)*1.5; //1.5=magic
  double order1; // 1/highest power in the polynomial, 2^period in case of mandelbrot set
  int maxm;
  {
    newtres.first_guess_lagu.assign(*r);
    newtres.first_guess_newt.assign(*r);
    double r_re=r->re.toDouble();
    double r_im=r->im.toDouble();
    newtres.first_fejer_re=r_re; newtres.first_fejer_im=r_im;
    newtres.first_naive1_re=r_re; newtres.first_naive1_im=r_im;
    newtres.first_naive2_re=r_re; newtres.first_naive2_im=r_im;
    newtres.first_naive_re=r_re; newtres.first_naive_im=r_im;
    newtres.naiveChoice=NewtonNaiveChoice::ncClose;
    newtres.first_neumaier1_re=r_re; newtres.first_neumaier1_im=r_im;
    newtres.first_neumaier2_re=r_re; newtres.first_neumaier2_im=r_im;
    newtres.first_lagu1_re=r_re; newtres.first_lagu1_im=r_im;
    newtres.first_lagu1o_re=r_re; newtres.first_lagu1o_im=r_im;
    newtres.firstMu_re=1; newtres.firstMu_im=0; //newtres.firstM=1;
    newtres.firstMum_re=1; newtres.firstMum_im=0;
    newtres.accy_tostop=1;
    newtres.accy_multiplier=1;
    newtres.accy_noise=1;

    //ideally, I want max(abs(r_re), abs(c_re)) then round up to next power of 2
    //but neither ilogb or frexp can do that so I round 1 to 2, 2 to 4
    /*int lor=std::ilogb(r_re); //3->1 2->1 1->0  0.75->-1  0->-max
    if (lor<-2)
      lor=-2;
    r_re=ldexp(2, lor); //1->4  0->2  -2->0.5
    lor=std::ilogb(r_im);
    if (lor<-2)
      lor=-2;
    r_im=ldexp(2, lor);
    r_mag_rough=r_re*r_re+r_im*r_im;*/
  }
  //double accyBound=r_mag_rough*currentWorker->eps2()*period; //eps*sqrt(period) as eps bleeds out// 3e-28/(period*period);
  if (period<5)
  {
    maxm=period+1; //actually for Mandelbrot it's at most p+1 roots nearby   ldexp(1, period-1); //in theory up to n-1 but for Mandelbrot that's rarely the case
    order1=ldexp(1, -period);
  }
  else if (period<1024)
  {
    maxm=15;
    order1=ldexp(1, -period);
  }
  else
  {
    maxm=15;
    order1=0;
  }
  //maxm=1;
  //int multiplicity1=1;
  int lastm=1;
  double lastm_err=0;
  int prevm=1;
  bool last_was_newton=true;
  //double prevm_err=10;
  bool triedZeroGzrm=false;
  /*struct //instead, implement the last case in periodCheck()
  {
    bool didfix;
    int mult;
  } clever; //improve accuracy around point where 2 bulbs touch
  clever.didfix=false;
  if (suggestedMultiplicity>1)
    clever.mult=suggestedMultiplicity;
  else
    clever.mult=1;*/
  newt.prevR.assign(*r); //preferably .zero(infinity, 0) but now we do if (newtonCycle>0)
  //prev_gz_mag=inf, prev_g_mag=inf
  for (int newtonCycle=0; newtonCycle<50; newtonCycle++)
  {
    newtres.cyclesNeeded=newtonCycle+1;
    if ((movedOff) && (newtonCycle>10) && (order1>=0))
    {                                    //  p m -> p
      order1=-1;                         //  2 2    1
      bestfm=1e10;                       //  4 3    2
      //multiplicity1=1;                   //  4 5    1
    };
    //TODO: can we skip computing fzz_r if order1<0? and remember last valid multiplicity or set it to 1
    //always half of eps_cumul10   double eps_cumul05=0.5;
    //double fc_re=0, fc_im=0;
    if (newtonCycle>0)
    {
      if (!loope.eval_zz(period, c, *r, true))
        return 0;
    }
    else if (currentParams.nth_fz<=0 || currentParams.nth_fz>=period)
    {
      if (!loope.eval_zz(period, c, *r, true))
        return 0;
      newtres.nth_fz.assign(loope.f_z);
      newtres.nth_fz.re.add_double(1);
    }
    else
    {
      if (!loope.eval_zz(currentParams.nth_fz, c, *r, false))
        return 0;
      newtres.nth_fz.assign(loope.f_z);
      if (!loope.eval_zz(period-currentParams.nth_fz, c, *r, true, false))
        return 0;
    }
    newt.f_r.assign(loope.f);
    newtres.fz_r.assign(loope.f_z);
    newt.fzz_r.assign(loope.f_zz);
    //g(r)=f(r)-r, gz(r)=fz(r)-1
    //newt.f_r.sub(r);
    //currentWorker->add_double(newtres.fz_r.re, -1);
    double g_r_mag=newt.f_r.getMag_double();
    double gz_r_mag=newtres.fz_r.getMag_double();
    //newtres.accy_tostop=eps_cumul;//r_mag_rough*eps_cumul10;
    newtres.accy_tostop=std::max(1.0, g_r_mag/c.eps2());
    newtres.accy_multiplier=std::max(1.0, 1/gz_r_mag);
    double f_z_mag=MandelMath::sqr_double(1+loope.f_z.re.toDouble())+
                   MandelMath::sqr_double(loope.f_z.im.toDouble());
    newtres.accy_noise=period*f_z_mag/
        loope.f_z_minmag.toDouble(); //the error at minimum f_z is magnified to f_z_mag
                                  //the process is continuous but
                                  //a reasonable conservative estimate of error at that time is sqrt(period)*eps
#if CLEVER_FIX
//c=-0.7499 p=2
//  ideally, r=-0.5+-0.01i who are repelling  (and +0.5+-sqrt(0.9999) who are repel and attr)
//  but we have f(-0.5001)=-0.49979999, f^2(-0.5001)=-0.5000999699959999
//  due to rounding errors, it looks as if we are at a root
//  and this point is attracting, so we have verified a false double period
//  but also c=0.25+0.5i p=4 r=0.5i: f(0.5001i)=-0.00010001+0.5i
//                                           f2=0.000000100020001+0.49989999i
//                                           f3=0.00009999999800000004000600040001+0.500000009999999499939998i
//                                           f4=6.00240023997099179944002000640056002800080001 × 10^-20 + 0.50009999999999999989999400280036997598939871991799719996... i
//  and this point seems to be repelling
//  there's really no way around this using finite precision
//  so we need something CLEVER
    if (!clever.didfix &&
        (g_r_mag<1e-16) &&
        (((gz_r_mag<5e-3) && !currentWorker->isle0(fz_r.re_s)) || //in bulb close to its base and at the wrong root
         (gz_r_mag<1e-9) || //so close to the base we don't know which root we have
         ((period==2) && !currentWorker->isle0(fz_r.re_s) && (currentWorker->toDouble(fz_r.re_s)<0.14)))) //we skip check for period=1 so special check for the point of attachment of bulb 1/2
    {
      clever.didfix=true;
      fixRnearBase(r, c, 0, period, &clever.mult);
      continue;
    };
#endif
    if (g_r_mag==0)  //7e-33..4e-40 does not need more; much..5e-38 needs more
    { //r is good enough already      (f_c.re*f_c.re+f_c.im*f_c.im)/(f_zc.re*f_zc.re+f_zc.im*f_zc.im)
      return lastm;
    };
    if (gz_r_mag==0)
    {
      if (triedZeroGzrm)
        return 0;
      triedZeroGzrm=true;
    }
    else if (newtonCycle>0)
    {
      //new conditions
      //the one legit reason to end: step<2^-53/|f'| (for |f'|<1) exactly because step*|f'|=2^-53
      //    |f|/|f'|<2^-53/|f'|
      //    |f|^2<2^-106=1.23e-32
      //if (g_r_mag<1.0*r_mag_rough*currentWorker->eps2()) //maybe up to (1+gz_r_mag)*r_mag*eps2*log2(period)
      //if (g_r_mag<eps_cumul10*r_mag_rough*currentWorker->eps2()) //maybe up to (1+gz_r_mag)*r_mag*eps2*log2(period)

      bool cond1=prevm==lastm && lastm_err<0.125;//we don't *really* need prevm_err small   && prevm_err<0.125;
      //even better: once step<eps^(3/4), choose between previous and current root and that's it
      double step_limit;
      if (lastm_err<=0 || //==0 really, just in case
          lastm!=1 || !cond1 || period<=1) //no need to compute step_limit
        step_limit=c.eps234();
      else if (last_was_newton)
      {
        step_limit=c.eps2()/(lastm_err)*4; //rel. error after Newton step ~ sqrt(lastm_err)/2  (for lastm==1)
        if (step_limit>c.eps234())
          step_limit=c.eps234();
        double step=newt.newtX.getMag_double();
        if ((step<step_limit) && !(newt.newtX.getMag_double()<c.eps234()))
          nop(); //we stop now and continued before
        else if (!(step<step_limit) && (newt.newtX.getMag_double()<c.eps234()))
          nop(); //we continue now and stopped before
      }
      else
      {
        step_limit=c.eps2()/MandelMath::sqr_double(lastm_err*4); //rel. error after Ostrowski step ~ lastm_err/4
        if (step_limit>c.eps234())
          step_limit=c.eps234();
        double step=newt.newtX.getMag_double();
        if ((step<step_limit) && !(newt.newtX.getMag_double()<c.eps234()))
          nop(); //we stop now and continued before
        else if (!(step<step_limit) && (newt.newtX.getMag_double()<c.eps234()))
          nop(); //we continue now and stopped before
      }
      bool stepgood=lastm==1 && cond1 &&
                    (newt.newtX.getMag_double()<step_limit ||
                     period<=1 || //with some luck, we did a Laguerre step and that hits the root in one step (period=1 .. degree=2)
                     //when gz_r_mag<eps^(1/4), we won't make steps smaller than eps234
                     //newt.newtX.getMag_double()<currentWorker->eps234()*newtres.accy_multiplier || //covers newt.newtX.getMag_double()<currentWorker->eps234()
                     prev_grmag<2*c.eps2()); //hopefully better than currentWorker->eps234()*newtres.accy_multiplier
      //for multiple root, step never gets small, but value does
      //bool ggood=prev_grmag<currentWorker->eps234();
      //if prevm==2 we arrive at the middle of the 2 roots and we're not quite done yet
      //for prevm>2 there is actually a root in the middle so all is good
      //before unleashing the full power of prevm>lastm, lemme try
      bool ggood=lastm>1 && cond1 && prev_grmag<c.eps234();
      //sometimes lastm says 3 but in reality it's 1; crosscheck by f/f' small (ideally just f small && f' small, but what limit?)
      //since lastm>1, we know f'' is not small; since lastm is of order ~1, we know f/f' ~ f'/f''
      //so check f'/f'' instead of f/f'
      bool ggood2=gz_r_mag*(lastm*lastm)<loope.f_zz.getMag_double()*c.eps234();
      if (ggood && !ggood2)
        ggood=false; //quite rare but helps
      if (newtonCycle>0 && (stepgood || ggood))
      {
        if (g_r_mag>prev_grmag)
        {
          r->assign(newt.prevR);
          newtres.accy_tostop=std::max(1.0, prev_grmag/c.eps2());
          newtres.accy_multiplier=prev_accy_multiplier;
          newtres.accy_noise=prev_accy_noise;
          newtres.fz_r.assign(newt.prevGz);
        };
        //if (!stepgood) ... ? set accy_multiplier?
        if (newtres.accy_multiplier>100)
        { //we care mostly about not letting infinity into accy_multiplier
          //max error=eps2^(1/4) -> max multiplier=eps2^(3/4)
          double max_accy_multiplier=1/c.eps234();
          if (newtres.accy_multiplier>max_accy_multiplier)
            newtres.accy_multiplier=max_accy_multiplier;
        };
        return lastm;
      }
    };
    /*if (r->isequal(&newt.bestr) &&
        (bestfm<1e10) && (newtonCycle>2)) //Lagu can cycle in first 2 cycles
    { //Laguerre can cycle (at least my version), e.g. per=2, c=-0.6640625-0.015625i, r=-0.614213552325963974-0,0179149806499481201i
      if (g_r_mag<6*eps_cumul10*r_mag_rough*currentWorker->eps2()) //should be tested above but maybe use different margin here?
        return lastm;
      return 0; //just fail and try again next time
    };*/
    if (g_r_mag<bestfm)
    {
      bestfm=g_r_mag;
      newt.bestr.assign(*r);
    };
    if (newtonCycle>0 && g_r_mag>=prev_grmag)//+1000*newtres.accy_tostop*currentWorker->eps2())
    { //g_r_mag didn't go down and it wasn't because of limited precision
      //no reason to continue from here, polynomials don't have poles
      //undo half of last step and try again
      newt.newtX.lshift(-1);
      //r->add(&newt.newtX);
      r->assign(newt.prevR); //Justin Case
      r->sub(newt.newtX);
      if (r->isequal(newt.prevR))
      {
        newt.newtX.zero();
        //prevm_err=0;
        lastm_err=0;
        prevm=lastm;
      }
      continue;
    };
    prev_grmag=g_r_mag; //almost like bestfm?
    prev_accy_multiplier=newtres.accy_multiplier;
    prev_accy_noise=newtres.accy_noise;
    newt.prevGz.assign(newtres.fz_r);

    /*
    see also: Gaston H. Gonne: A Study of Iteration Formulas for Root Finding,
              Where Mathematics, Computer Algebra and Software Engineering Meet
    https://www.math.uni-bielefeld.de/~rehmann/ECM/cdrom/3ecm/pdfs/pant3/gonnet.pdf
    eta=1-H/G^2

    see also jenkins-traub: https://github.com/jervisfm/JenkinsTraub/blob/master/poly.py
    (not used at all, requires manipulation of coefficients)

    derive Laguerre's method, multiplicity m!=1, order of poly=n
    assume roots A and B, distance a=z-A away, others b=z-B away
    f(z)=C (z-A)^m (z-B)^(n-m)
    take ln, diff twice
    ln f(z) = ln C + m*ln(z-A) + (n-m)*ln(z-B)
    d/dz ln f(z) = d/dz ln C + d/dz m*ln a + d/dz (n-m)*ln b
    G = f'(z)/f(z) = m/a + (n-m)/b
    d/dz (f'(z)/f(z)) = d/dz m/a + d/dz (n-m)/b     a=z-A, (1/a)'=-1/a^2
    -H = (f''(z)*f(z) - f'(z)*f'(z)) / f(z)^2 = - m/a^2 - (n-m)/b^2
    G = m/a + (n-m)/b
    H = m/a^2 + (n-m)/b^2
    solve for aa=1/a from G=f'/f, H=G^2-f''/f
    (G-m*aa)=(n-m)/b
    H*(n-m) = m*(n-m)*aa^2 + (G-m*aa)^2
    0 = m*n*aa^2 - 2*m*G*aa + G^2-H*(n-m)
    0 = aa^2 - 2*G/n*aa + (G^2-H*(n-m))/(m*n)
    aa*n = G +- sqrt( (n/m-1)*(n*H-G^2) )
    Newton's step is f/f' = 1/G
    Laguerre's step is a=1/aa=n/(n*aa)=n/(G +- sqrt( (n/m-1)*(n*H-G^2) ))
    a=1/G/(1/n +- sqrt( (1/m-1/n)*(H/G^2-1/n) ))
    H/G^2=1-f''*f/f'^2 = 1/M
>>  a=f/f'/(1/n +- sqrt( (1/m-1/n)*(1/M-1/n) ))
      but fails if f'=0
    a=f/ffix'
    ffix'=f'/n +- sqrt( (1/m-1/n)*(f'^2-f''*f-f'^2/n) )         M=f'^2/(f'^2-f''*f)
      only fails if f'=f''=0

    when f'=0, G=0, M=inf; ideal m=1
    G=f'/f   H=G^2-f''/f   M=1/(1-f''*f/f'^2)   Re(M) rounds to 1 except for f''*f/f'^2 in circle (c=2/3 r=1/3)
>>  a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))    M=G^2/(G^2-f''/f) (G can be 0, can't divide)
    x^2+1 @ 0: f=1 f'=0 f''=2 G=0 H=-2
    +-i=1/(+- sqrt( (1-2/m) ))  m=1

    let's not talk about the derivation but
    either f'^2-f''*f=a+bi    f'^2=c+di
    or     H=f'^2/f^2-f''/f=a+bi    G^2=f'^2/f^2=c+di
    then   m=Round( (a*c+b*d)/(a*a+b*b) ) = Re(G^2*H^T)/mag(H) = Re(G^2/H)   x^T = conj(x)



    for m=n
    G = m/a
    a=n/G
    same as above even though the derivation is invalid

    when m:=M, a=f/f'/(1-f''*f/f'^2)=f/f'*M

    when laguH=0; ideal m=1
    a=1/(G/n +- sqrt( (1/m-1/n)*(-G^2/n) ))
    a=1/G/(1/n +- i*sqrt( (1/m-1/n)*(1/n) )) = n/G/(1 +- i*sqrt(n/m-1))
    x^2+x+1/2  f(0)=1/2 f'=1 f''=2  1*1-2/2=0 G=2 f''/f=4     -1/2+-i/2
    a=1/(1 +- sqrt( (1-2/m) ))     m=1
    -1/2+-i/2=1/(1 +- i*sqrt(2/m-1))
    1=m

    (x-1)^3*(x+1)^2  : H=0 at x = -1/5 +- (2 i sqrt(6))/5
    ReplaceAll [{(x - 1)^3*(x + 1)^2, D[(x - 1)^3*(x + 1)^2, {x}], D[(x - 1)^3*(x + 1)^2, {x, 2}]}, {x->-1/5 + (2*I*Sqrt[6])/5}]
    {-17856/3125 + (2112 i sqrt(6))/3125, 864/125 + (672 i sqrt(6))/125, 144/5 - (48 i sqrt(6))/5}
    f*f''=(-17856/3125 + (2112 i sqrt(6))/3125)*(144/5 - (48 i sqrt(6))/5) = -1963008/15625 + (1161216 i sqrt(6))/15625
    f'^2 = -1963008/15625 + (1161216 i sqrt(6))/15625
    ->a=4/5+(2 I sqrt(6))/5 = 5/(864/125 + (672 i sqrt(6))/125)/(1 +- I*sqrt(5/m-1))
    solve 4/5+(2*I*Sqrt(6))/5 = 5*(-17856/3125 + (2112 i sqrt(6))/3125)/(864/125 + (672*I*Sqrt(6))/125)/(1 + I*Sqrt(5/x-1))
      m=x=2, no solutions for 1-I*Sqrt


    when f''=0, f'!=0
    a=f/f'/(1/n +- sqrt( (1/m-1/n)*(1-1/n) ))
    x+1=0 f''=0 f'=1 f=1 G=1
    1=1/1/(1/1 +- sqrt( (1/m-1/1)*(1-1/1) ))    any m

    ostrowski-1 (method 3 from [Gonnet]): delta/sqrt(1-eta) = sqrt(1/H) = sqrt(m)/G       H=f'^2/f^2-f''/f
    ostrowski-n (method 22 from [Gonnet]): 1/G/sqrt(1-eta-eta^2/(n-1)) = sqrt(1/H)/sqrt(1-(m+1/m-2)/(n-1))
      blows if H=0 or m+1/m=n+1: for n=1, m=1; for n=2, m=2.62; for n=3, m=3.73 so should never blow
      [petkovic] says choose sqrt that sqrt(H) is closer to G/m, means sqrt(1/H) is closer to m/G
      also that x2:=x-sign(f*f')/sqrt(H) at https://miodragpetkovic.com/publikacije/on-some-improvements-of-square-root-iteration-for-polynomial-complex-zeros/ page 1
    https://interval.louisiana.edu/reliable-computing-journal/volume-16/reliable-computing-16-pp-225-238.pdf
      Ostrowski-Like Method for the Inclusion of a Single Complex Polynomial Zero∗
      Mimica R. Milošević, Miodrag S. Petković
      eq(14) can't figure out what they mean, looks like division by 0 to me
      ...          sqrt(m)/sqrt(H-(N-m)/(x-x)^2)
      from Gonnet: sqrt(m/H)/sqrt(1+(n-m)(m-1)/(n-1)) = (G/H)/sqrt(1+(n-m)(m-1)/(n-1))
    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.723.1400&rep=rep1&type=pdf
    Beny Neta, Changbum Chun: On a family of Laguerre methods to find multiple roots of nonlinear equations
      lagu-m: 2m/G/(1+-sqrt(2m-1-2mff''/f'^2)) "Euler-Cauchy"
              2/(G/m+-sqrt(2H-G^2/m)/sqrt(m))   H=G^2/M
              2sqrt(1/H)/(1+1/sqrt(m)) if m=M   ; for m=10 about 1/2 step of ostro-m
      ostro-m: sqrt(m)/sqrt(H)
        notice that my lagu-m for n->inf becomes ostrowski-m

    CLAM formula z1=z0−(N*f/f')*(1−Q^(m/N))/(1−Q), where Q=[N*(1−f*f''/f'^2)−1]/(N/m −1)
        the formula has never been found to fail for general polynomials
        https://link.springer.com/chapter/10.1007/3-540-62598-4_83?noAccess=true
        Tien Chi Chen: Convergence in iterative polynomial root-finding
        z0−(N*f/f')*(1−Q^(m/N))*(1/m−1/N)/[1/m-(1−f*f''/f'^2)]   Q=[N*(1−f*f''/f'^2)−1]/(N/m −1)
        1-f''*f/f'^2 = 1/M
           z0−(n*f/f')*(1/m−1/n)/(1/m-1/M)*(1−Q^(m/n))   Q=[1/M−1/n]/(1/m−1/n)
      lagu z0-f/f'/(1/n +- sqrt( (1/m-1/n)*(1/M-1/n) ))
      lagu z0-n*f/f'/(1 +- sqrt( n^2*(1/m-1/n)^2*(1/M-1/n)/(1/m-1/n) ))
      lagu z0-(n*f/f')/(1 +- n*(1/m-1/n)*sqrt( Q ))
      clam z0−(n*f/f')*(1/m−1/n)*(1/(1/m-1/M)−(Q/(1/m-1/M)^(n/m))^(m/n))   Q=(1/M−1/n)/(1/m−1/n)  M=k*m  1/m-1/M=1/m-1/k/m=1/m*(1-1/k)   Q=(1/k-1)*(1/m)/(1/m−1/n)+1=(1/M-1/m)/(1/m−1/n)+1
      clam z0−(n*f/f')*(1/m−1/n)/(1/m-1/M)*(1−(Q)^(m/n))   Q=1-(1/m-1/M)/(1/m−1/n)   Q=1-q
      clam z0−(n*f/f')*(1−(1-q)^(m/n))/q    q=(1/m-1/M)/(1/m−1/n)
             C[q]=(1−(1-q)^(m/n))/q   C[0]~m/n C[1]~infinity  C[q]~m/n+1/2*m/n*(1-m/n)*q+1/6*m/n*(1-m/n)*(2-m/n)*q^2+...
             C[q]~m/n*-log(1-q)/q for small m/n  (untested)
      clam z0−(n*f/f')*(1−(1-q)^(m/n))/q    q=(1/m-1/M)/(1/m−1/n)
          clam z0−(m*f/f')*(1−(1-q)^a)/q/a             a=m/n
      lagu z0-(n*f/f')/(1 +- n*(1/m-1/n)*sqrt( 1-q ))
          lagu z0-(m*f/f')/(a +- (1-a)*sqrt( 1-q ))    a=m/n
          clam and lagu same for a=0.5 and a=1
          plot [(1−(1-x)^a)/x/a; 1/(a + (1-a)*sqrt( 1-x ))] for a=0.25 for x from 0 to 1
            similar but not the same, very similar for 0<=q<=0.5
      clam for small q: ~z0−(n*f/f')*m/n=z0−m*f/f'
      lagu for small q: ~z0-(f/f')/(1/m)=z0-m*f/f'
      clam z0−(m*f/f')*(1+1/2*(1-m/n)*q+1/6*(1-m/n)*(2-m/n)*q^2+1/24*(1-m/n)*(2-m/n)*(3-m/n)*q^3)
      lagu z0-(m*f/f')/(1+(1-m/n)*( - q/2 - q^2/8 - q^3/16 ))
         1/(1+(1-a)*(sqrt(1-x)-1))=1+1/2 (1 - a) x + 1/8 (2 a^2 - 5 a + 3) x^2 + 1/16 (-2 a^3 + 8 a^2 - 11 a + 5) x^3 + ...
      lagu z0-(m*f/f')*(1+1/2*(1-m/n)*q+1/8*(2*(m/n)^2-5*m/n+3)*q^2 + 1/16 (-2 a^3 + 8 a^2 - 11 a + 5) x^3 + ...)
         (1-m/n)*(2-m/n)=(m/n)^2-3*m/n/2+2  1/8*6/8*8/6=1/6*6/8
      lagu z0-(m*f/f')*(1+1/2*(1-m/n)*q+1/6*(6/4*(m/n)^2-15/4*m/n+9/4)*q^2 + 1/16 (-2 a^3 + 8 a^2 - 11 a + 5) x^3 + ...)





    ----- how to find m ----
    simple: f=x^m
    f=x^m f'=m x^(m-1)  f''=m(m-1) x^(m-2)
    f/x^(m-2)=x^2  f'/x^(m-2)=m x  f''/x^(m-2)=m(m-1)
    f''*f/f'^2=m(m-1) x^(m-2) x^m / m/m / x^(m-1)/ x^(m-1) = (m-1)/m
    // f'/f=m/x    f''/f=m(m-1)/x^2
    // f''/f/f'*f=f''/f'=(m-1)/x
    1/m=1-f''*f/f'^2

    (x-1)^3*(x+1)^2 at 0.99
    f=-3.9601×10^-6  f'=0.00118405  f''=-0.23522
    f''*f/f'^2=0.6644  1/(1-...)=2.98   1/m=1-f''*f/f'^2

    full:
    f=(x-a)^m (x-b)^(n-m)
    f'=(x-a)^m (x-b)^(n-m)= m (x-a)^(m-1) (x-b)^(n-m) + (n-m) (x-a)^m (x-b)^(n-m-1)
    f''=(m-1) m (x-a)^(m-2) (x-b)^(n-m) + 2 m (n-m) (x-a)^(m-1) (x-b)^(n-m-1) + (n-m-1) (n-m) (x-a)^m (x-b)^(n-m-2)
      limit of 1/(1-D[(x-a)^m (x-b)^(n-m),{x,2}]*(x-a)^m (x-b)^(n-m)/D[(x-a)^m(x-b)^(n-m),{x,1}]^2) as b goes to infinity
        ->m
    w.l.o.g. x=0   Z=a/b
      M=1/(1-f''*f/f'^2)=(b m + a (n-m))^2/(b^2 m + a^2 (n-m))=(m + Z (n-m))^2/(m + Z^2 (n-m)) ~ m + 2*Z (n-m)
      from afar (a~b): m=n

    find m such that newtX=laguX
    m/G=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
    G/m=G/n + sqrt( (1/m-1/n)*(H-G^2/n) )
    1/m = (H-G^2/n)/G^2+1/n
    1/m = H/G^2


    find b (bb) from Lagu:
    G = m/a + (n-m)/b
    H = m/a^2 + (n-m)/b^2
    solve for bb=1/b from G=f'/f, H=G^2-f''/f
    (G-(n-m)*bb)^2/m = m/a^2
    H = m/a^2 + (n-m)*bb^2
    H = (G-(n-m)*bb)^2/m + (n-m)*bb^2
    m*H = G*G-2*G*(n-m)*bb+(n-m)*(n)*bb*bb
    (n-m)*(n)*bb*bb-2*G*(n-m)*bb+(G*G-m*H) = 0
    bb1,2 = G/n*(1+-sqrt(1-1/(n-m)*(n)*(1-m*H/G^2)))
    b1,2 = n/G/(1+-sqrt(1-1/(n-m)*(n)*(1-m*H/G^2)))   H/G^2=1/M
    b1,2 = n/G/(1+-sqrt(1+n*(m/M-1)/(n-m))
    a = 1/G/(1/n + sqrt( (1/m-1/n)*(1/M-1/n) ))
    a/b=(1/n+-sqrt(1/n^2+(m/M-1)/n/(n-m))/(1/n + sqrt( (1/m-1/n)*(1/M-1/n) ))
    a/b=(1+-sqrt(1+(m/M-1)*n/(n-m))/(1 + sqrt( (n/m-1)*(n/M-1) ))
    a/b=(1+-sqrt(1+(m/M-1)*n/(n-m)))/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
    m ~ (M-2*n*a/b)/(1-2*a/b)  for small a/b
    a/b=(1-sqrt(1+(m/M-1)*n/(n-m)))/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
    a/b~-(m/M-1)*n/(n-m)/2/(1 + sqrt( 1+((n-m)/M-1)*n/m )) = -(m/M-1)*n/(n-m)/2/Q
    m=M+(m/M-1)*n/(1 + sqrt( 1+((n-m)/M-1)*n/m ))
      -> m:=M
      not gonna work
    a/b @ m=1: =(1-sqrt(1+(1/M-1)*n/(n-1)))/(1 + sqrt( 1+((n-1)/M-1)*n/1 ))
      m ~ (M-2*n*a/b)/(1-2*a/b)
      -> funny function that only depends on M and n; for n=5 (x=M):
        plot (x(1 + Sqrt[-(((x - 5) (-1 + 5))/x)]) + 2 5 (-1 + Sqrt[(x - 5)/(x - x 5)]))/(-1 + Sqrt[-(((x - 5) (-1 + 5))/x)] + 2 Sqrt[(x - 5)/(x - x 5)])




*/

    //1/f should be fine, or we'd be at the root
    newt.tmp1.assign(newt.f_r);
    newt.tmp1.recip();    //1/f
    newt.laguG.assign(newtres.fz_r);
    newt.laguG.mul(newt.tmp1); //laguG = f'/f
    newt.fzzf.assign(newt.fzz_r);
    newt.fzzf.mul(newt.tmp1); //f''/f

    // laguH=laguG^2-fzzf
    // m=Round( Re(G^2*H^T)/mag(H) )
    // a=1/(G/n +- sqrt( (1/m-1/n)*(G^2-fzzf-G^2/n) ))
    newt.laguG2.assign(newt.laguG);
    newt.laguG2.sqr();    //G^2
    newt.laguH.assign(newt.laguG2);
    newt.laguH.sub(newt.fzzf); //G^2-fzzf = f'^2/f^2-f''*f/f^2 = laguH
    double mum_re=1, mum_im=0;
    {
      /*double G2HT_re=currentWorker->toDouble(laguG2.mulreT(&laguH));
      double H_mag=currentWorker->toDouble(laguH.getMagTmp());
      //turns out that if mu=m then mu=m=G^2/H
      //1.5*mag(H)>Re(G^2*H^T) ... m=1
      //300*mag(H)<Re(G^2*H^T) ... m=300
      //mag(H)<Re(G^2*H^T)*order1 ... m=1/order1
      if (1.5*H_mag>=G2HT_re) //>= so we don't divide G^2/H if H=G=0
        m=1;
      else if ((clever.mult+0.5)*H_mag<=G2HT_re)
        m=1; //best practice is to use m=1 if H=0   clever.mult;
      else if (H_mag*(maxm-0.5)<G2HT_re)
        m=maxm;
      else
        m=qRound(G2HT_re/H_mag);*/

      //solve for m=mu:   m=G^2/H
      //double mum_re=1, mum_im=0; //better than mu? yes
      double h_re=newt.laguH.re.toDouble();
      double h_im=newt.laguH.im.toDouble();
      double h_mag=h_re*h_re+h_im*h_im;
      double g2_re=newt.laguG2.re.toDouble();
      double g2_im=newt.laguG2.im.toDouble();
      if (h_mag>0.01)
      { //h_mag ok
        mum_re=(g2_re*h_re+g2_im*h_im)/h_mag;
        mum_im=(g2_im*h_re-g2_re*h_im)/h_mag;
      };
      if (newtonCycle==0)
      {
        newtres.firstMum_re=mum_re;
        newtres.firstMum_im=mum_im;
      };
    }

  //can't quite remember what all this mu thing means :-(
    //seems to give ~sqrt(m)
  //m= some func of mu where mu is solution of ((1-1/n)*H/G^2-1/n) mu^2 + 2*mu/n -1=0
  //with m as input:                           ((1-m/n)*H/G^2-1/n) mu^2 + m/n 2*mu -m = 0
  double mu_re, mu_im;
  double G2_mag=newt.laguG2.getMag_double();
  if (G2_mag<0.01)
  { //G2_mag bad
    mu_re=1; mu_im=0;
    if (newtonCycle==0)
    {
      newtres.firstMu_re=1;
      newtres.firstMu_im=0;
    };
  }
  else
  {
    newt.laguX.assign(newt.laguG2);
    newt.laguX.im.chs();
    newt.laguX.mul(newt.laguH);
    double a_re=newt.laguX.re.toDouble()/G2_mag*(1-order1)-order1;
    double a_im=newt.laguX.im.toDouble()/G2_mag*(1-order1);
    //double mu_re, mu_im;
    MandelMath::complex_double_quadratic(&mu_re, &mu_im, a_re, a_im, order1, 0, -1, 0);
    if (newtonCycle==0)
    {
      newtres.firstMu_re=mu_re;
      newtres.firstMu_im=mu_im;
    };
  }
  int m=1;
#if 0 //0..use mum, 1..use mu
  if (!(mu_re>=1.3)) //also m=1 if mu_re is NaN    (mu_re<1.3)
    m=1;
  else {/*if (abs(mu_im)>mu_re/2)
    m=1;
  else
  {
    double mu_mag=mu_re*mu_re+mu_im*mu_im;
    m=qRound(sqrt(mu_mag)); //or just round mu_re?
    */
    m=qRound(mu_re);
    if (m>maxm)
      m=maxm;
  }
#else
  //using lastm_err instead if (!(mum_re>=0.5)) //also m=1 if mum_re is NaN    (mum_re<0.5)
  //  m=0;
  if (!(mum_re>=1.3)) //also m=1 if mum_re is NaN    (mum_re<1.3)
    m=1;
  else {/*if (abs(mum_im)>mum_re/2)
    m=1;
  else {*/
    m=qRound(mum_re); //some say qRound(sqrt(mum_re*mum_re+mum_im*mum_im))
    //TODO: try m=qRound(sqrt(mum_re*mum_re+mum_im*mum_im));
    //  does it have singularities at the same places?
    if (m>maxm)
      m=maxm;
  }
#endif

  if (newtonCycle==0)
      {
        //Fejer bound: smaller solution x of
        //fzz/(n-1) x^2+2 fz x + n f=0
        //x=y*n
        //fzz*n/(n-1) y^2+2 fz y + f=0

        double r_re=r->re.toDouble();
        double r_im=r->im.toDouble();
        //numbers are small but don't need precision so let's do it in double
        double a_re=newt.fzz_r.re.toDouble()/(1-order1);
        double a_im=newt.fzz_r.im.toDouble()/(1-order1);
        double fz_re=newtres.fz_r.re.toDouble();
        double fz_im=newtres.fz_r.im.toDouble();
        double f_re=newt.f_r.re.toDouble();
        double f_im=newt.f_r.im.toDouble();
        MandelMath::complex_double_quadratic(
              &newtres.first_fejer_re, &newtres.first_fejer_im,
              a_re, a_im,
              fz_re, fz_im,
              f_re, f_im);
        newtres.first_fejer_re=r_re+ldexp(newtres.first_fejer_re, period);
        newtres.first_fejer_im=r_im+ldexp(newtres.first_fejer_im, period);

        //Batra's bound https://www.tuhh.de/ti3/paper/rump/Ru03c.pdf theorem 3.8
          //but only for real coefficients
        //|fz r|-|f + fzz/2 r^2|=0, find r
        //sqrt(fz fz^T) r=sqrt((f + fzz/2 r^2)(f^T + fzz^T/2 r^2))
        //sqrt(fz fz^T) r=sqrt((|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4))
        //(a+bi)(c-di)+(a-bi)(c+di)=2ac+2bd=2 Re(f fzz^T)
        //|fz|^2 r^2=|f|^2+ Re(f^T fzz) r^2 + |fzz|^2/4 r^4
        //0=|f|^2+ (Re(f^T fzz)-|fz|^2) rr + |fzz|^2/4 rr^2    r=sqrt(rr)

        /*MandelMath::complex_double_quadratic(&newtres.first_batra, &a_im,
            currentWorker->toDouble(fzz_r.getMagTmp())/4, 0,
            (currentWorker->toDouble(f_r.mulreT(&fzz_r))-currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
            currentWorker->toDouble(f_r.getMagTmp()), 0);
        if (newtres.first_batra>=0)
          newtres.first_batra=sqrt(newtres.first_batra);*/

        //https://ur.booksc.eu/book/5736333/a5b588
        //ZAMM - Journal of Applied Mathematics and Mechanics / Zeitschrift für Angewandte Mathematik und Mechanik
        //1988 Vol. 68; Iss. 6
        //Dr. A. Neumaier: An Existence Test for Root Clusters and Multiple Roots
        //fi(c, r, alpha)=r abs(re((f(c+r e^ialpha)-f(c))/(c+r e^ialpha)))-abs(f(c))
        //  addition from https://ur.booksc.eu/book/5736333/a5b588 remark 3:
        //  f needs to be divided (or rotated) by f' first to make f' real
        //for all alpha, which r makes fi==0 ?
        //abs(re(f'*r+f''/2 r^2 e^ialpha))=abs(f)
        //for max re(f'*r+f''/2 r^2 e^ialpha), we need max re(f'+f''/2 r e^ialpha) because r is real
        //f'' e^ialpha=real
        //e^ialpha=f''^T/sqrt(f'' f''^T)=sqrt(f''^T/f'')
        //abs(re(f'*r+ r^2/2 sqrt(f'' f''^T)))-abs(f)=0
        //r*abs(re(f'))+ r^2/2 sqrt(f'' f''^T)-abs(f)=0
        /*if (currentWorker->isle0(fz_r.re_s))
          MandelMath::complex_double_quadratic(&newtres.first_batra, &a_im,
              +sqrt(currentWorker->toDouble(fzz_r.getMagTmp()))/2, 0,
              //currentWorker->toDouble(fz_r.re_s)/2, 0,
              sqrt(currentWorker->toDouble(fz_r.getMagTmp()))/2, 0,
              +sqrt(currentWorker->toDouble(f_r.getMagTmp())), 0);
        else*/
        MandelMath::complex_double_quadratic(&newtres.first_neumaier1_re, &newtres.first_neumaier1_im,
            -sqrt(newt.fzz_r.getMag_double())/2, 0,
            sqrt(newtres.fz_r.getMag_double())/2, 0,
            -sqrt(newt.f_r.getMag_double()), 0);

        /* Neumaier for k=2:
        Re(f''(z)/2) > |f| r^-2 + |f'| r^-1      r real, r>0
        f''(z)=f''+(z-z0)f'''+...
        |f''+r|f'''|+...|/2 r^2 > |f| + |f'| r
        ...+|f'''|r^3/2+|f''| r^2/2 > |f| + |f'|r
                        |f''| r^2/2 > |f| + |f'|r
        gives always r<0 but that's the wrong root
        |f''| r^2/2 - |f'|r - |f| =0
        ( |f'|+-sqrt(|f'|^2+4 |f''|/2*|f|) )/|f''|
        |f'|/|f''|+-sqrt(|f'|^2/|f''|^2+2*|f|/|f''|)
        works if |f'|^2/|f''|+2*|f|>0 i.e. always
        but r1 always<0, r2>2*|f'|/|f''|
        r2=( |f'|+sqrt(|f'|^2+2*|f|*|f''|) )/|f''|

        test x(x-1) at 2+i
        f=1+3i f'=2x-1=3+2i f''=2
        r2=(|3+2*I|+sqrt(|3+2*I|^2+2*2*|1+3*I|))/2
        4.33502318885498454
        correct is 2.236
        */
        double fm=sqrt(newt.f_r.getMag_double());
        double fzm=sqrt(newtres.fz_r.getMag_double());
        double fzzm=sqrt(newt.fzz_r.getMag_double());
        newtres.first_neumaier2_re=(fzm + sqrt(fzm*fzm+2*fm*fzzm))/fzzm;
        newtres.first_neumaier2_im=0;

        /* Ostrowski theorem from page 1 of
        A theorem on clusters of roots of polynomial equations - A.M.Ostrowski
        https://epubs.siam.org/doi/epdf/10.1137/0707046
        R1C=2f/f' R1X=f'/f''  if R1C<R1X, contains 1 root in R1C and no roots in R1C..R1X
        R2C=max(sqrt(8f/f''), 8f'/f'') R2X=3/2*f''/f'''  if R2C<R2X, contains 2 roots in R2C and not roots in R2C..R2X
        */
        newtres.ostrowski_r1c=2*fm/fzm;
        newtres.ostrowski_r1x=fzm/fzzm;
        newtres.ostrowski_r2c=std::max(sqrt(8*fm/fzzm), 8*fzm/fzzm);
        //no f''' for r2x

        /* naive: approximate f with c(x-a)^m
        m=f'^2/(f'^2-f f'') = f'^2/f^2/(f'^2/f^2-f''/f)=G^2/H
        x-a=m/(f'/f)=m/G=G/H    looks good if |m_im|<|m_re|
        m*(x-a)=G^3/H^2

        trouble: singularities when f f''=f'^2 -> m=infinity, iteration jumps too far
                                    f'=0 -> m=0, m/(f/f') jumps too little
        */
        /*double g_re=currentWorker->toDouble(laguG.re_s);
        double g_im=currentWorker->toDouble(laguG.im_s);
        double g_mag=g_re*g_re+g_im*g_im;
        if (1e6*H_mag<=g_mag*g_mag)
        {
          newtres.first_naive_re=currentWorker->toDouble(r->re_s);
          newtres.first_naive_im=currentWorker->toDouble(r->im_s);
        }
        else
        {
          double g2_re=currentWorker->toDouble(laguG2.re_s);
          double g2_im=currentWorker->toDouble(laguG2.im_s);
          double h_re=currentWorker->toDouble(laguH.re_s);
          double h_im=currentWorker->toDouble(laguH.im_s);

          double m_re=(g2_re*h_re+g2_im*h_im)/H_mag;
          double m_im=(g2_im*h_re-g2_re*h_im)/H_mag;
          //couldn't find smooth function that:
          //1->1 2->2 3->3... 0->1 -1->1 i->1 -i->1
          //esp. since we need to have 1->1 exact and in neigborhood too
          if ((m_re<abs(m_im)*2))
          //if ((m_re<0.9) || (m_re<abs(m_im)*2)) //for m~0, we need something like sqrt(m): m is too small, 1 is too large
          {
            m_re=1;
            m_im=0;
          };
          newtres.first_naive_re=currentWorker->toDouble(r->re_s)-(m_re*g_re+m_im*g_im)/g_mag;
          newtres.first_naive_im=currentWorker->toDouble(r->im_s)-(m_im*g_re-m_re*g_im)/g_mag;
        }*/

        /* even naiver: show the 2 roots of c(x-a)(x-b) that have the same f, f', f''
        w.l.o.g. x=0
        c(x^2-(a+b)x+ab)=f''x^2/2+f'x+f   just solve Ax^2+Bx+C where A=f''/2 B=f' C=f
        //cx^2-c(a+b)x+cab=f''x^2/2+f'x+f
        -f'/f''+-sqrt(f'^2/f''^2-2*f/f'')
        f/(-f'+-sqrt(f'^2-2f''f))

        if x1 close to x2 (relative to x), use (x1+x2)/2 else use x1
        at |x1|=|x2|, 90 degrees..mult~2, use (x1+x2)/2
        at |x1|=|x2|, 60 degrees..mult~1, use x1
        at |x1|=0.8|x2|, 80% weight from x1
        at |x1|=0.5|x2|, 90% weight from x1
        at |x1|=0.3|x2|, use x1
        when x1~x2, correct guess is actually around 0.7 x1
        */
        a_re=newt.fzz_r.re.toDouble()/2;
        a_im=newt.fzz_r.im.toDouble()/2;
        MandelMath::complex_double_quadratic2(&newtres.first_naive1_re, &newtres.first_naive1_im,
                                              &newtres.first_naive2_re, &newtres.first_naive2_im,
                                              a_re, a_im, fz_re/2, fz_im/2, f_re, f_im);
        double n2_rmag=1/(newtres.first_naive2_re*newtres.first_naive2_re+newtres.first_naive2_im*newtres.first_naive2_im);
        //d=naive1/naive2
        double d_re=(newtres.first_naive1_re*newtres.first_naive2_re+newtres.first_naive1_im*newtres.first_naive2_im)*n2_rmag;
        double d_im=(newtres.first_naive1_im*newtres.first_naive2_re-newtres.first_naive1_re*newtres.first_naive2_im)*n2_rmag;
        double d_mag=(d_re*d_re+d_im*d_im);
        double w1=1, w2=0;
        if (d_re<-0.5) //angle>120deg, even if close in magnitude
        { w1=1; w2=0; newtres.naiveChoice=NewtonNaiveChoice::ncWide; }
        else if (d_mag<0.3*0.3)
        { w1=1; w2=0; newtres.naiveChoice=NewtonNaiveChoice::nc03; }
        else if (d_mag<0.5*0.5)
        { w1=0.9; w2=0.1; newtres.naiveChoice=NewtonNaiveChoice::nc05; }
        else if (d_mag<0.8*0.8)
        { w1=0.8; w2=0.2; newtres.naiveChoice=NewtonNaiveChoice::nc08; } //or just 1;0
        else if (d_re<-0.1)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          w1=1; w2=0.0;
          newtres.naiveChoice=NewtonNaiveChoice::nc100;
        }
        else if (d_re<0)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          //don't trust M here
          w1=1; w2=0.0;
          newtres.naiveChoice=NewtonNaiveChoice::nc90;
        }
        else if (d_re<0.1)
        {
          //most problematic case, some times best is [1,0] other times [1,1]
          //don't trust M here
          w1=1; w2=0.0;
          newtres.naiveChoice=NewtonNaiveChoice::nc80;
        }
        else if (d_re<0.5)
        {
          //can (try) use M here
          w1=1; w2=0;
          newtres.naiveChoice=NewtonNaiveChoice::nc60;
        }
        else
        {
          //can (try) use M here
          w1=newtres.firstMum_re-1; w2=0;
          newtres.naiveChoice=NewtonNaiveChoice::ncClose;
        }
        newtres.first_naive_re=w1*newtres.first_naive1_re+w2*newtres.first_naive2_re;
        newtres.first_naive_im=w1*newtres.first_naive1_im+w2*newtres.first_naive2_im;
        newtres.first_naive1_re=r_re+newtres.first_naive1_re;
        newtres.first_naive1_im=r_im+newtres.first_naive1_im;
        newtres.first_naive2_re=r_re+newtres.first_naive2_re;
        newtres.first_naive2_im=r_im+newtres.first_naive2_im;
        newtres.first_naive_re=r_re+newtres.first_naive_re;
        newtres.first_naive_im=r_im+newtres.first_naive_im;

        //Lagu=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
        //Muller solution of a quadratic equation =-c/(b+-sqrt(b^2-a*c))
        //Laguerre is the solution of ax^2+bx+c=0 with
        //  c=-n  b=f'/f  a=f'^2/f^2*(1-n/m+1/m)-f''/f*(1-n/m)=H*(1-n/m)+G^2/m
        //  G=f'/f   H=G^2-f''/f
        //  c=-1  b=G/n  a=(1/m-1/n)*H-G^2/n/m =G^2/n^2-(1/m-1/n)*(H-G^2/n)
        //>>c=-m  b=m*G/n  a=H*(1-m/n)-G^2/n
        /*
        a_re=currentWorker->toDouble(laguH.re_s)*(1-m*order1)-currentWorker->toDouble(laguG2.re_s)*order1;
        a_im=currentWorker->toDouble(laguH.im_s)*(1-m*order1)-currentWorker->toDouble(laguG2.im_s)*order1;
        double b_re=currentWorker->toDouble(laguG.re_s)*m*order1;
        double b_im=currentWorker->toDouble(laguG.im_s)*m*order1;
        MandelMath::complex_double_quadratic(
              &newtres.first_lagum_re, &newtres.first_lagum_im,
              a_re, a_im,
              b_re, b_im,
              -m, 0);
        newtres.first_lagum_re=currentWorker->toDouble(r->re_s)-newtres.first_lagum_re;
        newtres.first_lagum_im=currentWorker->toDouble(r->im_s)-newtres.first_lagum_im;
        */
        a_re=newt.laguH.re.toDouble()*(1-order1)-newt.laguG2.re.toDouble()*order1;
        a_im=newt.laguH.im.toDouble()*(1-order1)-newt.laguG2.im.toDouble()*order1;
        double bn_re=newt.laguG.re.toDouble();
        double bn_im=newt.laguG.im.toDouble();
        MandelMath::complex_double_quadratic2(
              &newtres.first_lagu1_re, &newtres.first_lagu1_im,
              &newtres.first_lagu1o_re, &newtres.first_lagu1o_im,
              a_re, a_im,
              bn_re*order1, bn_im*order1,
              -1, 0);
        //sometimes this way, sometimes that way... if (order1==0 && bn_re*newtres.first_lagu1_re-bn_im*newtres.first_lagu1_im>0)
        if (order1==0 && bn_re*newtres.first_lagu1_re-bn_im*newtres.first_lagu1_im<0)
        { //if b is exactly 0, complex_double_quadratic can't choose the right root, so we use original b before b*order1==0
          //need smaller root
          // | -c/(b+sqrt(b^2-a*c)) | > | -c/(b-sqrt(b^2-a*c)) |
          // | 1/(b+sqrt(b^2-a*c)) | > | 1/(b-sqrt(b^2-a*c)) |
          // current lagu1 is wrong if
          //   way 1: b . sqrt() <0   =  b . 1/lagu1 <0   = b . lagu1^T <0
          //   way 2: b . sqrt() <0   =  b . lagu1*a <0      a*lagu1=1/lagu1
          // also newton direction 1/G = 1/b   lagu1 wrong if 1/b . lagu1 <0  =  b^T . lagu1 <0  = b  .lagu1^T <0
          //now it doesn't //and now it does not // funny that it works if b  .lagu1^T >0
          double tre, tim;
          tre=r_re-newtres.first_lagu1_re;
          tim=r_im-newtres.first_lagu1_im;
          newtres.first_lagu1_re=r_re-newtres.first_lagu1o_re;
          newtres.first_lagu1_im=r_im-newtres.first_lagu1o_im;
          newtres.first_lagu1o_re=tre;
          newtres.first_lagu1o_im=tim;
        }
        else
        {
          newtres.first_lagu1_re=r_re-newtres.first_lagu1_re;
          newtres.first_lagu1_im=r_im-newtres.first_lagu1_im;
          newtres.first_lagu1o_re=r_re-newtres.first_lagu1o_re;
          newtres.first_lagu1o_im=r_im-newtres.first_lagu1o_im;
        }
      };
    prevm=lastm;
    //prevm_err=lastm_err;
    lastm=m;
    lastm_err=(mum_re-m)*(mum_re-m)+mum_im*mum_im;
    bool lagu_valid=false;
    bool newt_valid=false;
    if (order1>=0)
    {
      // a=1/(G/n +- sqrt( (1/m-1/n)*(H-G^2/n) ))
      // all but last few cycles can be done just in double precision
      //   but the cost of this compared to evaluation of f,f',f'' is negligible
      newt.laguX.assign(newt.laguG2);
      newt.laguX.lshift(-period); //G^2/n
      newt.laguX.rsub(newt.laguH); //H-G^2/n
      newt_tmp2.zero(lastm);
      newt_tmp2.recip();
      newt_tmp2.add_double(-order1); //1/m-1/n
      newt.laguX.mul(newt_tmp2); //(1/m-1/n)*(H-G^2/n)
      newt.laguX.sqrt();
      //normally we need to choose +- to maximize mag(G/n+-sqrt...) but that's numerically unstable
      if (newt.laguX.mulreT_tmp(newt.laguG)->isl0())
      {
        newt.laguX.chs();
      };
      newt.laguG.lshift(-period); //G/n
      newt.laguX.add(newt.laguG);
      //if 1/n~0: a=1/(0 +- sqrt( (1/m)*(H) )), m can still be 1..max
      //   fine if H!=0:       a=1/( sqrt( (1/m)*(H) )), m can still be 1..max
      //   if H==0: 1/G/(1/n + sqrt( (1/300-1/n)*(-1/n) ))=1/G* -i*sqrt(300*n)
      //   if H=G=0: 1/0
      //if G=0: a=1/(+- sqrt( (1/m-1/n)*(H) ))     m=1
      //   fine if H!=0: a=+-(sqrt(n/(n-1))*sqrt(f/-f''))       x^2+9 at 0: f=9 f''=2 -> +-3i
      //   if H=0: a=1/0
      //if H=0: a=1/G*m*(1 - i*sqrt(n/m-1))  m~n -> a=n/G;  m~300 -> a=-i/G*sqrt(n*300)
      //        a=1/G*m*n*(1/n - i*sqrt(1/m/n-1/n^2))
      double X_mag=newt.laguX.getMag_tmp()->toDouble();
      if (X_mag>=1e-60)
      {
        newt.laguX.recip_prepared();
        lagu_valid=true;
      };
      //else
      //we should move the guess a little and try again, but
      //  we can leave this to the caller
      //return 0;
    };
    if (gz_r_mag!=0)
    { //laguG is no more, was multiplied by order1
      //newton near multiroot:
      //f=(x-a)^m   f'=m*(x-a)^(m-1)  f/f'=(x-a)/m
      //Newton corrected for multiroot = f/f'*m
      //1/M=1-f''*f/f'^2   x-a=1/(f'/f-f''/f')
      newt.newtX.assign(newtres.fz_r);
      newt.newtX.recip();
      newt.newtX.mul(newt.f_r); //f/f'
      if (lastm>1)
      {
        /*newt_tmp2.zero(lastm);
        currentWorker->mul(newt.newtX.re, inst_newt_tmp2.getptr());
        currentWorker->mul(newt.newtX.im, inst_newt_tmp2.getptr());*/
        newt.newtX.mul_double(lastm);
      };
      newt_valid=true;
    };
    if (newtonCycle==0)
    {
      newtres.first_guess_newt.assign(*r);
      if (newt_valid)
      {
        newtres.first_guess_newt.sub(newt.newtX);
      };

      newtres.first_guess_lagu.assign(*r);
      if (lagu_valid)
      {
        newtres.first_guess_lagu.sub(newt.laguX);
      };
    };
    last_was_newton=true;
    if (!newt_valid)
    {
      if (!lagu_valid)
      {
        return 0;
      };
      if (!movedOff)
      {
        movedOff=true;
        newt.newtX.assign(newt.laguX);
        last_was_newton=false;
      }
      else
        return 0;
    }
    else if (!lagu_valid)
    {
    }
    else
    {
      if (period==1 || //need to help the luck a little...
          (fastHoming && (newtonCycle<2) && (lastm>1)))
      {
        newt.newtX.assign(newt.laguX);
        last_was_newton=false;
      }
      else
      {//take the smaller of newton and laguerre to a) avoid Newton's jumps when moving off the real axis, b) avoid Laguerre skipping to far root
        double N_mag=newt.newtX.getMag_double();
        double L_mag=newt.laguX.getMag_double();
        if (N_mag*1.05>L_mag) //5% will do no harm, and switch to Lagu can speed up convergence
        {
          newt.newtX.assign(newt.laguX);
          last_was_newton=false;
        };
      }
    }

    newt.prevR.assign(*r);
    r->sub(newt.newtX);
  } //for newtonCycle

  r->assign(newt.prevR);
  newtres.accy_tostop=std::max(1.0, prev_grmag/c.eps2());
  newtres.accy_multiplier=prev_accy_multiplier;
  newtres.accy_noise=prev_accy_noise;
  newtres.fz_r.assign(newt.prevGz);
  return lastm;
}

//result=0 means the period check failed; -1 means the check failed and the root returned is invalid
template <typename BASE>
int MandelEvaluator<BASE>::periodCheck(int period/*must =eval.lookper_lastGuess*/, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &root_seed, bool exactMatch)
{
  if (period<1)
  {
    dbgPoint();
    return -1;
  };
  MandelMath::number<BASE> &eval_fz_mag=eval.fz_mag;

  /*int aroundCount; //estimate multiplicity (mult-1)
  if ((currentData.store->lookper_prevGuess_>0) &&
      ((currentData.store->lookper_lastGuess % currentData.store->lookper_prevGuess_)==0))
    aroundCount=currentData.store->lookper_lastGuess / currentData.store->lookper_prevGuess_;
  else
    aroundCount=0;
  if (aroundCount==0)
    aroundCount=1; //a fresh nearestIteration means this is a new atom, so mult=2
  */
  //look for root nearest to C - better stability of newton/laguerre
  //MandelMath::complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
  mandelData.root.assign(root_seed);
  //root.sqr();
  //root.add(c);

  /*checked before call to periodCheck if (currentWorker->toDouble(&eval.lookper_totalFzmag)>=MAGIC_MIN_SHRINK) //correct totalFZmag?
  { //TODO: is this correct? we're not evaluating at root, just some point around here...
    return -1;
  };*/

  if (exactMatch)
  {
    if (!loope.eval_zz(period, c, mandelData.root, true))
      return -1;
    newtres.fz_r.assign(loope.f_z);
    newtres.accy_tostop=1;//try 1
    newtres.accy_multiplier=std::max(1.0, 1/loope.f_z.getMag_double());//try 5
    double f_z_mag=MandelMath::sqr_double(1+loope.f_z.re.toDouble())+
                   MandelMath::sqr_double(loope.f_z.im.toDouble());
    newtres.accy_noise=period*f_z_mag/
        loope.f_z_minmag.toDouble(); //the error at minimum f_z is magnified to f_z_mag
  }
  else
  {
    if (period>MAX_PERIOD)
    {
      return -1; //special case
    };
    int newtRes=laguerre(period, c, &mandelData.root, true);
    totalNewtonIterations+=newtres.cyclesNeeded*period;
    if (newtRes<=0)
    { //this, of course, means that Newton() should be improved, not that there's a problem with the numbers
      return -1; //e.g. evaluating the initial guess mand.root leads to overflow immediately
    };
  }

  this->interior.fz.assign(newtres.fz_r);
  interior.fz.re.add_double(1);
  double fz_log=std::log(interior.fz.getMag_double())/2;
  double ori_over1=interior.fz.getMag1_tmp()->toDouble(); //could use (1+r)^2+i^2-1=2*r+r^2+i^2
  //double fz_mag=newtres.fz_r.getMag_double();
  if (ori_over1>0.01)//often we are at long root long period when correct is short/short (inner edge of a bulb) 0.002) //fz_mag>0 && fz_mag*fz_mag>25*currentWorker->eps2())
    return -1;
  else if (period==1)
  { //cannot shorten or extend
    if (ori_over1>0)
      return -1;
    else
      return period;
  }
  if (ori_over1<-0.001)
  {
    //we can test one way: estimate accuracy of root => estimate distance after p iterations
    //x=dist_around
    //f^p(r)=r+eps && f^p(r+x)=f^p(r)+f^p'*x => x=eps/f^p' so the principle should be good
    double dist_aroundN=4.192*16*(newtres.accy_multiplier*newtres.accy_tostop)*c.eps2();
      //maybe *period or something instead of 16?
      //1.097*16 for period=35 first=1
      //1.137*16 for period=43 first=1
      //1.1583*16 for period=4752 first=432 at doubledouble
      //1.1591*16 for period=9 first=1
      //1.2389*16 for period=98 first=1
      //1.2930*16 for period=34 first=1
      //1.6644*16 for period=46 first=1
      //2.0500*16 for period=7680 first=1536
      //3.0196*16 for period=68 first=1
      //3.2807*16 for period=512 first=64
      //4.1912*16 for period=832 first=64

    //or test another way: estimate distance to other roots
    //try neumaier bound for 2 roots
    //r2=( |f'|+sqrt(|f'|^2+2*|f|*|f''|) )/|f''|
    //for |f|=0: r2=2*|f'|/|f''|
    double dist_around2;
    {
      double fz=newtres.fz_r.getMag_double();
      double fzz=loope.f_zz.getMag_double();
      //for outer roots, seems to approximate nicely the root's pool (good)
      //for central root, it covers all the satellites (don't care if the central is repelling)

      //if we have f' and f'' after p iterations,
      //after k*p iterations we have f_z=(f'+1)^k-1, f_zz=f''*(f'+1)^(k-1)*((f'+1)^k-1)/(f')
      //f_z/f_zz=((f'+1)^k-1)/( f''*(f'+1)^(k-1)*((f'+1)^k-1)/(f') )
      //f_z/f_zz=f'/f'' /(f'+1)^(k-1) = f'/f''/(f'+1)^k*(f'+1)
      // //approx x^n where x in R, 0<=x<=1, n=period/i
      // //exp(n*ln(x)) ~ exp(n*(x-1))  but exp won't work -> just use std::pow
      //k=p/i  f_z/f_zz=f'/f''*(f'+1) /(f'+1)^(p/i) but (f'+1)^(p/i) is eval.fz_mag1+1
      //dist_around2=4*fz/fzz*(1+ori_over1);
      //it only works IF we are at the root, or the bound gets bigger and includes the false root we are examining
      //try to include some safety factor <- remove 4*
      dist_around2=fz/fzz*interior.fz.getMag_double();//1+ori_over1);
    }
    newt.f_r.assign(mandelData.root);
    newtres.fz_r.zero(1, 0);
    loope.f_zz.zero(0, 0);
    eval_fz_mag.zero(3); //ignore first pass

    int firstBelow1dividing=-1;
    for (int i=0; i<period; i++)
    { //fz_mag1 has to be (ori_over1+1)^(i/period)-1
      //upper bound: ori_over1*(i/period)
      //lower bound: x:=i/period;  x*(x*(ori_over1-ln(ori_over1+1))+ln(ori_over1+1))
      //                                 parabola that starts like (ori_over1+1)^x-1 and f(1)=ori_over1
      //               or ori_over1 if ori_over1<-0.7968
      if ((//currentWorker->isle0(eval.fz_mag1.ptr) ||
           //(MandelMath::sqr_double(eval.fz_mag1.toDouble())<=fz_e2_m)))
           (eval_fz_mag.toDouble()<=ori_over1*MandelMath::sqr_double(i/(double)period)+1+1e-8))) //should be <=(1+fzmag)^(i/period)-1
      {
        if (firstBelow1dividing<0)
          if ((period % i)==0)
          {
            //needs more checks than that, e.g. fz_mag^(period/i) <=~ final fz_mag
            //per-actual  per-found  root-found   |   status at short
            //  short       short       short         no long to loop over
            //  short       short        long         does not solve newton
            //  short        long       short         |f-r|<eps
            //  short        long        long         |fz|>1
            //  long        short       short         |fz|>1
            //  long        short        long         does not solve newton
            //  long         long       short         |fz|>1
            //  long         long        long         no short to test
            //if (currentWorker->isle(f_r.getMagTmp(), root.getMagTmp()))
            double dist=newt.f_r.dist2_double(mandelData.root);
            //double dist_around2_fixed=dist_around2*std::pow(1+eval.fz_mag1.toDouble(), period/i-1);
              //maybe find lower bound on that using 2^-n
            //or do it easy
            //double tmp=1+eval.fz_mag1.toDouble(); //sometimes fz_mag1 is too small, making dist_fixed too large
            double tmp=std::exp(i*fz_log/period);//1+ori_over1*MandelMath::sqr_double(i/(double)period);
            double dist_around2_fixed;
            if (tmp==0)
              dist_around2_fixed=dist_around2/tmp;
            else
              dist_around2_fixed=dist_around2/tmp; //maybe round down to 2^-n where n=(period/i-1)*fz_mag1
            double dist_around_direct=0.25* 4*newtres.fz_r.getMag_double()/loope.f_zz.getMag_double();
            (void)dist_around_direct;
            if (dist_around2_fixed==0)
              nop();
            else if (dist_aroundN>dist_around2_fixed) //happens when around2==0
              nop();
            if (dist<=dist_aroundN && dist<=dist_around2_fixed) //dist=0, limit=0 -> match -> use "<="
            {
              firstBelow1dividing=i; //short long short
              break;
            }
            else if (dist>dist_aroundN && dist>dist_around2_fixed)
            { }
            else if (dist<25*dist_aroundN)
            {
              nop(); //plenty
              firstBelow1dividing=i;
            }
            else if (dist>1e-7)
            {
              nop();
              firstBelow1dividing=i;
            }
            else
            {
              nop(); //at bigger zoom/exact match
              firstBelow1dividing=i;
            }
          };
      };
      //f_zz:=2*(f_z*f_z + f*f_zz)
      loope.f_zz.mul(newt.f_r);
      bulb.t1.assign(newtres.fz_r);
      bulb.t1.sqr();
      loope.f_zz.add(bulb.t1);
      loope.f_zz.lshift(1);
      //f_mag could hardly be >4 since it's tested in Newton (as well as fz_mag, BTW)
      newtres.fz_r.mul(newt.f_r);
      newtres.fz_r.lshift(1);
      eval_fz_mag.assign(*newtres.fz_r.getMag_tmp());
      if (eval_fz_mag.toDouble()>LARGE_FLOAT2)
      {
        dbgPoint();
        return -1; //so is it checked or not
      }
      //f:=f^2+c
      //use EXACTLY same code as the main loop: if temporaries are more precise at one place, it will never work!
      newt.f_r.sqr();
      //we don't need last f, otherwise should do if (i+1==period) { add(c-r); }
      newt.f_r.add(c);
    };
    /*bool fz_r_mag_over1=(newtres.fz_r.mag_cmp_1()>0);
    if (fz_r_mag_over1)
    {
      dbgPoint(); //mag_cmp_1 should be ==ori_over1 so never over 1
      return -1;//inevitably result:=-1
    }*/

    if (firstBelow1dividing<1)
      return period;
    else
      return firstBelow1dividing;
  }
  else
  {
    bool tested_case;
    double c_re=c.re.toDouble();
    double c_im=abs(c.im.toDouble());
    if ((c_re==0.25 && c_im==0.5) ||
        (c_re==-1.75 && c_im==0) ||
        (c_re==-1.25 && c_im==0) ||
        (c_re==-0.75 && c_im==0) ||
        (c_re==-1 && c_im==0.25))
      //0.237060546875 0.531494140625
      tested_case=true; //tested to work
    else
      tested_case=false;
    //first find short root and see if it is attractive
    if (!loope.eval_multi(period, c, mandelData.root, interior.fz, 0)) //TODO: tolerance ~ accy_noise
    {
      dbgPoint(); //should be detected earlier
      return -1;
    };
    if (loope.multi>1)
    {
      int short_period=period/loope.multi;
      int short_newtRes=laguerre(short_period, c, &loope.sumA, true);
      totalNewtonIterations+=newtres.cyclesNeeded*short_period;
      if (short_newtRes<=0)
      {
        return -1;
      };
      this->interior.fz.assign(newtres.fz_r);
      interior.fz.re.add_double(1);
      double short_over1=interior.fz.getMag1_tmp()->toDouble();
      //bool fz_r_mag_over1=;//bulb.bulbe.f_z.mag_cmp_1()>0);
      if (short_over1>0 && short_over1*short_over1>25*c.eps2())
      { //central surely repulsive
        nop();//return -1;
      }
      else if (short_over1<0 && short_over1*short_over1>25*c.eps2())
      { //central surely attractive
        mandelData.root.assign(loope.sumA);
        return short_period;
      }
      else
      {
        if (!tested_case)
          nop();
        if (short_over1>0) //bulb.bulbe.f_z.mag_cmp_1()>0);
          nop();//return -1;
        else
        {
          mandelData.root.assign(loope.sumA);
          return short_period;
        }
      }
    }
    //then check if we are at long period with long root and is attractive
    if (ori_over1>0 && ori_over1*ori_over1>25*c.eps2())
    { //long period with our root is repulsive
      nop();//return -1;
    }
    else if (ori_over1<0 && ori_over1*ori_over1>25*c.eps2())
      return period; //our root is attractive and we did check central root above so it's an outer root
    else
    { //at cardioid cusp, there's no multi to reduce but f_z=0
      nop();
      if (ori_over1>0) //bulb.bulbe.f_z.mag_cmp_1()>0);
        nop();//return -1;
      else
        return period;
    }
    //last we can be at long period, short root which is repulsive, find long root and check
    //TODO: findBulbBase...
    /* old notes:
       deltac = c - base_c
       central r = base_r + deltac * f_cc / (-2 f_zc)  because dc^2/2 f_cc + dc dz f_zc=0
       outer r = (deltac/(xc-base_c))^(1/m) * (xc-base_r) + base_r
         for c=xc we get: outer r = ((xc-base_c)/(xc-base_c))^(1/m) * (xc-base_r) + base_r = xc
         for c=base_c, we get: outer r = (0/(xc-base_c))^(1/m) * (xc-base_r) + base_r = base_r

       simple guess: approximate f with k(x-r)((x-r)^m-deltac)
         outer r = deltac^(1/m) + central_r
         but why deltac? seems to work...
       better guess: approximate f with a(x-r)((x-r)^m-b*deltac)
         such that xc is one of the outer roots
         (xc-rb)^m-b*(xc-cb)=0
           xc-cb=(+-) 1/(fzc-fc/fz*fzz)
           xc=c+(fz+1)/(fzc-fc/fz*fzz)
           deltac=fz/(fzc-fc/fz*fzz)
           (xc-rb)^m/(xc-cb)=b
           (x-r)^m-b*(fz+1)/(fzc-fc/fz*fzz)=0
           (x-r)^m=(xc-rb)^m*(fz+1)
        >> (x-r)=(xc-rb)*(fz+1)^(1/m)
             deltac/(xc-base_c) = fz/(fzc-fc/fz*fzz)/1*(fzc-fc/fz*fzz) = fz
             so with some renaming of fz+1 to fz we get the old result
    */
    if (loope.multi>1)
    { //sumA went through short newton already
      if (!loope.eval2(period, c, loope.sumA, false))
        return -1;
      //        0=f=0+C*fc+R*fz -> R=-C*fc/fz
      //target_fz-1=fz+C*fzc+R*fzz   fz/(fc/fz*fzz-fzc)=C=-deltaC    (fz-target_fz+1)/(fzc-fc/fz*fzz)=deltaC

      //0-1=fz+C*fzc+R*fzz   (fz-0+1)/(fzc-fc/fz*fzz)=deltaXC   xc-basec=(xc-C)-(basec-C)=deltac-deltaXC=-1/(fzc-fc/fz*fzz)
      //baser=centralr   xc=c+deltaXC=c+(fz+1)/(fzc-fc/fz*fzz)
      //outer r = (-fz)^(1/m) * ((fz+1)/(fzc-fc/fz*fzz)+c-base_r) + base_r
      bulb.t3.assign(loope.f_z);
      bulb.deltac.assign(loope.f_z);
      loope.f_z.re.add_double(-1);
      bulb.deltac.recip();
      bulb.deltac.mul(loope.f_c); // fc/fz
      //bulb.deltar.assign(&deltac);
      bulb.deltac.mul(loope.f_zz);
      bulb.deltac.rsub(loope.f_zc); //fzz fc/fz - fzc
      {
        double mag=bulb.deltac.getMag_double();
        if (mag<1e-8)
        {
          nop();
          return -1;
        }
        else if (mag<0.1)
          nop();
      }
      bulb.deltac.recip(); // 1/(fzc-fc/fz*fzz)
      bulb.t2.assign(bulb.deltac); // -"-
      bulb.deltac.mul(loope.f_z);  //deltac, cb=c-deltac
      bulb.deltac.chs();
      bulb.deltac.root_approx(loope.multi);
      bulb.deltac.chs();

      //s3=fz+1
      bulb.t2.mul(bulb.t3); //(fz+1)/(fzc-fc/fz*fzz)    c-s2/multi~xc?
      bulb.t1.zero(loope.multi, 0);
      bulb.t1.recip();
      bulb.t2.mul(bulb.t1);
      bulb.t2.add(c);
      bulb.t2.sub(loope.sumA); //xc-centralr
      bulb.t3.assign(loope.f_z);
      bulb.t3.root_approx(loope.multi);
      bulb.t3.mul(bulb.t2);
      bulb.t1.zero(std::sqrt(loope.multi), 0);
      bulb.t3.mul(bulb.t1); //((fz+1)/(fzc-fc/fz*fzz) /m+c-rb) *(fz)^(1/m) *sqrt(m)
      //     should have been  ((fz+1)/(fzc-fc/fz*fzz)   +c-rb) *(fz)^(1/m)

      bulb.deltac.add(loope.sumA); //TODO: it's actually around base_r, not central_r
      bulb.t3.add(loope.sumA);
      laguerre(period, c, &bulb.t3, true);
      newtres.fz_r.re.add_double(1);
      double outer_over1=newtres.fz_r.getMag1_tmp()->toDouble();
      bool big_enough1=outer_over1*outer_over1>25*c.eps2();
      bool big_enough2=outer_over1*outer_over1>newtres.accy_multiplier*c.eps2();
      if (big_enough1!=big_enough2)
        nop();
      if (outer_over1>0 && big_enough1 && big_enough2)
      { //outer surely repulsive
        nop();//return -1;
      }
      else if (outer_over1<0 && big_enough1 && big_enough2)
      { //outer surely attractive
        mandelData.root.assign(bulb.t3);
        return period;
      }
      else
      {
#if 0 //worse than just comparing outer_over1>0
        //last desperate try: Re(f''^m)<0
        bulb.t1.assign(&newt.fzz_r);
        bulb.t1.pow_int(loope.multi, &tmp);
        if (!bulb.t1.re.isl0())
          nop();//return -1;
        else
        {
          mandelData.root.assign(&bulb.t3);
          return period;
        }
#else
        if (outer_over1>0) //bulb.bulbe.f_z.mag_cmp_1()>0);
          nop();//return -1;
        else
        {
          mandelData.root.assign(bulb.t3);
          return period;
        }
#endif
      }
    }

    return -1;
  }







#if 0
  //complex f_r(currentWorker, &newt.f_r_re, &newt.f_r_im, true);
  //complex fz_r(currentWorker, &newtres.fz_r_re_, &newtres.fz_r_im_, true);
  //double dist_around=5*currentWorker->eps2()/currentWorker->toDouble(fz_r.getMagTmp());
  //need at least safety factor of 14.2
  //correct root maps to .to_stop, then also 2*error can map to .to_stop
  //and we oscillate +-error so that's 4*to_stop, or 16 in dist squared
  double dist_around=16*(newtres.accy_multiplier*newtres.accy_tostop)*currentWorker->eps2();
  newt.f_r.assign(&currentData.root);
  double f_e_re=sqrt(dist_around/6), f_e_im=0;
  double f_e2_re=sqrt(dist_around), f_e2_im=0;
  newtres.fz_r.zero(1, 0);
  double fz_e_m=0, fz_e2_re=0, fz_e2_im=0, fz_e2_m=0;
  eval.fz_mag1.zero(1);
  aroundCount=0;
  bool someBelow1=false;

  //TODO: try to guess the radius of "around" roots   double dist_around=newt.fz*2/fzz  f(r+x)=r+x=r+f'x+f''x^2/2

  //int firstBelow1=-1;
  int firstBelow1dividing=-1;
  for (int i=0; i<period; i++)
  {
    //if (fz_mag1 && currentWorker->isl0(fz_mag1)) //I think we intentionally skip last fz_mag
    if (//replaced with fz_mag1.zero(1) fz_mag1.asf64 &&  //I think we intentionally skip last fz_mag
        (currentWorker->isle0(eval.fz_mag1.ptr) ||
         //(MandelMath::sqr_double(eval.fz_mag1.toDouble())<=fz_e2_m)))
         (eval.fz_mag1.toDouble()<=fz_e2_m)))
    {
      someBelow1=true;
      //if (firstBelow1<0)
        //firstBelow1=i;
      if (firstBelow1dividing<0)
        if ((period % i)==0)
        {
          //needs more checks than that, e.g. fz_mag^(period/i) <=~ final fz_mag
          //per-actual  per-found  root-found   |   status at short
          //  short       short       short         no long to loop over
          //  short       short        long         does not solve newton
          //  short        long       short         |f-r|<eps
          //  short        long        long         |fz|>1
          //  long        short       short         |fz|>1
          //  long        short        long         does not solve newton
          //  long         long       short         |fz|>1
          //  long         long        long         no short to test
          //if (currentWorker->isle(f_r.getMagTmp(), root.getMagTmp()))
          double dist=newt.f_r.dist2_double(&currentData.root);
          if (dist<dist_around)//3.4e-28) //related to newton's accyBound=3e-28/period^2
            firstBelow1dividing=i; //short long short
          else if (dist>1e-7)
          { }
          else
          {
            nop();
          }
          //does "else" even exist? should be always true? 0.00804 vs 1e-31 nope
        };
    };
    //fm could hardly be >4 since it's tested in Newton (as well as f_z_r, BTW)
    double fre=currentWorker->toDouble(newt.f_r.re);
    double fim=currentWorker->toDouble(newt.f_r.im);
    //fz:=2*f*fz
    {
      double fzre=currentWorker->toDouble(newtres.fz_r.re);
      double fzim=currentWorker->toDouble(newtres.fz_r.im);
      double re=2*(fzre*f_e_re-fzim*f_e_im); //fz_e=2*fz_r*f_e
      double im=2*(fzim*f_e_re+fzre*f_e_im);
      fz_e_m=re*re+im*im;
      //ez:=2*(fz*e+ez*f +ez*e)
      re=2*(fzre*f_e2_re-fzim*f_e2_im+fz_e2_re*fre-fz_e2_im*fim +fz_e2_re*f_e2_re-fz_e2_im*f_e2_im);
      im=2*(fzre*f_e2_im+fzim*f_e2_re+fz_e2_re*fim+fz_e2_im*fre +fz_e2_re*f_e2_im+fz_e2_im*f_e2_re);
      fz_e2_re=re;
      fz_e2_im=im;
      fz_e2_m=(fz_e2_re*fz_e2_re+fz_e2_im*fz_e2_im);
    }
    newtres.fz_r.mul(&newt.f_r);
    newtres.fz_r.lshift(1);
    eval.fz_mag1.assign(newtres.fz_r.getMag1_tmp());
    if (currentWorker->toDouble(eval.fz_mag1.ptr)>LARGE_FLOAT2)
      return -1; //so is it checked or not
    {
      double re=fre*f_e_re-fim*f_e_im; //f_e=f_r*f_e
      f_e_im=fim*f_e_re+fre*f_e_im;
      f_e_re=re;
      //e:=2*f*e +e^2
      re=2*(fre*f_e2_re-fim*f_e2_im +f_e2_re*f_e2_re-f_e2_im*f_e2_im);
      f_e2_im=2*(fre*f_e2_im+fim*f_e2_re +2*f_e2_re*f_e2_im);
      f_e2_re=re;
    }
    //f:=f^2+c
    //use EXACTLY same code as the main loop: if temporaries are more precise at one place, it will never work!
    newt.f_r.sqr();
    newt.f_r.add(c);
  };
  //normally we would check |ff|<1 but it needs to be precise
  //first reduce fz_re, fz_im to 2nd octant:
  //  fz.re>=0, fz.im>=0, fz.im>=fz.re
  bool fz_r_mag_over1=(newtres.fz_r.mag_cmp_1()>0);
  if (fz_r_mag_over1)
  {
    /*don't need right now if (newt.f_r.isequal(&currentData.root))
      return -2;
    else*/
      return -1;//inevitably result:=-1
  }

  //evaluate F_c^period at r ; its abs must be below 1 for the point to attract

  if (!someBelow1)
    return period; //seems to work
  //if (firstBelow1!=firstBelow1dividing)
    //dbgPoint();
  if (firstBelow1dividing<1)
    return period;
  else
    return firstBelow1dividing;
#endif
}

template <typename BASE>
int MandelEvaluator<BASE>::estimateInterior(int period, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &root)//, InteriorInfo *interior)
{
  /*MandelMath::complex &f=newt.f_r;
  MandelMath::complex &fz=interior.fz;
  MandelMath::complex &fc=newt.laguG;
  MandelMath::complex &fzz=newt.fzz_r;
  MandelMath::complex &fzc=newt.laguG2;
  MandelMath::complex &fcc=newt.laguH;
  // Initial values:  f = r;  fc = 0;  fz = 1;  fzz = 0;  fzc = 0;
  f.assign(root); //z   z^2+c  (z^2+c)^2+c
  interior.fz.zero(1, 0);           //1   2z     4z^3+4cz=2*2z(z^2+c)
  fc.zero(0, 0);           //0   1      2z^2+2c+1=2*1*(z^2+c)+1
  newt.fzz_r.zero(0, 0);          //0   2      12z^2+4c=2(4z^2+2z^2+2c)=2*((2z)^2+(z^2+c)*2)
  fzc.zero(0, 0);          //0   0      4z=2*(2z*1+(z^2+c)0)
  fcc.zero(0, 0);*/
  //interior.fz_mag.zero(1);
  //MandelMath::number_instance_fix<WORKER_MULTI, Newt::iiw_tmp2> inst_newt_tmp2(currentWorker, &newt.tmp2);
  //MandelMath::number_instance_fix<WORKER_MULTI, InteriorInfo::iiw_inte_abs> inst_interior_inte_abs(currentWorker, &interior.inte_abs);
  if (!loope.eval2_mag(period, c, root))
  {
    //nobody checked f_zc until now so it can overflow dbgPoint(); //does it? //yes sometimes it does...
    interior.inte_abs.zero(-1);
    return -1;
  };
  interior.fz.assign(loope.f_z); //for later read

  //imma gonna skippa another test here
  //  of derivatives<1 -> would refine period

  newt.tmp2.assign(loope.f_z_mag);
  //newt.tmp2.sqrt(); //(1-|fz|)/(...) same at centers, touches near edges; official formula overlaps 2x near edges
  newt.tmp2.add_double(-1); //|fz|^2-1
  /*if (abs(newt.tmp2.toDouble())<6e-18) //for c=0.25+0.5i we find r=0.0002615+0.5000231 -> f_z_mag=1.00009284... (per=4/1)
  { //parabolic point
    dbgPoint();
    interior.inte_abs.zero(0);
    interior.inte.zero(0, 0);
    return period;
  };*/
  //                    1-|fz|^2          .   (1-fz)(1-fz fz^T)/(fzc*(1-fz) + fzz fc)
  // interior=  -----------------------   .
  //            | fzc + fzz fc/(1-fz) |   .
  //dropping d signs: (1-1/z^2)/(1/z/c+1/z/z/c*z)=(1-1/z^2)/(1/z/c+1/z/c)~(c/z) or units of (d/dz)/(d/dc)
#if 0
  //(|fz|^2-1)/(fzz fc/(fz-1) - fzc)
  currentWorker->assign(interior.inte.re, newt.tmp2.ptr);
  currentWorker->zero(interior.inte.im, 0); //|fz|^2-1
  newt.tmp1.assign(&bulb.bulbe.f_z);
  currentWorker->add_double(newt.tmp1.re, -1); //fz-1
  //skip this step because it's just not right: interior:=newt.tmp2 * tmp1/|tmp1|
  newt.tmp1.recip();
  newt.tmp1.mul(&bulb.bulbe.f_c);
  newt.tmp1.mul(&bulb.bulbe.f_zz);
  newt.tmp1.sub(&bulb.bulbe.f_zc); //fzz fc/(fz-1) - fzc
  if (newt.tmp1.is0())
  { //probably wrong period, should not happen
    interior.inte_abs.zero(5);
    interior.inte.zero(5, 0);
    return period;
  };
  newt.tmp1.recip();
  interior.inte.mul(&newt.tmp1);
#else
  //(fz-1)(|fz|^2-1)/(fzz fc - (fz-1)*fzc)
  loope.f_z.re.add_double(-1); //fz-1
  newt.tmp1.assign(loope.f_z);
  newt.tmp1.mul(loope.f_zc);
  interior.inte.assign(loope.f_zz);
  interior.inte.mul(loope.f_c);
  interior.inte.sub(newt.tmp1); //inte = fzz*fc - (fz-1)*fzc
  if (interior.inte.is0())
  { //parabolic point
    dbgPoint();
    interior.inte_abs.zero(0);
    interior.inte.zero(0, 0);
    return period;
  };
  interior.inte.recip();
  interior.inte.mul(loope.f_z); //(fz-1)/(...)
  interior.inte.mul(newt.tmp2); //(fz-1)(|fz|^2-1)/(...)
#endif
  interior.inte_abs.assign(*interior.inte.getMag_tmp());
  interior.inte_abs.sqrt();
  if (!newt.tmp2.isle0())
  {
    interior.inte_abs.chs();
  }
  return period;
}

template <typename BASE>
void MandelEvaluator<BASE>::mandel_until_bailout()
{
  for (int i=0; i<100; i++) //should be enough to reach 10000^2 except around (-2, 0)
  {
    double f_mag=mandelData.f.getMag_double();
    if (f_mag>1e8)
      return;
    //fc_c:=2*f*fc_c+1
    mandelData.fc_c.mul(mandelData.f);
    mandelData.fc_c.lshift(1);
    mandelData.fc_c.re.add_double(1);
    double fc_c_mag=mandelData.fc_c.getMag_double();
    if (fc_c_mag>LARGE_FLOAT2)
    {
      mandelData.store->rstate=MandelPointStore::ResultState::stBoundary;
      return;
    };
    //fz_c_mag:=4*fz_c_mag*f.mag
    mandelData.fz_c_mag.mul(*mandelData.f.getMag_tmp());
    mandelData.fz_c_mag.lshift(2);
    double fz_c_mag=mandelData.fz_c_mag.toDouble();
    if (fz_c_mag>LARGE_FLOAT2)
    {
      mandelData.store->rstate=MandelPointStore::ResultState::stDiverge;
      return;
    };
    //f:=f^2+c
    mandelData.f.sqr();
    mandelData.f.add(currentParams.mandel.c);
    mandelData.store->iter++;
  };
}

template <typename BASE>
void MandelEvaluator<BASE>::julia_until_bailout()
{
  for (int i=0; i<100; i++) //should be enough to reach 10000^2 except around (-2, 0)
  {
    double f_mag=juliaData.f.getMag_double();
    if (f_mag>1e8)
      return;
    //fz_z:=2*f*fz_z
    juliaData.fz_z.mul(juliaData.f);
    juliaData.fz_z.lshift(1);
    //juliaData.fz_z.re.add_double(1);
    double fz_z_mag=juliaData.fz_z.getMag_double();
    if (fz_z_mag>LARGE_FLOAT2)
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stBoundary;
      return;
    };
    //fz_c_mag:=4*fz_c_mag*f.mag
    juliaData.fz_z_mag.mul(*juliaData.f.getMag_tmp());
    juliaData.fz_z_mag.lshift(2);
    double fz_c_mag=juliaData.fz_z_mag.toDouble();
    if (fz_c_mag>LARGE_FLOAT2)
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stDiverge;
      return;
    };
    juliaData.since0fzm.mul(*juliaData.f.getMag_tmp());
    juliaData.since0fzm.lshift(2);
    //f:=f^2+c
    juliaData.f.sqr();
    juliaData.f.add(currentParams.julia.c);
    juliaData.nearmr.tap(currentParams.julia, &newt.newtX);
    juliaData.store->iter++;
    if (juliaData.store->near0iter_1==juliaData.store->iter)
    {
      juliaData.since0fzm.zero(1);
    };
  };
}

template <typename BASE>
void MandelEvaluator<BASE>::evaluateMandel()
{
  /*{
    eval.near0fmag.assign(*currentData.near0f.getMag_tmp(&tmp)); //TODO: update on changing near0f
  }*/

  for (; (mandelData.store->iter<currentParams.maxiter) &&
         (mandelData.store->rstate==MandelPointStore::ResultState::stUnknown) &&
         this->workIfEpoch==this->busyEpoch;
       mandelData.store->iter++)
  { //? near0iter gives cycle offset, near_first_z_iter gives period; for mandel near0iter+1=near_first_z_iter //not really ;also we know julia period in advance
    // peri   1   2   3  jul   jul
    // near0  1   2   3   1     2
    // nearz  1   2   3   3     3
    // test   0       0   0     0
    // test   3       9  7,10  8,11
    // test   6      18 16,19 17,20
    // test  12      36 34,37 35,38
    //no worky, can't really find period from nearz and near0
    int rem;
    int quo;
    if (mandelData.store->iter==0)
    {
      rem=0;
      quo=1;
    }
    else
    {
      int per=mandelData.store->nearziter_0;
      rem=(mandelData.store->iter-mandelData.store->near0iter_1)%(3*per);
      quo=(mandelData.store->iter-mandelData.store->near0iter_1)/(3*per);
    }
    if (rem==0)
    {
      if ((quo&(quo-1))==0) //also at iter==0  //TODO: maybe better 3*(2^k-1)*near not 3*(2^k)*near
      { // //need k*iter for f' to start at the worst moment to reduce false positives; need k*iter-1 for good near0 -> switch to nearc
        mandelData.store->lookper_startiter=mandelData.store->iter;
        mandelData.lookper_startf.assign(mandelData.f);
        mandelData.lookper_nearr.assign(mandelData.f);
        if (mandelData.store->iter<=1)
          mandelData.lookper_nearr_dist.assign(*mandelData.f.getMag_tmp());
        else
          mandelData.store->lookper_nearr_dist_touched=false;//currentWorker->assign(&currentData.lookper_nearr_dist, f.dist2_tmp(&c));
        mandelData.store->lookper_lastGuess=0;
        mandelData.lookper_totalFzmag.zero(1.0);
      };
    }
    const MandelMath::number<BASE> *f_mag1=mandelData.f.getMag_tmp();
    if (!f_mag1->isle(currentParams.mandel.bailout)) //exactly 4..not yet
    {
      mandelData.store->rstate=MandelPointStore::ResultState::stOutside;
      //theory says the relative error in estimate is less than 3/bailout for large bailout
      //so lets move out a bit
      mandel_until_bailout();//, &currentData.f, &currentData.fc_c); //may switch state to stBoundary
      if (mandelData.store->rstate!=MandelPointStore::ResultState::stOutside)
      {
        //currentWorker->zero(&currentData.exterior_avoids, 0);
        //currentWorker->zero(&currentData.exterior_hits, 0);
        mandelData.store->exterior_avoids=0;
        mandelData.store->exterior_hits=0;
      }
      else
      {
        // !! but that's for Julia not Mandelbrot - will try fc[c] instead of fc[z]
        //https://www.evl.uic.edu/hypercomplex/html/book/book.pdf p17, p29
        //G=ln(sqrt(f_mag))/2^iter   G'=sqrt(fc_c_mag)/(2^iter*sqrt(f_mag))       sinh(G)=(exp(G)-exp(-G))/2
        //sinh(G)/(2*exp(G)*G') < exterior < 2*sinh(G)/G'
        //(1-exp(-2G))/(4*G') < exterior < (exp(G)-exp(-G))/G'
        //for high iter (small G) we can use exp(x)=1+x, exp(-x)=1-x, sinh(G)=(1+G-(1-G))/2=G
        //  (2*G)/(4*G') < exterior < (2*G)/G'
        //  ln(f_mag)*sqrt(f_mag/fc_c_mag)/4 < exterior < ln(f_mag)*sqrt(f_mag/fc_c_mag)
        //  E/4 < exterior < E
        //otherwise, 1/G'=sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-exp(-2*ln(sqrt(f_mag))/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(ln(sqrt(f_mag))/2^iter)-exp(-ln(sqrt(f_mag))/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-1/exp(ln(f_mag)/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(ln(f_mag)/2/2^iter)-1/exp(ln(f_mag)/2/2^iter))*sqrt(f_mag/fc_c_mag)*2^iter
        //  define X=ln(f_mag)/2/2^iter
        //  (1-1/exp(X)^2)*sqrt(f_mag/fc_c_mag)*2^iter/4 < exterior < (exp(X)-1/exp(X))*sqrt(f_mag/fc_c_mag)*2^iter
        //  (1-1/exp(X)^2)/X*E/8 < exterior < (exp(X)-1/exp(X))/X*E/2
        //  exp(X)=1+A*X 1/exp(X)=1-B*X 1/exp(X)^2=1-2*C*X, A,B,C~1
        //  C*ln(f_mag)*sqrt(f_mag/fc_c_mag)/4 < exterior < (A+B)/2*ln(f_mag)*sqrt(f_mag/fc_c_mag)
        //  C*E/4 < exterior < (A+B)/2*E
        //  A=(exp(X)-1)/X   B=(1-1/exp(X))/X   (A+B)/2=(exp(X)-1/exp(X))/X/2    C=(1-1/exp(X)^2)/X/2
        //  (A+B)/2 = 1 + x^2/6 + x^4/120 + x^6/5040 + x^8/362880 + x^10/39916800 + x^12/6227020800 + O(x^13)
        //  C = 1 - x + (2 x^2)/3 - x^3/3 + (2 x^4)/15 - (2 x^5)/45 + O(x^6)
        //  assuming f_mag<10000^2, approx up to x^2 should be accurate to 1e-20 with iter>26
        //  1 should be accurate to 1e-20 with iter>71

        //again sinh(G)/(2*exp(G)*G') < exterior < 2*sinh(G)/G'   G=ln(f)/2^iter   G'=f'/(2^iter*f)       sinh(G)=(exp(G)-exp(-G))/2
        // (exp(ln(f)/2^iter)-exp(-ln(f)/2^iter))/2/(2*exp(ln(f)/2^iter)*f'/(2^iter*f)) < exterior < (exp(G)-exp(-G))/f'*(2^iter*f)
        // X=ln(f)/2^iter (=G)   E=2*ln(f)*f/f'
        // (1-1/exp(X)^2)*(2^iter)*f/f'/4 < exterior < (exp(G)-exp(-G))/f'*(2^iter*f)
        // (1-1/exp(X)^2)/X/2*E/4 < exterior < (exp(X)-exp(-X))/X/2*E

        //f0(z)=phi(fc(antiphi(z)))    z=phi(fc(antiphi(u)))
        //f0(f0(u))=phi(fc(fc(antiphi(u))))    u=phi(fc(antiphi(v)))
        //f0(f0(f0(v)))=phi(fc(fc(fc(antiphi(v)))))  phi(large)=large
        //f0(f0(f0(v)))=fc(fc(fc(antiphi(v))))    v=phi(w)
        //f0(f0(f0(phi(w))))=fc(fc(fc(w)))
        //phi(w)^(2^iter)=fc(fc(fc(w)))
        //phi(w)=fc(fc(fc(w)))/(2^iter)
        //G(w)=ln(phi(w))=ln(fc(fc(fc(w))))/(2^iter)
        //G=ln(f)/(2^iter)   G'=f'/f/(2^iter)
        //sinh(G)/(2*exp(G)*G') < exterior < 2*sinh(G)/G'
        //approx RHS: exp(G)/G'=exp(ln(f)/(2^iter))/f'*f*2^iter
        //(1-1/exp(G)^2)*f/f'*(2^iter)/4 < exterior < (exp(G)-exp(-G))/f'*f*(2^iter)
        double fm=mandelData.f.getMag_double();
        double fcm=mandelData.fc_c.getMag_double();
        double x=log(fm);
        mandelData.store->exterior_hits=x*sqrt(fm/fcm);
        mandelData.store->exterior_avoids=mandelData.store->exterior_hits*0.25;
        if (mandelData.store->iter>71)
        { }
        else
        {
          double G=ldexp(x, -1-mandelData.store->iter);
          if (mandelData.store->iter>26)
          {
            mandelData.store->exterior_hits+=G*G/6*mandelData.store->exterior_hits;
            mandelData.store->exterior_avoids+=G*(G*2/3-1)*mandelData.store->exterior_avoids;
          }
          else
          {
            //double ex=exp(x);
            //currentData.store->exterior_hits*=(ex-1/ex)/x/2;
            //currentData.store->exterior_avoids*=(1-1/(ex*ex))/x/2;
            mandelData.store->exterior_hits=ldexp((exp(G)-1/exp(G))*sqrt(fm/fcm), mandelData.store->iter);
            mandelData.store->exterior_avoids=ldexp((1-1/exp(2*G))*sqrt(fm/fcm), mandelData.store->iter-2);
          }
        }
      }
      //already there currentData.fc_c.assign(&fc_c);
      mandelData.store->period=mandelData.store->lookper_lastGuess; //preliminary
      if (mandelData.store->period<1)
        mandelData.store->period=1;
      break;
    };
    double fc_c_mag=mandelData.fc_c.getMag_double();
    if (fc_c_mag>1e57)
    {
      mandelData.store->rstate=MandelPointStore::ResultState::stBoundary;
      mandelData.store->exterior_avoids=0;
      mandelData.store->exterior_hits=0;
      break;
    };
    double fz_c_mag=mandelData.fz_c_mag.toDouble();
    if (fz_c_mag>1e60)
    {
      mandelData.store->rstate=MandelPointStore::ResultState::stDiverge;
      mandelData.store->exterior_avoids=0;
      mandelData.store->exterior_hits=0;
      break;
    };
    //TODO: similar to eval_until_bailout
    //fc_c:=2*f*fc_c+1
    mandelData.fc_c.mul(mandelData.f);
    mandelData.fc_c.lshift(1);
    mandelData.fc_c.re.add_double(1);
    /* TODO: copy test here from above?
    fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>LARGE_FLOAT2)
    {
      currentData.state=MandelPoint::State::stBoundary;
      break;
    };*/
    //const MandelMath::number<BASE> *f_mag2=currentData.f.getMag_tmp(&tmp);
    tmp.tmp2.assign(*mandelData.f.getMag_tmp());
    tmp.tmp2.lshift(2);
    //f'=2*f'*f, f'_mag=4*f'_mag*f_mag
    mandelData.fz_c_mag.mul(tmp.tmp2); //TODO: can use f_mag from above? would need storage, tmp can't survive this long
    //currentData.fz_c_mag.lshift(2);
    mandelData.sure_fz_mag.mul(tmp.tmp2); //TODO: can use f_mag from above? would need storage, tmp can't survive this long
    //currentData.sure_fz_mag.lshift(2);
    mandelData.lookper_totalFzmag.mul(tmp.tmp2);
    //currentData.lookper_totalFzmag.lshift(2);
    //f:=f^2+c
    mandelData.f.sqr();
    mandelData.f.add(currentParams.mandel.c);
    //currentData.iter++;

    //from here, iter should be iter+1

    const MandelMath::number<BASE> *f_mag3=mandelData.f.getMag_tmp();
    if (!mandelData.near0m.isle(*f_mag3)) //f_mag<near0m
    {
      mandelData.store->near0iter_1=mandelData.store->iter+2;
      mandelData.near0m.assign(*f_mag3);
      mandelData.sure_fz_mag.zero(1);
    };
    const MandelMath::number<BASE> *f_distfirst=mandelData.f.dist2_tmp(currentParams.mandel.first_z);
    if (!mandelData.nearzm.isle(*f_distfirst))
    {
      mandelData.store->nearziter_0=mandelData.store->iter+1;
      mandelData.nearzm.assign(*f_distfirst);
    };

    //lookper_nearr = "near(est) to lookper_startf"
    const MandelMath::number<BASE> *lpdiff=mandelData.lookper_startf.dist2_tmp(mandelData.f);
    int new_guess_now_=0; //0..not, 1..maybe, 2..exact match
    switch (MandelMath::strong_ordering_cast(lpdiff->compare(mandelData.lookper_nearr_dist))) //|f-r|<best
    {
      case MandelMath::less:
      {
        mandelData.lookper_nearr.assign(mandelData.f);
        mandelData.lookper_nearr_dist.assign(*lpdiff);
        mandelData.store->lookper_nearr_dist_touched=false;
        mandelData.store->lookper_prevGuess=mandelData.store->lookper_lastGuess;
        mandelData.store->lookper_lastGuess=(mandelData.store->iter+1-mandelData.store->lookper_startiter);
        new_guess_now_=1;
      } break;
      case MandelMath::equal:
      {
        if (mandelData.f.isequal(mandelData.lookper_nearr))
        { //cycle from nearr - may be earlier than cycle from startf
          int prev_nearr=mandelData.store->lookper_lastGuess+mandelData.store->lookper_startiter;
          int next_guess=mandelData.store->iter+1-prev_nearr;
          if (next_guess>mandelData.store->lookper_lastGuess)
          {
            mandelData.store->lookper_prevGuess=mandelData.store->lookper_lastGuess;
          };
          mandelData.store->lookper_lastGuess=next_guess;
          new_guess_now_=2;
        }
        else if (lpdiff->is0())
        { //cycle from startf
          int next_guess=MandelMath::gcd(mandelData.store->iter+1-mandelData.store->lookper_startiter, mandelData.store->nearziter_0);
          if (next_guess>mandelData.store->lookper_lastGuess)
          {
            mandelData.store->lookper_prevGuess=mandelData.store->lookper_lastGuess;
          };
          mandelData.store->lookper_lastGuess=next_guess;
          new_guess_now_=2;
        }
        else if (!mandelData.store->lookper_nearr_dist_touched) //we need to stop increasing lastGuess and restart search
        {                                                   //and also retest period once after new lookper_start
          mandelData.lookper_nearr.assign(mandelData.f);
          mandelData.lookper_nearr_dist.assign(*lpdiff);
          mandelData.store->lookper_nearr_dist_touched=true;
          mandelData.store->lookper_prevGuess=mandelData.store->lookper_lastGuess;
          mandelData.store->lookper_lastGuess=(mandelData.store->iter+1-mandelData.store->lookper_startiter);
          new_guess_now_=1;
        }
        else
        {        //Misiurewicz, e.g. c=i: i^2+i=-1+i, (-1+i)^2+i=-i, (-i)^2+i=-1+i
          nop(); //or just lucky hit of same dist of a few eps2
        }
      } break;
    };

#if SUREHAND_CHOICE>0
    int sure_cycle=(currentData.store->iter+1)/currentData.store->near0iter;
    if ((currentData.store->iter+1)==sure_cycle*currentData.store->near0iter) //iter mod near0==0
    {
      bool exact_match=sure_cycle>1 && currentData.sure_startf.isequal(&currentData.f);
      if ((sure_cycle==1 &&
          currentWorker->toDouble(currentData.fz_c_mag.ptr)<1) ||
          (sure_cycle>1 &&
           //MandelMath::is2kof(currentData.store->iter+1, currentData.store->near0iter) &&
           (sure_cycle&(sure_cycle-1))==0 &&
            //(currentData.store->iter+1)%currentData.store->near0iter==0 &&
            currentWorker->toDouble(currentData.sure_fz_mag.ptr)<1) ||
          exact_match)
      {
        int testperiod=currentData.store->near0iter;
        currentData.store->surehand=testperiod; //does happen that it's increased and that's right: bulb (1/2)*(1/5)

        //ideally, we'd .zero at nearest*2^k and test at nearest*(2^k+1)
        //but we can do both at nearest*2^k
        //if ((currentData.store->iter+1)==currentData.store->near0iter)
        currentData.sure_fz_mag.zero(1);

        int foundperiod;
        foundperiod=periodCheck(testperiod, &currentParams.c, &currentData.f, exact_match); //updates iter, f, f_c, root
        if (foundperiod>0)
        {
          currentData.store->rstate=MandelPointStore::ResultState::stPeriod2;
          currentData.store->period=foundperiod;
          currentData.store->newton_iter=newtres.cyclesNeeded;
          currentData.store->period=estimateInterior(foundperiod, &currentParams.c, &currentData.root);
          if (currentWorker->isl0(interior.inte_abs.ptr))
            currentData.store->rstate=MandelPointStore::ResultState::stMisiur;
          else
          {
            currentData.store->interior=interior.inte_abs.toDouble();
            currentData.fz_r.assign(&interior.fz); //d/dz F_c(r)
            if (testperiod!=currentData.store->period)
              currentData.store->rstate=MandelPointStore::ResultState::stPeriod3;
          }
          //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
          //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
          break;
        };
      };
      currentData.sure_startf.assign(&currentData.f);
      if (currentParams.breakOnNewNearest)
      {
        currentData.store->iter++;
        break;
      }
    };
#endif
    //int jp=currentData.store->nearziter_0;
    if (mandelData.store->iter+2>mandelData.store->near0iter_1 &&
        mandelData.store->near0iter_1!=mandelData.store->nearziter_0)
      nop(); //may be 1 period off, probably because of rounding (seen nearziter_0==near0iter_1+period)
    int jp=mandelData.store->near0iter_1; //trying to abandon nearziter
#if SUREHAND_CHOICE!=1
    if (//(currentData.store->lookper_lastGuess>0) &&
        new_guess_now_>1 ||
        (
          new_guess_now_>0 && //(currentData.store->lookper_lastGuess==(currentData.store->iter+1-currentData.store->lookper_startiter)) && //just found new guess
#if USE_GCD_FOR_CHECKPERIOD
#else
          ((jp % mandelData.store->lookper_lastGuess)==0) && //  period divides nearest, that's a fact
#endif
          (mandelData.store->iter>=3*jp)))  //speedup - don't check period too eagerly
    {
#if USE_GCD_FOR_CHECKPERIOD
      int testperiod=MandelMath::gcd(currentData.near0iter, currentData.lookper_lastGuess);//currentData.lookper_lastGuess
#else
      int testperiod=mandelData.store->lookper_lastGuess;
#endif
      int foundperiod=-1;
      //just assigned, too late to test here if (currentData.f.isequal(&eval.lookper_nearr) && !currentData.f.isequal(&currentData.lookper_startf))
      //  nop();
      if (!mandelData.f.isequal(mandelData.lookper_nearr) && mandelData.f.isequal(mandelData.lookper_startf))
        nop();
      if (mandelData.f.isequal(mandelData.lookper_startf))
      { //exact match - misiurewicz or converged after too many steps or just a lucky hit
        int testperiod2=mandelData.store->iter+1-mandelData.store->lookper_startiter;//currentData.store->lookper_lastGuess;
        if (testperiod2!=mandelData.store->lookper_lastGuess)
          nop(); //at preperiodic, test=2 last=1 we want last=2 at 0+i but last=1 at -2
        //testperiod=testperiod2; //gcd works but I want to debug periodCheck with large period
        testperiod=MandelMath::gcd(testperiod2, mandelData.store->nearziter_0);
        if (testperiod<mandelData.store->nearziter_0)
          nop(); //happens a lot
        foundperiod=periodCheck(testperiod, currentParams.mandel.c, mandelData.lookper_nearr, true);
        if (foundperiod<0) //exact match but repelling
          foundperiod=testperiod;
        //should not do I think, for misiur or period currentData.store->lookper_lastGuess=foundperiod;
        //done in periodCheck currentData.root.assign(&currentData.f);
        //TODO: still needs period cleanup... I think. Near 0+0I
      }
      else if (mandelData.lookper_totalFzmag.toDouble()<MAGIC_MIN_SHRINK)
      {
        foundperiod=periodCheck(testperiod, currentParams.mandel.c, mandelData.lookper_nearr, false); //updates iter, f, f_c, root
        /* does not clean period any more, so no point in calling it
        if ((foundperiod>0) && (foundperiod<testperiod))
        {
          //complex root(currentWorker, &currentData.root_re, &currentData.root_im, true);
          foundperiod=estimateInterior(foundperiod, &currentParams.c, &currentData.root); //-4.7e-22
            //foundperiod=-1; //the cycle can be exact but |f_z|>1 due to mistaken period or misplaced (rounding err) root
        }*/
      };
      if (foundperiod>0)
      {
        mandelData.store->rstate=MandelPointStore::ResultState::stPeriod2;
        mandelData.store->period=foundperiod;
        mandelData.store->newton_iter=newtres.cyclesNeeded;
        mandelData.store->period=estimateInterior(foundperiod, currentParams.mandel.c, mandelData.root);
        if (interior.inte_abs.isl0())
          mandelData.store->rstate=MandelPointStore::ResultState::stMisiur;
        else
        {
          mandelData.store->interior.hits=interior.inte_abs.toDouble();
          mandelData.store->interior.hits_re=interior.inte.re.toDouble();
          mandelData.store->interior.hits_im=interior.inte.im.toDouble();
          mandelData.fz_r.assign(interior.fz); //d/dz F_c(r)
          if (testperiod!=mandelData.store->period)
            mandelData.store->rstate=MandelPointStore::ResultState::stPeriod3;
        }
        //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
        //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
        break;
      };
      if (currentParams.breakOnNewNearest)
      {
        mandelData.store->iter++;
        break;
      }
    };
#endif //SUREHAND_CHOICE
  }
  //data.state=MandelPoint::State::stMaxIter;
  if (!mandelData.store->has_fc_r && currentParams.want_fc_r &&
      ((mandelData.store->rstate==MandelPointStore::ResultState::stPeriod2) ||
       (mandelData.store->rstate==MandelPointStore::ResultState::stPeriod3)))
  {
    loope.eval2(mandelData.store->period, currentParams.mandel.c, mandelData.root, false);
    mandelData.fc_c.assign(loope.f_c);
    mandelData.store->has_fc_r=true;
  };
  if (currentParams.want_extangle &&
      mandelData.store->rstate==MandelPointStore::ResultState::stOutside)
  {
    if (mandelData.store->iter>10000)
      mandelData.extangle.zero(ExtAngle::SPECIAL_VALUE_DEEP);
    else
      extangle.computeMJ(&mandelData.extangle, true, mandelData.store->iter, currentParams.mandel.c, currentParams.mandel.first_z);
      //extangle.compute3(&currentData.extangle, currentData.store->iter, &currentParams.c, &tmp);
    mandelData.store->rstate=MandelPointStore::ResultState::stOutAngle;
  };
}

template <typename BASE>
void MandelEvaluator<BASE>::evaluateJulia()
{
  /*{
    eval.near0fmag.assign(*currentData.near0f.getMag_tmp(&tmp)); //TODO: update on changing near0f
  }*/

  for (; (juliaData.store->iter<currentParams.maxiter) &&
         (juliaData.store->rstate==JuliaPointStore::ResultState::stUnknown) &&
         this->workIfEpoch==this->busyEpoch;
       juliaData.store->iter++)
  {
    //can't really find period from nearz and near0
    //so we use the result of mandel()
    if (juliaData.store->iter==0)
    {
      juliaData.store->lookper_startiter=juliaData.store->iter;
      juliaData.lookper_fz.zero(1.0, 0);
      juliaData.lookper_distr.assign(juliaData.f);
      juliaData.lookper_distr.sub(currentParams.julia.root);
      //newt.newtX.assign(&currentParams.julia.root0);
      //newt.newtX.chs();
      //juliaData.nearmrm.assign(*newt.newtX.dist2_tmp(&currentParams.julia.first_z, &tmp));
      //juliaData.nearmrm.zero(0);
      //juliaData.nearmrm.assign(*currentParams.julia.first_z.getMag_tmp(&tmp));
      juliaData.nearmrm.assign(*currentParams.julia.c.dist2_tmp(currentParams.julia.first_z));
      /*const MandelMath::number<BASE> *distm=newt.newtX.dist2_tmp(&currentParams.julia.first_z, &tmp);
      if (distm->isle(juliaData.nearmrm))
      {
        juliaData.nearmrm.assign(*distm);
        juliaData.store->nearmriter=1;
      };*/
      juliaData.nearmr_f.assign(currentParams.julia.first_z);
      juliaData.nearmr_fz.assign(currentParams.julia.alpha);//zero(1, 0);
      juliaData.store->nearmriter=currentParams.julia.period; //ideally 0 but need to avoid discontinuity at r1
      juliaData.store->nearmr_difficult=false;
    };
    const MandelMath::number<BASE> *f_mag1=juliaData.f.getMag_tmp();
    if (!f_mag1->isle(currentParams.julia.bailout))
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stOutside;
      //theory says the relative error in estimate is less than 3/bailout for large bailout
      //so lets move out a bit
      julia_until_bailout();//, &currentData.f, &currentData.fc_c); //may switch state to stBoundary
      if (juliaData.store->rstate!=JuliaPointStore::ResultState::stOutside)
      {
        //currentWorker->zero(&currentData.exterior_avoids, 0);
        //currentWorker->zero(&currentData.exterior_hits, 0);
        juliaData.store->exterior.zero(0);
      }
      else
      { //from https://en.wikipedia.org/wiki/Julia_set#Using_DEM/J
        //standard exterior distance estimate
        double fm=juliaData.f.getMag_double();
        double fcm=juliaData.fz_z_mag.toDouble();
        //double fcm=juliaData.fz_z.getMag_double();
        double x2=std::log(fm);//ln(|f|)*2
        //from https://www.evl.uic.edu/hypercomplex/html/book/book.pdf par 3.3

        // sinh(G(z))/2/exp(G(z))/G'(z) <= dist <= 2*sinh(G(z))/G'(z)    exact
        //     phi~z^(1/2^iter)  G=ln(phi(z)) ~ ln(z)/2^iter      G' ~ z'/z/2^iter
        // sinh(ln(z)/2^iter)/2/z^(1/2^iter)/z'*z*2^iter <= dist <= 2*sinh(ln(z)/2^iter)/z'*z*2^iter    with G from high iter
        //     sinh(x)=(exp(x)-exp(-x))/2 ~ (1+x-1+x)/2 ~ x
        // (ln(z)/2^iter)/2/z^(1/2^iter)/z'*z*2^iter <= dist <= 2*ln(z)/2^iter/z'*z*2^iter     with sinh(x)~x
        // ln(z)*z/z'/z^(1/2^iter)/2 <= dist <= 2*ln(z)*z/z'
        // ln(z)*z/z'/2 <= dist <= ln(z)*z/z'*2      with z^(1/2^iter)~1 happens close to fractal because high iter, small z

        juliaData.store->exterior.verify_mixed=x2*sqrt(fm/fcm);
        juliaData.store->exterior.verify_empty=juliaData.store->exterior.verify_mixed*0.25;
        if (juliaData.store->iter>71)
        { }
        else
        {
          double G=ldexp(x2, -1-juliaData.store->iter); //log(fm)/2/2^iter
          if (juliaData.store->iter>26)
          {
            juliaData.store->exterior.verify_mixed+=G*G/6*juliaData.store->exterior.verify_mixed;
            juliaData.store->exterior.verify_empty+=G*(G*2/3-1)*juliaData.store->exterior.verify_empty;
          }
          else
          {
            //double ex=exp(x);
            //currentData.store->exterior_hits*=(ex-1/ex)/x/2;
            //currentData.store->exterior_avoids*=(1-1/(ex*ex))/x/2;
            //(exp(x*(2^-iter))-1/exp(x*(2^-iter)))/x*2^iter   (exp(a)-1/exp(a))/a~(1+a-1+a)/a=2
            juliaData.store->exterior.verify_mixed=ldexp((exp(G)-1/exp(G))*sqrt(fm/fcm), juliaData.store->iter);
            juliaData.store->exterior.verify_empty=ldexp((1-1/exp(2*G))*sqrt(fm/fcm), juliaData.store->iter-2);
          }
        }

#if 1 //using glue at near0
        bool use_glue=(currentParams.julia.patchSizeExterior>0 && juliaData.near0m.toDouble()<currentParams.julia.patchSizeExterior);
        //1. find d.e. at near0+1
        int iter_p2=juliaData.store->iter-juliaData.store->near0iter_1;;
        double fcm_p2=juliaData.since0fzm.toDouble();
        double fm_p2=juliaData.f.getMag_double();
        double x2_p2=std::log(fm_p2);
        double exterior_p2_mixed=x2_p2*sqrt(fm_p2/fcm_p2);
        double exterior_p2_empty=exterior_p2_mixed*0.25;
        if (iter_p2>71)
        { }
        else
        {
          double G=ldexp(x2_p2, -1-iter_p2); //log(fm)/2/2^iter
          if (iter_p2>26) //G^4/120=1e-19 for G<6e-5 for iter_p2>19; G^3/3=1e-19 for iter_p2>=26
          {
            exterior_p2_mixed+=G*G/6*exterior_p2_mixed; //2G+G^3/3+G^5/60+...=2G*(G^2/6+ G^4/120+...)
            exterior_p2_empty+=G*(G*2/3-1)*exterior_p2_empty; //2G-2G^2+4G^3/3 -2G^4/3...=2G*(1-G+2G^2/3 -G^3/3...)
          }
          else
          {
            //double ex=exp(x);
            //currentData.store->exterior_hits*=(ex-1/ex)/x/2;
            //currentData.store->exterior_avoids*=(1-1/(ex*ex))/x/2;
            //(exp(x*(2^-iter))-1/exp(x*(2^-iter)))/x*2^iter   (exp(a)-1/exp(a))/a~(1+a-1+a)/a=2
            exterior_p2_mixed=ldexp((exp(G)-1/exp(G))*sqrt(fm_p2/fcm_p2), iter_p2);
            exterior_p2_empty=ldexp((1-1/exp(2*G))*sqrt(fm_p2/fcm_p2), iter_p2-2);
          }
        }

        //2. walk 1 iteration backward
        double near0=std::sqrt(juliaData.near0m.toDouble());
        juliaData.store->exterior.none_mixed=exterior_p2_mixed/2/near0; //preview of simple scaling
        juliaData.store->exterior.none_empty=exterior_p2_empty/2/near0;

        // r=dist estimate of previous step "d.e.before"
        // (x+r exp(i t))^2+c=x*x+2*x*r exp(i t)+r*r*exp(2 i t)+c
        //     2|x|r+r^2=d.e.after
        //     for r << |x| :
        //     d.e.before=d.e.after/2/|x|
        MandelMath::real_double_quadratic(&juliaData.store->exterior.always_mixed, 1, near0, -exterior_p2_mixed); //actual quad step
        MandelMath::real_double_quadratic(&juliaData.store->exterior.always_empty, 1, near0, -exterior_p2_empty);

        //trying to track the sinh during glue step, seems to actually hurt
        double glue_scale, glue_scale_empty;
        if (iter_p2>68) //(x2>>68)+1 ~ 1
        {
          glue_scale=1;
          glue_scale_empty=1;
        }
        else if (iter_p2>35) //with z<37 (bailout 10^8), z^2*(1/a^2)/6=1e-6 for a=2^14; 1e-19 for a=2^36 (actually should check 2nd term but who cares)
        { //sinh(x/a)/sinh(x/b)*a/b=1+z^2*(1/a^2-1/b^2)/6+O(z^4);
          double glue_g1=ldexp(x2_p2, -1-(iter_p2));
          double glue_g2=ldexp(x2_p2, -1-(iter_p2+1));
          glue_scale=1+(glue_g1*glue_g1-glue_g2*glue_g2)/6;
          double glue_t1=-ldexp(1, -iter_p2-1);//ldexp(1, -iter_p2-1)-ldexp(1, -iter_p2);  err~(x2_p2/2*glue_t1)^2/2 -> glue_t1<2^-36 for err~1e-19
          glue_scale_empty=1+x2_p2/2*glue_t1;
        }
        else
        {
          double glue_g1=ldexp(x2_p2, -1-(iter_p2));
          double glue_g2=ldexp(x2_p2, -1-(iter_p2+1));
          glue_scale=std::sinh(glue_g1)/std::sinh(glue_g2)/2;
          double glue_t1=-ldexp(1, -iter_p2-1);//ldexp(1, -iter_p2-1)-ldexp(1, -iter_p2);
          glue_scale_empty=exp(x2_p2/2*glue_t1);
        }
        glue_scale_empty*=glue_scale;
        glue_scale_empty=1/glue_scale_empty;
        glue_scale=1/glue_scale;

        juliaData.store->exterior.none_mixed*=glue_scale;
        juliaData.store->exterior.none_empty*=glue_scale_empty;
        //did I forget to scale always_ or was it on purpose? who knows...
        //for small iter_p2~1..5 it sucks both ways
        //for large iter_p2 it's 0.999 and similar
        //seems best to scale empty but not mixed
        //juliaData.store->exterior.always_mixed*=glue_scale;
        //juliaData.store->exterior.always_empty*=glue_scale_empty;

        if (!use_glue)
        {
          juliaData.store->exterior.condi_mixed=juliaData.store->exterior.none_mixed;
          juliaData.store->exterior.condi_empty=juliaData.store->exterior.none_empty;
          juliaData.store->exterior.blend_mixed=juliaData.store->exterior.none_mixed;
          juliaData.store->exterior.blend_empty=juliaData.store->exterior.none_empty;
        }
        else
        {
          juliaData.store->exterior.condi_mixed=juliaData.store->exterior.always_mixed;
          juliaData.store->exterior.condi_empty=juliaData.store->exterior.always_empty;
#if 0 //blend by varying "a"
          double weight=juliaData.near0m.toDouble()/currentParams.patchSizeExterior;
          MandelMath::real_double_quadratic(&juliaData.store->exterior.blend_mixed, 1-(weight), near0, -exterior_p2_mixed); //actual quad step
          MandelMath::real_double_quadratic(&juliaData.store->exterior.blend_empty, 1-(weight), near0, -exterior_p2_empty);

#elif 1 //(p2/(2*near0)*(near0m/patch)+p1*(1-near0m/patch)) = (p2/2*near0+p1*(patch-near0m))/patch
//looks best of the 3
          double help=currentParams.julia.patchSizeExterior-juliaData.near0m.toDouble();
          //MandelMath::real_double_quadratic(&exterior_p1_hits, 1, near0, -exterior_p2_hits); //actual quad step
          //MandelMath::real_double_quadratic(&exterior_p1_avoids, 1, near0, -exterior_p2_avoids);
          juliaData.store->exterior.blend_mixed=(exterior_p2_mixed/2*near0*glue_scale+juliaData.store->exterior.always_mixed*help)/currentParams.julia.patchSizeExterior;
          juliaData.store->exterior.blend_empty=(exterior_p2_empty/2*near0*glue_scale_empty+juliaData.store->exterior.always_empty*help)/currentParams.julia.patchSizeExterior;

#else   //(p2/(2*near0)*(near0m/patch)^2+p1*(1-(near0m/patch)^2)) = (p2/2*near0*near0m+p1*(patch^2-near0m^2))/patch^2
          double help=MandelMath::sqr_double(currentParams.patchSizeExterior)-MandelMath::sqr_double(juliaData.near0m.toDouble());
          juliaData.store->exterior.blend_mixed=(exterior_p2_mixed/2*near0*juliaData.near0m.toDouble()*glue_scale + juliaData.store->exterior.always_mixed*help)/(currentParams.patchSizeExterior*currentParams.patchSizeExterior);
          juliaData.store->exterior.blend_empty=(exterior_p2_empty/2*near0*juliaData.near0m.toDouble()*glue_scale_empty + juliaData.store->exterior.always_empty*help)/(currentParams.patchSizeExterior*currentParams.patchSizeExterior);
#endif
        }

        /*juliaData.store->exterior.always_hits*=glue_scale;
        juliaData.store->exterior.always_avoids*=glue_scale_avoids;
        juliaData.store->exterior.condi_hits*=glue_scale;
        juliaData.store->exterior.condi_avoids*=glue_scale_avoids;
        juliaData.store->exterior.blend_hits*=glue_scale;
        juliaData.store->exterior.blend_avoids*=glue_scale_avoids;*/

        //3. scale to the beginning
        //  full d.e. at iter==I with end at z, iter K; z'=z'[0..K]
        //  avoids = sinh(ln(z)/2^(K-I))/z'[I..K]*z*2^(K-I)/z^(1/2^(K-I))/2
        //  hits = 2*sinh(ln(z)/2^(K-I))/z'[I..K]*z*2^(K-I)
        // avoids[J]/avoids[I]=sinh(ln(z)/2^(K-J))/z'[J..K]*z*2^(K-J)/z^(1/2^(K-J))/2 / ( sinh(ln(z)/2^(K-I))/z'[I..K]*z*2^(K-I)/z^(1/2^(K-I))/2 )
        //                    =sinh(ln(z)/2^(K-J))/sinh(ln(z)/2^(K-I)) *z'[I..K]/z'[J..K] *2^(K-J)/2^(K-I) *z^(1/2^(K-I))/z^(1/2^(K-J))
        //                    =sinh(ln(z)/2^(K-J))/sinh(ln(z)/2^(K-I)) *z'[I..J] *2^(I-J) *z^(1/2^(K-I)-1/2^(K-J))
        //                    =sinh(ln(z)/2^(K-J))/sinh(ln(z)/2^(K-I)) *2^(I-J) *z'[I..J] *z^(2^(I-K)-2^(J-K))
        //iter_p2--;
        double phase1_scale, phase1_scale_avoids;
        if (iter_p2>68)
        {
          phase1_scale=std::sqrt(juliaData.near0fzm.toDouble());//     *juliaData.near0m.toDouble())*2;
          phase1_scale_avoids=1;
        }
        else if (iter_p2>35)
        {
          double avoids_g1=ldexp(x2_p2, -1-(iter_p2+1));
          double avoids_g2=ldexp(x2_p2, -1-juliaData.store->iter);
          phase1_scale=1+(avoids_g1*avoids_g1-avoids_g2*avoids_g2)/6;
          phase1_scale*=std::sqrt(juliaData.near0fzm.toDouble());//     *juliaData.near0m.toDouble())*2;
          double avoids_t1=ldexp(1, -juliaData.store->iter)-ldexp(1, -iter_p2-1);
          phase1_scale_avoids=1+x2_p2/2*avoids_t1;
        }
        else
        {
          double avoids_g1=ldexp(x2_p2, -1-(iter_p2+1));
          double avoids_g2=ldexp(x2_p2, -1-juliaData.store->iter);
          phase1_scale=std::sinh(avoids_g1)/std::sinh(avoids_g2)*ldexp(1, 1-juliaData.store->near0iter_1);
          phase1_scale*=std::sqrt(juliaData.near0fzm.toDouble());//     *juliaData.near0m.toDouble())*2;
          double avoids_t1=ldexp(1, -juliaData.store->iter)-ldexp(1, -iter_p2-1);
          phase1_scale_avoids=exp(x2_p2/2*avoids_t1);
        }
        phase1_scale_avoids=1/(phase1_scale*phase1_scale_avoids);
        phase1_scale=1/phase1_scale;
        juliaData.store->exterior.none_mixed*=phase1_scale;
        juliaData.store->exterior.none_empty*=phase1_scale_avoids;
        juliaData.store->exterior.always_mixed*=phase1_scale;
        juliaData.store->exterior.always_empty*=phase1_scale_avoids;
        juliaData.store->exterior.condi_mixed*=phase1_scale;
        juliaData.store->exterior.condi_empty*=phase1_scale_avoids;
        juliaData.store->exterior.blend_mixed*=phase1_scale;
        juliaData.store->exterior.blend_empty*=phase1_scale_avoids;
        //juliaData.store->exterior_hits=juliaData.store->exterior_hits/exterior_hits-1;
        //juliaData.store->exterior_avoids=juliaData.store->exterior_avoids/exterior_avoids-1;
#endif
      }
      break;
    };
    //fzz=2*(f*fzz+fz^2)
    juliaData.fzz_z.mul(juliaData.f);
    eval.fz_r.assign(juliaData.fz_z);
    eval.fz_r.sqr();
    juliaData.fzz_z.add(eval.fz_r);
    juliaData.fzz_z.lshift(1);

    double fz_z_mag1=juliaData.fz_z.getMag_double();
    if (fz_z_mag1>1e57)
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stBoundary;
      juliaData.store->exterior.zero(0);
      break;
    };
    double fz_z_mag2=juliaData.fz_z_mag.toDouble();
    if (fz_z_mag2>1e60)
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stDiverge;
      juliaData.store->exterior.zero(0);
      break;
    };
    //fz_z:=2*f*fz_z
    juliaData.fz_z.mul(juliaData.f);
    juliaData.fz_z.lshift(1);
    //juliaData.fz_z.re.add_double(1);
    /* TODO: copy test here from above?
    fc_c_mag=currentWorker->toDouble(fc_c.getMagTmp());
    if (fc_c_mag>LARGE_FLOAT2)
    {
      currentData.state=MandelPoint::State::stBoundary;
      break;
    };*/
    //const MandelMath::number<BASE> *f_mag2=currentData.f.getMag_tmp(&tmp);
    tmp.tmp2.assign(*juliaData.f.getMag_tmp());
    tmp.tmp2.lshift(2);
    //f'=2*f'*f, f'_mag=4*f'_mag*f_mag
    juliaData.fz_z_mag.mul(tmp.tmp2);
    juliaData.since0fzm.mul(tmp.tmp2);
    /*if (juliaData.biggest_fzm.isle(juliaData.fz_z_mag))
    {
      juliaData.biggest_fzm.assign(juliaData.fz_z_mag);
      juliaData.store->bigfzm_iter=juliaData.store->iter;
    }*/
    //currentData.fz_c_mag.lshift(2);
    juliaData.lookper_fz.mul(juliaData.f);
    juliaData.lookper_fz.lshift(1);
    //currentData.lookper_totalFzmag.lshift(2);
    //f:=f^2+c
    juliaData.f.sqr();
    juliaData.f.add(currentParams.julia.c);
    //currentData.iter++;

    //from here, iter should be iter+1

    //after test for near0 juliaData.nearmr.tap(currentParams.julia, &newt.newtX, &tmp);

    //find iter where f'/f'' is biggest
    bulb.t1.re.assign(*juliaData.fzz_z.getMag_tmp());
    bulb.t1.re.recip();
    bulb.t1.re.mul(juliaData.fz_z_mag);
    if (juliaData.bigfzfzzm.isle(bulb.t1.re))
    {
      juliaData.bigfzfzzm.assign(bulb.t1.re);
      juliaData.store->bigfzfzziter=juliaData.store->iter+1;
    };

    //find cumulative shrink factor S*=1+f_z^2/(f*f_zz)
    bulb.t1.assign(juliaData.fzz_z);
    bulb.t1.mul(juliaData.f);
    bulb.t1.recip();
    bulb.t1.mul(juliaData.fz_z);
    bulb.t1.mul(juliaData.fz_z);
    bulb.t1.re.add_double(1);
    juliaData.shrinkfactor.mul(bulb.t1);

    int need_reset_check=0; //1..check and reset; 2..reset
    if (juliaData.store->near0iter_1==juliaData.store->iter+1)
    {
      need_reset_check=2; //TODO: finding nearciter_0 instead of near0iter_1 could make things easier
      juliaData.since0fzm.zero(1);
    };

    const MandelMath::number<BASE> *f_mag3=juliaData.f.getMag_tmp();
    if (!juliaData.near0m.isle(*f_mag3)) //f.mag<near0.mag
    {
      juliaData.store->near0iter_1=juliaData.store->iter+2;
      juliaData.near0m.assign(*f_mag3);
      juliaData.near0fzm.assign(juliaData.fz_z_mag);
      juliaData.since0fzm.zero(1); //next iteration really
      //need_reset_check=2;
    }
    else if (currentParams.julia.period>0 && juliaData.store->iter-juliaData.store->lookper_startiter==currentParams.julia.period)
    {
      need_reset_check=1;
    }

    juliaData.nearmr.tap(currentParams.julia, &newt.newtX);

    const MandelMath::number<BASE> *f_distfirst=juliaData.f.dist2_tmp(currentParams.julia.first_z);
    if (!juliaData.nearzm.isle(*f_distfirst))
    {
      juliaData.store->nearziter_0=juliaData.store->iter+1;
      juliaData.nearzm.assign(*f_distfirst);
    };

    if (need_reset_check==1)
    {
      if (juliaData.store->iter+1<3*currentParams.julia.period+juliaData.store->near0iter_1 ||
          juliaData.lookper_fz.mag_cmp_1()!=std::strong_ordering::less)
        need_reset_check=2;
    };
    /* has false positives with alpha~1
    if (need_reset_check==1)
    { //test if (f-r) shrinks alpha times over last juliaPeriod
      bulb.t1.assign(&juliaData.f);
      bulb.t1.sub(&currentParams.juliaRoot_);
      bulb.t2.assign(&juliaData.lookper_distr);
      bulb.t2.mul(&currentParams.juliaAlpha, &tmp);
      bulb.t2.recip(&tmp);
      bulb.t1.mul(&bulb.t2, &tmp); // (f-r)/((old_f-r)*alpha)
      bulb.t1.re.add_double(-1);
      if (bulb.t1.getMag_double()>0.1)
        need_reset_check=2;
    };*/
    if (need_reset_check==1)
    {
      bulb.t1.assign(juliaData.lookper_fz);
      bulb.t1.sub(currentParams.julia.alpha);
      double safety=std::sqrt(bulb.t1.getMag_double())/-currentParams.julia.alpha.getMag1_tmp()->toDouble(); // |fz-alpha|/(1-|alpha|)
      if (safety>1.0) //1.0 seems enough //0.5 seems enough
        need_reset_check=2;
    };

    if (need_reset_check==1)
      {
        juliaData.store->rstate=JuliaPointStore::ResultState::stPeriod2;
        loope.f.assign(juliaData.f);
        loope.f_z.assign(juliaData.fz_z);
        loope.f_zz.assign(juliaData.fzz_z);
        newt.newtX.zero(0, 0);
        //need (iter-near0iter_1)%juliaPeriod=0  (n-i)%p=-(i-n)%p=(p-(i-n)%p)%p=(p-1-(i-n+1)%p)%p=p-1-(i-n+1)%p
        // //if (-1)%2==1 int need_iters=(juliaData.store->near0iter_1-juliaData.store->iter)%juliaPeriod;
        //int need_iters=juliaPeriod-1-(juliaData.store->iter-juliaData.store->near0iter_1+1)%juliaPeriod; //if (-1)%2==-1
        //actually we need (iter+1-near0iter_1)%juliaPeriod=0 because iter is lagging behind juliaData
        int need_iters=currentParams.julia.period-1-(juliaData.store->iter-juliaData.store->near0iter_1)%currentParams.julia.period; //if (-1)%2==-1
        if (need_iters!=0)
           dbgPoint(); //need to advance juliaData.f,fz,fzz too
        if (!loope.eval_zz(need_iters, currentParams.julia.c, newt.newtX, false, false))
        {
          interior.inte.zero(0, 0);
          interior.inte_abs.zero(0);
        }
        else
        {
          juliaData.store->iter+=need_iters;
          //phi=alpha^iter/(f-r)    alpha1=2*r0*alpha0  phi0=alpha0/(f0-r0)  phi1=2*f0*alpha0/(f1-r1)   f1=f0^2+c  r1=r0^2+c
          //phi1=2*r0*alpha0/(f0^2-r0^2)=2*r0*alpha0/(f0-r0)/(f0+r0)=2*r0/(f0+r0)*phi0    2*r0/(f0+r0)=r0/((f0-r0)/2+r0)=1/(1+(f0-r0)/r0/2)=2/(1+q)=1/(1+w)  q=f0/r0   w=(q-1)/2  w2=f0/r0-1
          //phi0'=-alpha0*f0'/(f0-r0)^2   phi1'=-alpha1*f1'/(f1-r1)^2   f1'=2*f0*f0'
          //phi1'=-2*r0*alpha0*2*f0*f0'/(f0+r0)^2/(f0-r0)^2=4*r0*f0/(f0+r0)^2*phi0'=4/(f0/r0+1)/(1+r0/f0)*phi0'=4*q/(1+q)^2*phi0'=(1+w2)/(1+w2/2)^2*phi0'
          //phi1/phi1'=2*r0/(f0+r0)*phi0 / -2/r0/alpha0/2/f0/f0'*(f0+r0)^2*(f0-r0)^2=
          //          =r0/(f0+r0)*phi0 /r0/2/f0*(f0+r0)^2/-alpha0/f0'*(f0-r0)^2=1/2*(f0+r0)/f0*phi0/phi0'
          //phi1/phi1'=1/2*(1+(r0-f0)/f0+f0/f0)*phi0/phi0'=(1-(f0-r0)/f0/2)*phi0/phi0'   =(1+1/q)/2*phi0/phi0'=(1+w)/(1+2*w)*phi0/phi0'

          //phi0''=-alpha0*(f0''/(f0-r0)^2-2*f0'^2/(f0-r0)^3)  f1''=2*(f0'^2+f0*f0'')
          //phi1''=-alpha1*(f1''/(f1-r1)^2-2*f1'^2/(f1-r1)^3) = -2*r0*alpha0*(2*(f0'^2+f0*f0'')/(f0^2-r0^2)^2-2*(2*f0*f0')^2/(f0^2-r0^2)^3)
          //nope phi1''/phi0''=4*r0*((r0^2-f0^2)*(f0*f0''+f0'^2)+4*f0^2*f0'^2) )/((r0+f0)^3*(f0''*(r0-f0)+2*f0'^2))
          //phi1''/(-2*r0*alpha0)=2*(f0'^2+f0*f0'')/(f0^2-r0^2)^2-2*(2*f0*f0')^2/(f0^2-r0^2)^3
          //                     =2*f0'^2/(f0^2-r0^2)^2 + 2*f0*f0''/(f0^2-r0^2)^2 -8*f0^2*f0'^2/(f0^2-r0^2)^3=
          //                     =2*f0*(f0''/(f0^2-r0^2)^2) +4*f0^2*(-2*f0'^2/(f0^2-r0^2)^3) +2*f0'^2/(f0^2-r0^2)^2
          //       phi0'' nowhere to be found
          //phi1''/phi0''=4*r0*( (1-f0/r0)*(r0+f0)*(f0*f0''+f0'^2)+4*f0^2*f0'^2/r0) )/((r0+f0)^3*(f0''*(1-f0/r0)+2*f0'^2/r0))
          //phi1''/phi0''=4*r0*( (1-f0/r0)*((r0+f0)*(f0*f0''+f0'^2)-4*f0*f0'^2)+4*f0*f0'^2) )/((r0+f0)^3*(f0''*(1-f0/r0)+2*f0'^2/r0))
          //phi1''/phi0''=4*r0*( 4*f0*f0'^2 + (1-f0/r0)*((r0+f0)*(f0*f0''+f0'^2)-4*f0*f0'^2)) )/((r0+f0)^3*(2*f0'^2/r0 + (1-f0/r0)*f0''))
          //lim[f0->r0] phi1''/phi0''=8*(f0/r0)/(1+f0/r0)^3
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0-(f0-r0)^2*f0'^2/(f0+r0)/((f0-r0)*f0''-2*f0'^2) )
          //    (f0-r0)^2=(f0+r0)^2-4*f0*r0
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0 -(f0+r0)*f0'^2/((f0-r0)*f0''-2*f0'^2) +4*f0*r0*f0'^2/(f0+r0)/((f0-r0)*f0''-2*f0'^2) )
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0*(1 -f0'^2/((f0-r0)*f0''-2*f0'^2))  +r0*f0'^2*(4*f0/(f0+r0)-1)/((f0-r0)*f0''-2*f0'^2) )
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0*((f0-r0)*f0''-3*f0'^2)  +r0*f0'^2*(3*f0-r0)/(f0+r0) )/((f0-r0)*f0''-2*f0'^2)
          //v=(f0-r0)*f0''/f0'^2
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0*(v-3)  +r0*(3*f0-r0)/(f0+r0) )/(v-2)
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0*(v-2)  +4*f0*r0/(f0+r0) -f0-r0 )/(v-2)
          //phi1''/phi0''=4*r0/(f0+r0)^2 * ( f0  +(f0^2+r0^2)/(f0+r0)/(2-v) )

          //again:
          //phi1''/phi0''=4 r0 ((r0^2+f0^2) - r0 v f0 + (2-v) f0^2)/((2 - v) (r0 + f0)^3)
          //u=f0/r0=1+w   u-1=w
          //phi1''/phi0''=4 * ((1+u^2) - v u + (2-v) u^2)/((2 - v) (1 + u)^3)
          //phi1''/phi0''=(4 (u^2 (v - 3) + u v - 1))/((u + 1)^3 (v - 2))
          //phi1''/phi0''=4 ((3w-v*(1+w)) (2+w) + 4)/((2-v)(2+w)^3)
          interior.alphak.assign(currentParams.julia.alpha);
          //actually, if we get the power off by one it doesn't matter because it cancels on phi/phi'
          interior.alphak.pow_uint((juliaData.store->iter+1)/currentParams.julia.period);

          // phi_z =-alpha^iter*f_z*(f-r)^-2
          // phi_zz=-alpha^iter*(f_zz*(f-r)^-2-2*f_z^2*(f-r)^-3)
          // phi_z/phi_zz=1 / ((f_zz/f_z-2*f_z/(f-r)))
          //estimate root: where phi=inf with phi=a/(x-b)+c phi'=-a/(x-b)^2 phi''=2*a/(x-b)^3
          //phi'/phi''=-a/(x-b)^2/2/a*(x-b)^3=-(x-b)/2
          // -(x-b)=2*phi_z/phi_zz=1/ ((f_zz/f_z/2-f_z/(f-r))) so if we neglect fzz/fz, we're back to (f-r)/f'
          //f_z at estimated r: f_z+f_zz*(x-b) = f_z+2*f_z/ ((1-2*f_z^2/(f-r)/f_zz))
          interior.step_to_root.assign(loope.f_z);
          interior.step_to_root.recip();
          interior.step_to_root.mul(loope.f_zz);
          bulb.t2.assign(loope.f);
          bulb.t2.sub(currentParams.julia.root);
          bulb.t2.recip();
          bulb.t2.mul(loope.f_z);
          interior.step_to_root.lshift(-1);
          interior.step_to_root.sub(bulb.t2);
          interior.step_to_root.recip(); // 1/(f_zz/f_z/2-f_z/(f-r1))

          interior.alphak_other.assign(interior.step_to_root);
          interior.alphak_other.mul(loope.f_zz);
          interior.alphak_other.add(loope.f_z); //verified experimetally against f_z-f_zz*step_to_root

          interior.alphak.assign(interior.alphak_other);


          //actually^2, it's way more complicated
          //interior.alphak.pow_int((juliaData.store->iter-juliaData.store->nearmriter+1)/currentParams.julia.period, &tmp);
          //interior.alphak.pow_int((juliaData.store->iter+1-juliaData.store->near0iter_1)/juliaPeriod, &tmp);
          //interior.alphak.pow_int((juliaData.store->iter+1-juliaData.store->bigfzm_iter)/juliaPeriod, &tmp);
          MandelMath::complex<BASE> &phi=bulb.t1;
          MandelMath::complex<BASE> &phi_z=bulb.t2;
          phi.assign(loope.f);
          phi.sub(currentParams.julia.root);
          phi.recip();
          phi_z.assign(phi);
          phi.mul(interior.alphak); //phi=alpha^iter/(f-r)
          phi_z.mul(phi);
          phi_z.mul(loope.f_z);
          phi_z.chs(); //phi'=-alpha^iter*f'/(f-r)^2
          interior.inte.assign(phi_z);
          interior.inte.recip();
          interior.inte.mul(phi); //d.e.=phi/phi'
          MandelMath::complex<BASE> &w2=bulb.t3;
          MandelMath::complex<BASE> &recipr=newt.laguG;
          MandelMath::complex<BASE> &root=interior.fz;
          root.assign(currentParams.julia.root);
          int didcycles=0;
          for (;didcycles<5*currentParams.julia.period;)
          {
            recipr.assign(root);
            recipr.recip();
            w2.assign(loope.f);
            w2.sub(root);
            w2.mul(recipr);   //(f-r)/r
            double test=w2.getMag_double();
            if (!(test*test>loope.f.re.eps2())) //break if test is nan
              break;
            recipr.assign(w2);
            recipr.lshift(-1);
            recipr.re.add_double(1);
            recipr.recip();
            phi.mul(recipr); //phi*=1/(1+w2/2)
            recipr.sqr();
            phi_z.mul(recipr);
            recipr.assign(w2);
            recipr.re.add_double(1);
            phi_z.mul(recipr); //phi'*=(1+w2)/(1+w2/2)^2 = 1/(1+w2^2/4/(1+w2)) ~ 1
            recipr.assign(w2);
            recipr.re.add_double(1);
            recipr.recip();
            recipr.re.add_double(1);
            recipr.lshift(-1);
            interior.inte.mul(recipr); //inte*=(1+1/(1+w2))/2

            //fzz=2*(f*fzz+fz^2)
            /* completely f...mixed up juliaData, eval and loope
            juliaData.fzz_z.mul(&loope.f, &tmp);
            eval.fz_r.assign(&juliaData.fz_z);
            eval.fz_r.sqr(&tmp);
            juliaData.fzz_z.add(&eval.fz_r);
            juliaData.fzz_z.lshift(1);
            //fz_z:=2*f*fz_z
            juliaData.fz_z.mul(&loope.f, &tmp);
            juliaData.fz_z.lshift(1);
            //f:=f^2+c
            loope.f.sqr(&tmp);
            loope.f.add(&currentParams.c);*/
            loope.eval_zz(1, currentParams.julia.c, newt.newtX, false, false);
            //r:=r^2+c
            root.sqr();
            root.add(currentParams.julia.c);
            didcycles++;
          }
#if 0
          //keep phi/phi'
#elif 0
          //recompute phi/phi'
          interior.inte.assign(&phi_z);
          interior.inte.recip(&tmp);
          interior.inte.mul(&phi, &tmp);
#elif 0
          //linearization through exp(x)-1: phi/phi' -> (exp(phi)-1)/(exp(phi)-1)'=(exp(phi)-1)/(exp(phi)*phi')=(1-exp(-phi))/phi'
          interior.inte.assign(&phi_z);
          interior.inte.recip(&tmp);
          w2.assign(&phi);
          w2.chs();
          w2.exp_approx();
          w2.chs();
          w2.re.add_double(1);
          interior.inte.mul(&w2, &tmp);
#elif 0
          //linearization through exp(-x)-1: phi/phi' -> (exp(-phi)-1)/(exp(-phi)-1)'=(exp(-phi)-1)/(-exp(-phi)*phi')=(exp(phi)-1)/phi'
          interior.inte.assign(&phi_z);
          interior.inte.recip(&tmp);
          w2.assign(&phi);
          w2.exp_approx();
          w2.re.add_double(-1);
          interior.inte.mul(&w2, &tmp);
#elif 0
          //linearization through exp(-|x|)-1: phi/phi' -> (exp(-|phi|)-1)/(exp(-|phi|)-1)'=(exp(|phi|)-1)/|phi'|
          interior.inte.sign(&tmp);
          double phiabs=std::sqrt(phi.getMag_double());
          interior.inte.mul_double((std::exp(phiabs)-1)/std::sqrt(phi_z.getMag_double()));
#elif 0
          //linearization through 1-1/(1+|x|) = 1/(1+1/|x|): phi/phi' -> (1-1/(1+|phi|))/(1-1/(1+|phi|))'=(|phi|+1)*|phi|/|phi'|
          double phiabs=std::sqrt(phi.getMag_double());
          interior.inte.mul_double(1+phiabs);
#elif 1
          //lin(0)=0 lin(eps)=eps lin(phi(r+eps))/(lin(phi(r+eps)))'=1/jas
          //phi(r+eps)~1/eps     and use lin(x)=Rx/(R+x)
          //need lin(1/eps)/(lin(1/eps))'=1/jas   -(1-R*eps)/R=-1/R  R=-jas
          //linearization through Rx/(R+x): phi/phi' -> (1+phi/R)*phi/phi'(x)              near r1, phi*phi/phi' ~ -1; phi/phi'-1/R
          //for first_z in juliaRoot's pool, we want d.e.=1/jas^2
          //more precisely, d.e.[r0]=(1/jas-2*r0)   d.e.[r1]=(1/jas)*(1/jas-2*r0)
          double phiabs=phi.getMag_double();
          //otherwise, scale by one more juliaAlpha. Don't ask why, it's complicated
          //if ((juliaData.store->iter+1)%juliaPeriod !=0)
          //  phiabs*=currentParams.juliaAlpha.getMag_double();

          //scale such that 2*f'*P+f''/2*P^2=R=1/jas^2   ... with f'~0 ... bigger of the 2 roots so find 1/P
          // p=-1/P   2*f'*p-f''/2+p^2/jas^2=0
          double r1de, rfde;
          //f'/f'' seems to converge to a value, what is it?
          //  2*f'/f'' = first_z-first_r
          //rde=2*sqrt(juliaData.fz_z.getMag_double()/juliaData.fzz_z.getMag_double());
          //rde*=currentParams.juliaAlphaShort.getMag_double();
          //MandelMath::real_double_quadratic(&rde, 1/currentParams.juliaAlphaShort.getMag_double(),
          //MandelMath::real_double_quadratic(&rde, interior.alphak.getMag_double()/currentParams.juliaAlphaShort.getMag_double(),
          if (didcycles%currentParams.julia.period!=0)
            dbgPoint(); //who knows what to do... more cycles until it's 0?
          //MandelMath::real_double_quadratic(&rde, interior.alphak.getMag_double()*std::pow(currentParams.juliaAlpha.getMag_double(), didcycles/juliaPeriod)/currentParams.juliaAlphaShort.getMag_double(),
          //                                        sqrt(loope.f_z.getMag_double()),
          //                                        -sqrt(loope.f_zz.getMag_double())/2);
          /*MandelMath::real_double_quadratic(&rde, interior.alphak.getMag_double()*std::pow(currentParams.julia.alpha.getMag_double(), didcycles/currentParams.julia.period)/currentParams.julia.alphaShort.getMag_double(),
                                                  sqrt(loope.f_z.getMag_double()),
                                                  -sqrt(loope.f_zz.getMag_double())/2);*/

          //d.e.[root]=(1/jas)*(1/jas-2*|root0|)
          //rde=1/std::sqrt(currentParams.julia.alphaShort.getMag_double());
          //rde=rde*(rde-2*std::sqrt(currentParams.julia.root0.getMag_double()));
          //rde/=std::sqrt(juliaData.nearmr_fz.getMag_double());

          //actually 1/d.e. is easier
          //1/d.e.[root]=(jas^2)/(1-2*jas*|root0|)
          r1de=currentParams.julia.alphaShort.getMag_double();
          r1de=r1de/(1-2*std::sqrt(r1de*currentParams.julia.root0.getMag_double())); //1/de at r1
#if 1
          rfde=r1de*std::sqrt(juliaData.nearmr_fz.getMag_double()); // 1/de at first r
#else
          double r0de=r1de/std::sqrt(currentParams.julia.alphaShort.getMag_double());
          rfde=r0de*std::sqrt(juliaData.nearmr_fz_.getMag_double()/currentParams.julia.c.dist2_double(&juliaData.nearmr_f, &tmp))/2;
#endif
#if 0 //something old
          rde=r1de
          //rde*=std::sqrt(juliaData.nearmr_fz.getMag_double());  //rde at first_z
          double v=std::sqrt(interior.inte.getMag_double()*juliaData.nearmr_fz_.getMag_double())*rde;
            //v.. relative distance to r1, for small v (<~0.3)
            //    relative distance to pool=1 edge, for large v (~1)
          recipr.assign(&juliaData.nearmr_f);
          recipr.sub(&currentParams.julia.root);
          double u=std::sqrt(recipr.getMag_double())*rde;
            //u..relative distance from nearmr to r - measures how much correction is needed
          //need rde*(1-v) for v~u; rde*v for v<<u; actually 1/v and
          //rde*v*(1+v/u^2*(1-2*u))
          //interior.inte.mul_double(1+v/(u*u)*(1-2*u));
          //another func: rde*v*(u+1-2*v)/(u+1)*(u^2+v)/(u^2)
          //  w=v/u ->    rde*v* (1+u*(1-2*w))/(u+1)*(u+w)/(u)
          //if (u>0.3) u=1;
          interior.inte.mul_double((u+1-2*v)/(u+1)*(u*u+v)/(u*u));
#elif 0 //weighted average of rfde and phi/phi'
          recipr.assign(&juliaData.nearmr_f);
          recipr.sub(&currentParams.julia.root);
          double de=std::sqrt(interior.inte.getMag_double());
          if (1/rfde-de>de)
          {
            double u;
            constexpr double MAGIC1=0.25;
            constexpr double MAGIC2=2.0;
            u=std::sqrt(recipr.getMag_double())*r1de;
            //make it such that the correction becomes 0 halfway between root0 and [0,0]
            //find coef such that (root0/2)^2*k=1    (root0/2)^2*k=4/root0^2
            u=u*u *MAGIC1/MandelMath::sqr_double(3.0/4.0*currentParams.julia.root0.getMag_double()*r1de);
            if (u>1)
              u=1;
            de=(MAGIC2/rfde-1*de)*(1-u)+(de)*(u);
          };
          interior.inte.zero(de, 0);
#else // phi/phi'*(1+phi*rfde)
          //double corr=1+std::sqrt(phiabs*interior.alphak.getMag_double())/rfde; //"rfde"=1/actual rfde
          //now phi~1/(first_z-fr) so skip alphak
          bulb.t1.assign(interior.inte);
          bulb.t1.recip();
          bulb.t1.mul(interior.step_to_root);
          //bulb.t1.chs();
          double weight=bulb.t1.getDist1_tmp()->toDouble();
          if (weight>1)
            weight=1;
          juliaData.store->interior.f1_re=bulb.t1.re.toDouble();
          juliaData.store->interior.f1_im=bulb.t1.im.toDouble();

          //need to find rfde as a complex number
          //for now, just use direction c->r1
          bulb.t2.assign(currentParams.julia.root);
          bulb.t2.sub(currentParams.julia.c);
          bulb.t2.mul(juliaData.nearmr_fz);
          bulb.t2.sign();
          bulb.t2.mul_double(1/rfde);
          bulb.t2.sub(interior.inte);

          bulb.t3.assign(bulb.t2);
          bulb.t3.mul_double(1-weight);
          interior.inte.mul_double(weight);
          interior.inte.add(bulb.t3);

          //double corr=1+std::sqrt(phiabs)/rfde; //"rfde"=1/actual rfde
          //interior.inte.mul_double(corr);
#endif

          //scale from root's pool back to first_z
          /*already included in inte    recipr.assign(&juliaData.nearmr_fz_);
          recipr.recip(&tmp);
          interior.inte.mul(&recipr, &tmp);*/

          //interior.inte.mul_double(1+std::sqrt(phiabs)/currentParams.juliaAlphaShort.getMag_double());
          //interior.inte.mul_double((1+std::sqrt(phiabs)/rde));// /std::sqrt(juliaData.nearmr_fz.getMag_double()));
#else
          //w2.assign(&juliaData.f);
          //w2.ln_approx();
          //interior.inte.mul(&w2, &tmp);
#endif

          interior.inte_abs.assign(*interior.inte.getMag_tmp());
          interior.inte_abs.sqrt();

          juliaData.store->interior.mixed=interior.inte_abs.toDouble();
          juliaData.store->interior.full=juliaData.store->interior.mixed/4;
          juliaData.store->interior.mixed_re=interior.inte.re.toDouble();
          juliaData.store->interior.mixed_im=interior.inte.im.toDouble();
          juliaData.store->interior.phi_re=phi.re.toDouble();
          juliaData.store->interior.phi_im=phi.im.toDouble();
          juliaData.store->interior.phi1_re=phi_z.re.toDouble();
          juliaData.store->interior.phi1_im=phi_z.im.toDouble();

          //bulb.t1.chs();
          juliaData.store->interior.phi2_re=interior.step_to_root.re.toDouble();//the_fzz.re.toDouble();
          juliaData.store->interior.phi2_im=interior.step_to_root.im.toDouble();//the_fzz.im.toDouble();

          /*same bulb.t1.assign(&loope.f);
          bulb.t2.assign(&loope.f_z);
          bulb.t3.assign(&loope.f_zz);
          loope.eval_zz(&tmp,juliaData.store->iter+didcycles+1, &currentParams.c, &currentParams.first_z, false, true);
          bulb.t1.sub(&loope.f);
          bulb.t2.sub(&loope.f_z);
          bulb.t3.sub(&loope.f_zz);*/

          /*juliaData.shrinkfactor.mul(&currentParams.julia.alphaShort, &tmp);
          juliaData.shrinkfactor.recip(&tmp);
          juliaData.store->interior.phi2_re=juliaData.shrinkfactor.re.toDouble();
          juliaData.store->interior.phi2_im=juliaData.shrinkfactor.im.toDouble();*/
          juliaData.store->interior.cycles_until_root=didcycles;
        }
        #if 0 //disable interior completely
        /*
        loope.eval_zz(&tmp, juliaPeriod, &currentParams.c, &currentParams.first_z, false, true); //this is the same for all points...
        interior.inte_abs.assign(*loope.f_z.getMag_tmp(&tmp));
        interior.inte_abs.sqrt();
        interior.fz.assign(&loope.f_z); //later moved to currentData.fz_r
        //interior.inte_abs.chs();
        //interior.inte_abs.add_double(1);
        loope.eval_zz(&tmp, juliaPeriod, &currentParams.c, &currentData.root, false, true); //this is the same for all points...
        interior.fz_mag.assign(*loope.f_zz.getMag_tmp(&tmp));
        interior.fz_mag.recip();
        interior.inte_abs.mul(interior.fz_mag);
        //currentData.store->period=estimateInterior(foundperiod, &currentParams.c, &currentData.root);
        */
#if 1 //estimate distance from full f,f',f'' in 1 step
        //find the root we will converge to
        //loope.eval_zz(&tmp, (currentData.store->nearziter_0+juliaPeriod-1)%juliaPeriod, &currentParams.c, &currentData.root, false, true);
        //loope.eval_zz(&tmp, (currentData.store->near0iter_1)%juliaPeriod, &currentParams.c, &currentData.root, false, true);
        loope.eval_zz(&tmp, juliaPeriod-1-(juliaData.store->near0iter_1+juliaPeriod-1)%juliaPeriod, &currentParams.c, &currentParams.juliaRoot, false, true);
        newt.tmp1.assign(&loope.f);
        interior.inte.assign(&loope.f);
        //alpha=1/f'^period[r0]
        loope.eval_zz(&tmp, juliaPeriod, &currentParams.c, &currentParams.juliaRoot, false, true); //this is the same for all points...
        interior.alpha.assign(&loope.f_z);
        interior.alpha.recip(&tmp);
        interior.alphak.zero(1.0, 0);
        //and let's go
        interior.zoom.zero(1, 0);
        loope.eval_zz(&tmp, 0, &currentParams.c, &currentParams.first_z, false, true);
#else //phase 1: iterate until near0iter; phase 2: as above then scale
      //it has promise but need better near0iter
        //find the root we will converge to - the one nearest to 0
        //sqrt(root-c) should work but choosing the right sign will need some effort
          //but we need the full cycle to find alpha anyway
        loope.eval_zz(&tmp, juliaPeriod-1, &currentParams.c, &currentData.root, false, true);
        newt.tmp1.assign(&loope.f);
        interior.inte.assign(&loope.f);
        //alpha=1/f'^period[r0]
        loope.eval_zz(&tmp, 1, &currentParams.c, &currentData.root, false, false);
        interior.alpha.assign(&loope.f_z);
        interior.alpha.recip(&tmp);
        interior.alphak.zero(1.0, 0);
        //and let's go
        loope.eval_zz(&tmp, currentData.store->near0iter_1-1, &currentParams.c, &currentParams.first_z, false, true);
        interior.zoom.assign(&loope.f_z);
        interior.zoom.recip(&tmp);
        //restart for phase 2
        loope.eval_zz(&tmp, 0, &currentParams.c, &loope.f, false, true);
#endif
        double stop_at=interior.alpha.getMag_double()*loope.f.re.eps234();
        int didcycles=0;
        juliaData.store->interior.f1_re=0;
        juliaData.store->interior.f1_im=0;
        juliaData.store->interior.first_under_1=-1;
        for (int i=0; i<juliaData.store->iter*100/juliaPeriod; i++)
        //for (int i=0; i<5; i++)
        {
          //if ((currentData.store->iter*10-i)%juliaPeriod==0)
          {
            if (interior.inte.dist2_double(&loope.f, &tmp)<stop_at)
              break;
          };
          // //new_de = old_de*(1-d/(z0+d)/2)    z0=root d=f-root
          //new_de = old_de*(1-(z0-r0)/z0/2)
          /*bulb.t1.assign(&loope.f);
          bulb.t1.sub(&currentData.fz_r);
          bulb.t2.assign(&loope.f);
          bulb.t2.recip(&tmp);
          bulb.t2.mul(&bulb.t1, &tmp);
          bulb.t2.lshift(-1);
          bulb.t2.chs();
          bulb.t2.re.add_double(1);
          interior.inte.mul(&bulb.t2, &tmp);
          //advance root
          currentData.fz_r.sqr(&tmp);
          currentData.fz_r.add(&currentParams.c);*/
          //advance f
          //loope.f.sqr(&tmp);
          //loope.f.add(&currentParams.c);
          loope.eval_zz(&tmp, juliaPeriod, &currentParams.c, &currentParams.juliaRoot, false, false);
          if (juliaData.store->interior.first_under_1<0)
          {
            if (loope.f_z.getMag1_tmp(&tmp)->toDouble()<0)
              juliaData.store->interior.first_under_1=i;
          }
          if (i==5-1)
          {
            juliaData.store->interior.f1_re=loope.f_z.re.toDouble();
            juliaData.store->interior.f1_im=loope.f_z.im.toDouble();
          };
          interior.alphak.mul(&interior.alpha, &tmp);
          didcycles++;
          //what would be the effect on f_zz if I multiplied f_z by alpha?
        }
        //currentData.store->interior.f1_re=loope.f_z.re.toDouble();
        //currentData.store->interior.f1_im=loope.f_z.im.toDouble();
        /*
        in summary
          P=b+(f-r)*alpha^k   where |b|=1, direction away from 0   //could use 2*P
          G=-1/P*f'*alpha^k
          Q=alpha^k*(2*f''*P - 3*alpha^k*f'^2)/P^2   QQ=alpha^k*(2*f''*P - 3*alpha^k*f'^2)=2*f''*P*alpha^k - 3*(alpha^k*f')^2
          d.e.=-2/(G+-sqrt(Q)) = 2*P/(f'*alpha^k+-sqrt(QQ))
        */

        interior.inte.rsub(&loope.f);
        interior.inte.mul(&interior.alphak, &tmp); //
        //interior.fz.assign(&interior.inte); //(f-r)*alpha^k for painting
        bulb.t3.assign(&interior.inte);
        bulb.t3.sign(&tmp); //b
        MandelMath::complex<BASE> &the_b=newt.newtX;
        MandelMath::complex<BASE> &the_coef=newt.f_r;
        MandelMath::complex<BASE> &the_fz=newt.fzzf;
        MandelMath::complex<BASE> &the_fzz=newt.fzz_r;
        the_b.assign(&bulb.t3);
        interior.inte.add(&bulb.t3); //P=b+(f-r)*alpha^k
        //interior.fz.assign(&interior.inte); //(f-r)*alpha^k for painting
        //interior.fz.recip(&tmp);

        //f'=-(a b fz)/(a (f-r) + b)^2=-(a b fz)/P^2
        //f''=a b (2 a fz^2 - fzz*P)/P^3=-a*(2 f' fz +b*fzz/P)/P=2 f'^2/b^2*P*b -a*b*fzz/P^2
        bulb.t3.assign(&interior.inte);
        bulb.t3.recip(&tmp);
        the_fz.assign(&bulb.t3);
        the_fz.sqr(&tmp);
        //the_fz.mul(&the_b, &tmp);
        the_fz.chs(); //-1/P^2
        the_fzz.assign(&loope.f_zz);
        the_fzz.mul(&the_fz, &tmp); //-fzz/P^2
        the_fzz.mul(&the_b, &tmp); //-b*fzz/P^2
        the_fzz.mul(&interior.alphak, &tmp); //-a*b*fzz/P^2

        bulb.t3.assign(&loope.f_z);
        bulb.t3.mul(&interior.alphak, &tmp); //alpha^k*f'=a*fz
        the_fz.mul(&bulb.t3, &tmp); //-(a fz)/P^2
        bulb.t1.assign(&the_fz);
        bulb.t1.sqr(&tmp);
        bulb.t1.mul(&the_b, &tmp);
        bulb.t1.mul(&interior.inte, &tmp); //f'^2/b^2*P*b
        bulb.t1.lshift(1);
        the_fzz.add(&bulb.t1); //f''=2 f'^2/b^2*P*b -a*b*fzz/P^2
        the_fz.mul(&the_b, &tmp); //f'=-a*b*fz/P^2
#if 1 //using Laguerre to find lin(phi)=0 on linearized (and inverted) phi,phi',phi''
        bulb.t1.assign(&loope.f_zz);
        bulb.t1.lshift(1);
        bulb.t1.mul(&interior.alphak, &tmp); //
        //interior.fz.assign(&bulb.t1); //2*f''*alpha^k for painting
        bulb.t1.mul(&interior.inte, &tmp); //2*f''*P*alpha^k
        bulb.t2.assign(&bulb.t3);
        bulb.t2.sqr(&tmp);
        bulb.t2.mul_double(3); //3*(alpha^k*f')^2
        bulb.t1.sub(&bulb.t2); //QQ
        bulb.t1.sqrt(&tmp);
        bulb.t2.assign(&bulb.t1);
        if (bulb.t1.mulreT_tmp(&bulb.t3, &tmp)->toDouble()<0)
        {
          bulb.t1.chs();
        }
        else
        {
          bulb.t2.chs();
        }
        bulb.t2.add(&bulb.t3);
        bulb.t2.recip(&tmp);
        bulb.t2.mul(&interior.inte, &tmp);
        bulb.t2.lshift(1);
        bulb.t1.add(&bulb.t3);
        bulb.t1.recip(&tmp);
        interior.inte.mul(&bulb.t1, &tmp);
        interior.inte.lshift(1);
#else
#if 1 //find point where lin(phi')=1
        //d.e.=2/(|fz|+sqrt(|fz|^2+2*|fzz|))
        //looks similar to the full thing but worse

        double discr=the_fz.getMag_double()+2*sqrt(the_fzz.getMag_double());
        discr=sqrt(discr)+sqrt(the_fz.getMag_double());
        discr=2/discr;
        //we have length, point in direction fz/fzz
        interior.inte.assign(&the_fzz);
        interior.inte.recip(&tmp);
        interior.inte.mul(&the_fz, &tmp);
        interior.inte.mul_double(discr/sqrt(interior.inte.getMag_double()));
        interior.inte_abs.zero(discr);

        //find point where f'=1
        //d.e.=2/(|fz|+sqrt(|fz|^2+2*|fzz|))
        //fails when iter is too high because f'~0 f''~0 -> dist ~sqrt(1/f'') which is very large
        /*
        double discr=loope.f_z.getMag_double()+2*sqrt(loope.f_zz.getMag_double());
        discr=sqrt(discr)+sqrt(loope.f_z.getMag_double());
        discr=2/discr;
        //we have length, point in direction f'/f''
        interior.inte.assign(&loope.f_zz);
        interior.inte.recip(&tmp);
        interior.inte.mul(&loope.f_z, &tmp);
        interior.inte.mul_double(discr/sqrt(interior.inte.getMag_double()));
        interior.inte_abs.zero(discr);
        */
#else
#if 0
        //d.e.=fz/fzz
        //seems to range from 0 to infinity
        //interior.inte.assign(&the_fzz);
        //interior.inte.recip(&tmp);
        //interior.inte.mul(&the_fz, &tmp);
#else
        //d.e.=f'/f'' pretty much the same as fi'/fi''
        //interior.inte.assign(&loope.f_zz);
        //interior.inte.recip(&tmp);
        //interior.inte.mul(&loope.f_z, &tmp);
        //d.e.=(fz/|fz|-fz)/fzz=fz*(1/|fz|-1)/fzz
        interior.inte.assign(&loope.f_zz);
        interior.inte.recip(&tmp);
        interior.inte.mul(&loope.f_z, &tmp);
        interior.inte_abs.assign(*loope.f_z.getMag_tmp(&tmp));
        interior.inte_abs.recip();
        interior.inte_abs.sqrt();
        interior.inte_abs.add_double(-1);
        interior.inte.mul(interior.inte_abs);
#endif
#endif
#endif
        //check if the_coef*(x-root1)(x-root2) == fi(x) = 1-1/(1+b/((f-r)*alpha^k))
        bulb.t3.assign(&currentParams.first_z);
        loope.eval_zz(&tmp, juliaPeriod*didcycles, &currentParams.c, &bulb.t3, false, true);
        newt.laguG.assign(&bulb.t3);
        newt.laguG2.assign(&bulb.t3);
        newt.laguG.sub(&bulb.t2);
        newt.laguG2.sub(&interior.inte);
        newt.laguG.mul(&newt.laguG2, &tmp);
        loope.f.sub(&newt.tmp1);
        loope.f.mul(&interior.alphak, &tmp);
        loope.f.recip(&tmp);
        loope.f.mul(&the_b, &tmp);
        loope.f.re.add_double(1);
        loope.f.recip(&tmp);
        loope.f.chs();
        loope.f.re.add_double(1);
        newt.laguG.recip(&tmp);
        newt.laguG.mul(&loope.f, &tmp);
        the_coef.assign(&newt.laguG);
        double f0_re=loope.f.re.toDouble(), f0_im=loope.f.im.toDouble();

//#define DERIVATIVE_H 0.000001
#define DERIVATIVE_H 0.01
        bulb.t3.assign(&currentParams.first_z);
        bulb.t3.re.add_double(DERIVATIVE_H);
        loope.eval_zz(&tmp, juliaPeriod*didcycles, &currentParams.c, &bulb.t3, false, true);
        newt.laguG.assign(&bulb.t3);
        newt.laguG2.assign(&bulb.t3);
        newt.laguG.sub(&bulb.t2);
        newt.laguG2.sub(&interior.inte);
        newt.laguG.mul(&newt.laguG2, &tmp);
        newt.laguG.mul(&the_coef, &tmp);
        loope.f.sub(&newt.tmp1);
        loope.f.mul(&interior.alphak, &tmp);
        loope.f.recip(&tmp);
        loope.f.mul(&the_b, &tmp);
        loope.f.re.add_double(1);
        loope.f.recip(&tmp);
        loope.f.chs();
        loope.f.re.add_double(1);
        double fp_re=loope.f.re.toDouble(), fp_im=loope.f.im.toDouble();
        newt.laguG.sub(&loope.f);

        bulb.t3.assign(&currentParams.first_z);
        //bulb.t3.im.add_double(0.000001);
        bulb.t3.re.add_double(-DERIVATIVE_H);
        loope.eval_zz(&tmp, juliaPeriod*didcycles, &currentParams.c, &bulb.t3, false, true);
        newt.laguH.assign(&bulb.t3);
        newt.laguX.assign(&bulb.t3);
        newt.laguH.sub(&bulb.t2);
        newt.laguX.sub(&interior.inte);
        newt.laguH.mul(&newt.laguX, &tmp);
        newt.laguH.mul(&the_coef, &tmp);
        loope.f.sub(&newt.tmp1);
        loope.f.mul(&interior.alphak, &tmp);
        loope.f.recip(&tmp);
        loope.f.mul(&the_b, &tmp);
        loope.f.re.add_double(1);
        loope.f.recip(&tmp);
        loope.f.chs();
        loope.f.re.add_double(1);
        double fm_re=loope.f.re.toDouble(), fm_im=loope.f.im.toDouble();
        newt.laguH.sub(&loope.f);

        //f'=(fp-f0)/0.000001;
        //f''=((fp-f0)/h-(f0-fm)/h)/h=(fp-2*f0+fm)/h^2
        fm_re=(fp_re-2*f0_re+fm_re)/(DERIVATIVE_H*DERIVATIVE_H);
        fm_im=(fp_im-2*f0_im+fm_im)/(DERIVATIVE_H*DERIVATIVE_H);
        fp_re=(fp_re-f0_re)/DERIVATIVE_H;
        fp_im=(fp_im-f0_im)/DERIVATIVE_H;


        /* what is f'/f'' ? quite similar to f/f'
        interior.inte.assign(&loope.f_zz);
        interior.inte.recip(&tmp);
        interior.inte.mul(&loope.f_z, &tmp);
        */
        interior.inte.mul(&interior.zoom, &tmp);
        interior.inte_abs.assign(*interior.inte.getMag_tmp(&tmp));
        interior.inte_abs.sqrt();

        //for painting
        //(f-r)*alpha^k above
        //interior.fz.assign(&bulb.t3); //f'*alpha^k
        //2*f''*alpha^k above
        //interior.fz.assign(&interior.inte); //complex distance estimate



        //interior.inte_abs.recip();
        //de=1/(f''[r0]/f'[r0]-f'[z0]/(z2-r2))
        /*interior.inte.recip(&tmp);
        loope.eval_zz(&tmp, currentData.store->iter+10, &currentParams.c, &currentData.root, false, true);
        loope.f_z.recip(&tmp);
        loope.f_zz.mul(&loope.f_z, &tmp);
        interior.inte.rsub(&loope.f_zz);
        interior.inte.recip(&tmp);
        interior.inte_abs.assign(*interior.inte.getMag_tmp(&tmp));
        interior.inte_abs.sqrt();*/
        juliaData.store->interior.hits=interior.inte_abs.toDouble();
        juliaData.store->interior.hits_re=interior.inte.re.toDouble();
        juliaData.store->interior.hits_im=interior.inte.im.toDouble();
        juliaData.store->interior.phi_re=f0_re;
        juliaData.store->interior.phi_im=f0_im;
        juliaData.store->interior.phi1_re=the_fz.re.toDouble();
        juliaData.store->interior.phi1_im=the_fz.im.toDouble();
        juliaData.store->interior.phi2_re=the_fzz.re.toDouble();
        juliaData.store->interior.phi2_im=the_fzz.im.toDouble();
        juliaData.store->interior.cycles_until_root=didcycles;
        //currentWorker->assign(&currentData.fc_c_re, &eval.fz_r_re);
        //currentWorker->assign(&currentData.fc_c_im, &eval.fz_r_im);
        #endif //disable interior completely
        break;
      };
    if (need_reset_check==1+100)
    {
      juliaData.store->rstate=JuliaPointStore::ResultState::stPeriod2;
      interior.inte.zero(0.3, 0);
      interior.inte_abs.zero(0.3);
      juliaData.store->interior.mixed=0.3;
      juliaData.store->interior.full=0.25;
      juliaData.store->interior.mixed_re=0.3;
      juliaData.store->interior.mixed_im=0;
      juliaData.store->interior.phi_re=0.3;
      juliaData.store->interior.phi_im=0;
      juliaData.store->interior.phi1_re=0.3;
      juliaData.store->interior.phi1_im=0;
      juliaData.store->interior.phi2_re=0.3;
      juliaData.store->interior.phi2_im=0;
    }
    if (need_reset_check>=1)
    {
      juliaData.store->lookper_startiter=juliaData.store->iter;
      juliaData.lookper_fz.zero(1.0, 0);
      juliaData.lookper_distr.assign(juliaData.f);
      juliaData.lookper_distr.sub(currentParams.julia.root);
      if (currentParams.breakOnNewNearest)
      {
        juliaData.store->iter++;
        break;
      };
    };

  };

  //data.state=MandelPoint::State::stMaxIter;
  /*if (!juliaData.store->has_fc_r && currentParams.want_fc_r &&
      ((juliaData.store->rstate==JuliaPointStore::ResultState::stPeriod2) ||
       (juliaData.store->rstate==JuliaPointStore::ResultState::stPeriod3)))
  {
    loope.eval2(&tmp, juliaData.store->period, &currentParams.c, &currentParams.first_z, false);
    juliaData.fz_z.assign(&loope.f_z);
    juliaData.store->has_fc_r=true;
  };*/
  if (currentParams.want_extangle &&
      juliaData.store->rstate==JuliaPointStore::ResultState::stOutside)
  {
    if (juliaData.store->iter>10000)
      juliaData.extangle.zero(ExtAngle::SPECIAL_VALUE_DEEP);
    else
      extangle.computeMJ(&juliaData.extangle, false, juliaData.store->iter, currentParams.julia.c, currentParams.julia.first_z);
      //extangle.compute3(&currentData.extangle, currentData.store->iter, &currentParams.c, &tmp);
    juliaData.store->rstate=JuliaPointStore::ResultState::stOutAngle;
  };
}


template<typename BASE>
MandelEvaluator<BASE>::ComputeParams::ComputeParams(MandelMath::number<BASE>::Scratchpad *spad):
  mandel(spad), julia(spad),
  epoch(-1), pixelIndex(-1), maxiter(1), breakOnNewNearest(false), want_fc_r(false), want_extangle(false),
  nth_fz(0)
{
}

template<typename BASE>
MandelEvaluator<BASE>::NewtRes::NewtRes(MandelMath::number<BASE>::Scratchpad *spad):
  cyclesNeeded(-1),
  fz_r(spad), nth_fz(spad),
  first_guess_lagu(spad), first_guess_newt(spad)
{
  //working_assert(self_allocator.checkFill());
}

template<typename BASE>
MandelEvaluator<BASE>::Eval::Eval(MandelMath::number<BASE>::Scratchpad *spad):
  fz_r(spad), fz_mag(spad)
{
  //working_assert(self_allocator.checkFill());
}

template<typename BASE>
MandelEvaluator<BASE>::Newt::Newt(MandelMath::number<BASE>::Scratchpad *spad):
  bestr(spad), f_r(spad), fzz_r(spad), tmp1(spad),
  laguH(spad), laguG(spad), laguG2(spad),
  laguX(spad), newtX(spad), prevR(spad), prevGz(spad),
  fzzf(spad), tmp2(spad)
{
  //working_assert(self_allocator.checkFill());
}

template<typename BASE>
MandelEvaluator<BASE>::InteriorInfo::InteriorInfo(MandelMath::number<BASE>::Scratchpad *spad):
  inte(spad), inte_abs(spad), fz(spad), fz_mag(spad), step_to_root(spad), alphak_other(spad), alphak(spad), zoom(spad)
{
  //working_assert(self_allocator.checkFill());
}

template<class BASE>
MandelEvaluator<BASE>::ExtAngle::ExtAngle(MandelMath::number<BASE>::Scratchpad *spad, MandelLoopEvaluator<BASE> *loope, LaguerreStep<BASE> *lagus, MandelEvaluator<BASE> *owner):
  owner(owner), ntype(spad->ntype), loope(loope), lagus(lagus),
  z(spad), r(spad), angleC(spad), angle(spad), x(spad), target(spad),
  dldlz(spad), d2ldlz2(spad), dbg(spad)
{
}

template<>
MandelEvaluator<MandelMath::number_any>::ExtAngle::ExtAngle(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, MandelLoopEvaluator<MandelMath::number_any> *loope, LaguerreStep<MandelMath::number_any> *lagus, MandelEvaluator<MandelMath::number_any> *owner):
  owner(owner), ntype(spad->ntype), loope(loope), lagus(lagus),
  z(spad), r(spad), angleC(spad), angle(spad), x(spad), target(spad),
  dldlz(spad), d2ldlz2(spad), dbg(spad)
{
  switch (ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      dbgPoint();
      goto lolwut;
    case MandelMath::NumberType::typeDouble: lolwut:
      break;
    case MandelMath::NumberType::typeFloat128:
      break;
    case MandelMath::NumberType::typeDDouble:
      break;
    case MandelMath::NumberType::typeQDouble:
      break;
    case MandelMath::NumberType::typeReal642:
      break;
  }
}

template<class BASE>
MandelEvaluator<BASE>::ExtAngle::~ExtAngle()
{
  loope=nullptr;
}

template<>
MandelEvaluator<MandelMath::number_any>::ExtAngle::~ExtAngle()
{
  switch (ntype)
  {
    case MandelMath::NumberType::typeEmpty:
      dbgPoint();
      break;
    case MandelMath::NumberType::typeDouble:
    {
    } break;
    case MandelMath::NumberType::typeFloat128:
    {
    } break;
    case MandelMath::NumberType::typeDDouble:
    {
    } break;
    case MandelMath::NumberType::typeQDouble:
    {
    } break;
    case MandelMath::NumberType::typeReal642:
    {
    } break;
  }
}

#if 0
template<class BASE>
void MandelEvaluator<BASE>::ExtAngle::compute3(MandelMath::number<BASE> *result, int iter, const MandelMath::complex<BASE> *c, typename MandelMath::complex<BASE>::Scratchpad *tmp)
{
  //on input, we should have iter such that |f_c^iter(c)|>10000
  //need to limit step to less than ~ exp(pi) to eliminate skips; exp(pi)~23 so we need radius~sqrt(sqrt(10000))
  iter-=2;

  //use "quick path" only
  //iter=1;

  //mix of both
  //iter-=2;
  //if (iter>61)
  //  iter=61;

  z.assign(c);
  angleC.zero(0);
  angle.zero(0);
  bool veryfirst=true;
  dbg.first_guess_valid=0;
  dbg.last_guess_valid=0;
  constexpr int ANGLE_ZOOM=0; //interferes with SAFE_RADIUS (shrinks it)
  for (; iter>=1+ANGLE_ZOOM; iter--)
  {
    for (int i=0; i<10; i++)
    {
      if (!loope->eval_ext_mand(tmp, &z, iter-1))
      {
        result->zero(SPECIAL_VALUE_EDGE);
        return;
      };
      if (i==0)
      {
        if (veryfirst)
        {
          veryfirst=false;
          /*
z^2+c
  |x    x|      |      |
 ---    ---    ---    ---
  |      |     x|      |x
z (. means ?)
 .|0    0|0    +|+    +|.
+ | 0  . | 0  + | .  + | 0
-----  -----  -----  -----
- | 0  - | .  . | 0  - | 0
 -|.    -|-    0|0    .|0
z, fill ? by making each quadrant uniform
 +|0    0|0    +|+    +|0
+ | 0  0 | 0  + | +  + | 0
-----  -----  -----  -----
- | 0  - | -  0 | 0  - | 0
 -|0    -|-    0|0    -|0
          */
          r_.assign(&loope->f);
          double angleC_bits=0;
          int m=0;
          while (m<60) //60=-log2(eps)
          {
            int way=0;
            x.assign(&r_);
            if (r_.isNegative()) way|=0x01; //.im<0 is not enough to tell -3 from +3
            r_.sqr(tmp);
            if (r_.im.isl0()) way|=0x04;
            target.assign(&r_);
            r_.add(&z);
            if (r_.isNegative()) way|=0x08;
            if (r_.ccw_tmp(&target, tmp)->isl0()) way|=0x02;
            double step=0;
            switch (way)
            {
              case  0: step=0; break; //++++ a
              case  1: step=-2; break; //+++- a
              case  2: step=0; break; //++-+ a
              case  3: step=-2; break; //++-- a
              case  4: step=0; break; //+-++ c,d
              case  5: step=-2; break; //+-+- c,d
              case  6: step=2; break; //+--+ b,e
              case  7: step=0; break; //+--- b,e
              case  8: step=0; break; //-+++ c,d
              case  9: step=-2; break; //-++- c,d
              case 10: step=2; break; //-+-+ b,e
              case 11: step=0; break; //-+-- b,e
              case 12: step=2; break; //--++ a
              case 13: step=0; break; //--+- a
              case 14: step=2; break; //---+ a
              case 15: step=0; break; //---- a
            }
            if (step!=0)
              angleC_bits+=ldexp(step, -m-1);

            /*if (r_.re.isl0())
            { //case 2 and 3
              if (r_.im.isl0())
              { //case 3
                if (!x.im.isl0())
                  angleC_bits+=ldexp(1, -m);
              }
              else
              {
                if (x.im.isl0())
                  angleC_bits-=ldexp(1, -m);
              }
            }
            else
            { //case 1 and 4
              if (x.re.isl0())
              {
                if (x.im.isl0())
                  angleC_bits-=ldexp(1, -m);
                else
                  angleC_bits+=ldexp(1, -m);
              };
            }*/
            m++;
            if (r_.getMag_double()*r_.re.eps2()>1.0)
              break;
          }
          r_.arctan2(&angleC, tmp);
          angleC.lshift(-m);
          angleC.add_pi(angleC_bits);
          target.assign(&loope->f);
        }
        else
        { //when squaring f to get target, are we adding a turn or not? We do if f.re<0
          int way=0;
          if (target.isNegative()) way|=0x08; //negative() needed over im<0 during forward iteration
          r_.assign(&target);
          target.sub(&z);
          if (target.im.isl0()) way|=0x04;
          //r_.recip(tmp);
          //r_.mul(&target, tmp);
          //if (r_.im.isl0()) way|=0x02;
          if (r_.ccw_tmp(&target, tmp)->isl0()) way|=0x02;
          target.sqrt(tmp);
          if (target.mulreT_tmp(&loope->f, tmp)->isl0())
            target.chs();
          if (target.isNegative()) way|=0x01; //.im<0 is not enough to tell -3 from +3
          /*static const int table[16]={
             0, //++++ a
            -2, //+++- a
             0, //++-+ a
            -2, //++-- a
             0, //+-++ c,d
            -2, //+-+- c,d
             2, //+--+ b,e
             0, //+--- b,e
             0, //-+++ c,d
            -2, //-++- c,d
             2, //-+-+ b,e
             0, //-+-- b,e
             2, //--++ a
             0, //--+- a
             2, //---+ a
             0  //---- a
          };
          double step=table[way];*/
          double step=0;
          switch (way)
          {
            case  0: step=0; break; //++++ a
            case  1: step=-2; break; //+++- a
            case  2: step=0; break; //++-+ a
            case  3: step=-2; break; //++-- a
            case  4: step=0; break; //+-++ c,d
            case  5: step=-2; break; //+-+- c,d
            case  6: step=2; break; //+--+ b,e
            case  7: step=0; break; //+--- b,e
            case  8: step=0; break; //-+++ c,d
            case  9: step=-2; break; //-++- c,d
            case 10: step=2; break; //-+-+ b,e
            case 11: step=0; break; //-+-- b,e
            case 12: step=2; break; //--++ a
            case 13: step=0; break; //--+- a
            case 14: step=2; break; //---+ a
            case 15: step=0; break; //---- a
          }
          if (step!=0)
          {
            angle.add_double(step);
          };
          angleC.lshift(-1);
          angle.lshift(-1);
        }
        double scale=10/std::sqrt(target.getMag_double());
        target.mul_double(scale);
      };

      //in log-log space: ln(z)+=(ln(target)-ln(f))/ dlog/dlogz f
      //dlog/dlogz f=z f'/f
      //z*=exp( ln(target/f)/ (z f'/f) ) = exp( ln(target/f)*f/(z f') )
      x.assign(&loope->f);
      x.recip(tmp);
      x.mul(&target, tmp);
      x.ln_approx();
      if (x.is0())
        break;
      {
        static double maxerr=1.8;
        double x_im=x.im.toDouble();
        if (x_im<-maxerr)//1.7)
        {
          maxerr=-x_im+0.001;
        }
        if (x_im>maxerr)//1.7)
        {
          maxerr=x_im+0.001;
        }
      }

      //dlog/dlogx f(x): d/du (log(f(exp(u)))), exp(u)=x ->
      //dlog/dlogx f(x) = x f'(x)/f(x) a.k.a. dldlz
      //d^2log/dlogx^2 f(x): d^2/du^2 log(f(exp(u))), exp(u)=x ->
      //d^2log/dlogx^2 f(x) = x^2 f'(x)^2/f(x)^2 (f''(x)*f(x)/f'(x)^2 - 1) + x f'(x)/f(x)
      //                    = ( (f''*f/f'^2 - 1)*dldlz + 1)*dldlz
      //find root of function with this d^2/dx^2, d/dx, = log(target/loope->f)
      //using Ostrowski iteration
      // 1/sqrt(F'^2/F^2-F''/F)=1/sqrt(dldlz^2/logx^2+d2ldlogz2/logx)=
      // 1/sqrt(A^2+((f_cc*f/f_c^2-1)*A*logx+1)*A) where A=z*fc/f/logx
      // 1/sqrt((((f_cc*f/f_c^2-1)*logx+1)*A+1)*A) where A=z*fc/f/logx
      dldlz.assign(&loope->f);
      dldlz.mul(&x, tmp);
      dldlz.recip(tmp);
      dldlz.mul(&z, tmp); // z/f/logx
      dldlz.mul(&loope->f_c, tmp); // z*f'/f/logx=A
#if 1 //use second derivative
      d2ldlz2.assign(&loope->f_cc);
      d2ldlz2.mul(&loope->f, tmp);
      r_.assign(&loope->f_c);
      r_.sqr(tmp);
      r_.recip(tmp);
      d2ldlz2.mul(&r_, tmp);
      d2ldlz2.re.add_double(-1.0);
#else //newton iteration only
      d2ldlz2.zero(-1, 0);
#endif
      d2ldlz2.mul(&x, tmp);
      d2ldlz2.re.add_double(1.0);
      d2ldlz2.mul(&dldlz, tmp);
      d2ldlz2.re.add_double(1.0);
      d2ldlz2.mul(&dldlz, tmp);    // (((f_cc*f/f_c^2-1)*logx+1)*A+1)*A
      d2ldlz2.sqrt(tmp);
      //the classical problem of choosing the right sqrt
      //officially should be the one closer in direction to f/f', i.e. lnx*f/f'/z = 1/A
      //1/sqrt() close to 1/A    means    sqrt() close to A
      if (d2ldlz2.mulreT_tmp(&dldlz, tmp)->isl0())
        d2ldlz2.chs();
      x.assign(&d2ldlz2);
      x.recip(tmp);

      //bool brake=false;//(x.getMag_double()<x.re.eps234());
      bool enough=x.getMag_double()<0.01; //0.0001 looks good, 0.01 looks good
      x.exp_approx();
      z.mul(&x, tmp); //log(z)+=x ... z*=exp(x)
      if (i==0 && dbg.first_guess_valid<3)
      {
        if (dbg.first_guess_valid==0)
          dbg.first_guess_0.assign(&z);
        else if (dbg.first_guess_valid==1)
          dbg.first_guess_1_.assign(&z);
        else
          dbg.first_guess_2.assign(&z);
        dbg.first_guess_valid++;
      };
      if (enough)
        break;
    }
    if (dbg.last_guess_valid<3)
    {
      if (dbg.last_guess_valid==0)
        dbg.last_guess_0.assign(&z);
      else if (dbg.last_guess_valid==1)
        dbg.last_guess_1_.assign(&z);
      else
        dbg.last_guess_2.assign(&z);
      dbg.last_guess_valid++;
    };
    if (owner->workIfEpoch!=owner->busyEpoch)
    {
      result->zero(0);
      return;
    };
  }
  if (veryfirst)
  {
    /*veryfirst=false;
    c->arctan2(&angleC, tmp);
    angleC.lshift(ANGLE_ZOOM);*/
    if (!loope->eval_ext_mand(tmp, c, ANGLE_ZOOM))
    {
      c->arctan2(&angleC, tmp);
      angleC.lshift(ANGLE_ZOOM);
    }
    else
      loope->f.arctan2(&angleC, tmp);
    //angleC.zero(SPECIAL_VALUE_DEEP);
  };
  result->assign(angleC);
  result->add_pi(angle.toDouble());
  //result->add(angle); //*2*pi
}
#endif

template<class BASE>
void MandelEvaluator<BASE>::ExtAngle::computeMJ(MandelMath::number<BASE> *result, bool mandel, int iter, MandelMath::complex<BASE> const &c, MandelMath::complex<BASE> const &zz)
{
  //on input, we should have iter such that |f_c^iter(c)|>10000
  //need to limit step to less than ~ exp(pi) to eliminate skips; exp(pi)~23 so we need radius~sqrt(sqrt(10000))
  iter-=2;

  //use "quick path" only
  //iter=1;

  //mix of both
  //iter-=2;
  //if (iter>61)
  //  iter=61;

  z.assign(zz);
  //life would be easier with: if (mandel) c=&z;
  angleC.zero(0);
  angle.zero(0);
  bool veryfirst=true;
  dbg.first_guess_valid=0;
  dbg.last_guess_valid=0;
  constexpr int ANGLE_ZOOM=0; //interferes with SAFE_RADIUS (shrinks it)
  for (; iter>=1+ANGLE_ZOOM; iter--)
  {
    for (int i=0; i<10; i++)
    {
      if (!loope->eval_ext_mandMJ(mandel?1.0:0.0, mandel?z:c, z, iter-1))
      {
        result->zero(SPECIAL_VALUE_EDGE);
        return;
      };
      if (i==0)
      {
        if (veryfirst)
        {
          veryfirst=false;
          /*
z^2+c
  |x    x|      |      |
 ---    ---    ---    ---
  |      |     x|      |x
z (. means ?)
 .|0    0|0    +|+    +|.
+ | 0  . | 0  + | .  + | 0
-----  -----  -----  -----
- | 0  - | .  . | 0  - | 0
 -|.    -|-    0|0    .|0
z, fill ? by making each quadrant uniform
 +|0    0|0    +|+    +|0
+ | 0  0 | 0  + | +  + | 0
-----  -----  -----  -----
- | 0  - | -  0 | 0  - | 0
 -|0    -|-    0|0    -|0
          */
          r.assign(loope->f);
          double angleC_bits=0;
          int m=0;
          while (m<60) //60=-log2(eps)
          {
            int way=0;
            x.assign(r);
            if (r.isl0()) way|=0x01; //.im<0 is not enough to tell -3 from +3
            r.sqr();
            if (r.isl0()) way|=0x04; //.im.isl0 creates thin strips in Julia with c=0
            target.assign(r);
            r.add(mandel?z:c);
            if (r.isl0()) way|=0x08;
            if (r.ccw_tmp(target)->isl0()) way|=0x02;
            double step=0;
            switch (way)
            {
              case  0: step=0; break; //++++ a
              case  1: step=-2; break; //+++- a
              case  2: step=0; break; //++-+ a
              case  3: step=-2; break; //++-- a
              case  4: step=0; break; //+-++ c,d
              case  5: step=-2; break; //+-+- c,d
              case  6: step=2; break; //+--+ b,e
              case  7: step=0; break; //+--- b,e
              case  8: step=0; break; //-+++ c,d
              case  9: step=-2; break; //-++- c,d
              case 10: step=2; break; //-+-+ b,e
              case 11: step=0; break; //-+-- b,e
              case 12: step=2; break; //--++ a
              case 13: step=0; break; //--+- a
              case 14: step=2; break; //---+ a
              case 15: step=0; break; //---- a
            }
            if (step!=0)
              angleC_bits+=ldexp(step, -m-1);

            /*if (r_.re.isl0())
            { //case 2 and 3
              if (r_.im.isl0())
              { //case 3
                if (!x.im.isl0())
                  angleC_bits+=ldexp(1, -m);
              }
              else
              {
                if (x.im.isl0())
                  angleC_bits-=ldexp(1, -m);
              }
            }
            else
            { //case 1 and 4
              if (x.re.isl0())
              {
                if (x.im.isl0())
                  angleC_bits-=ldexp(1, -m);
                else
                  angleC_bits+=ldexp(1, -m);
              };
            }*/
            m++;
            if (r.getMag_double()*r.re.eps2()>1.0)
              break;
          }
          r.arctan2(&angleC);
          angleC.lshift(-m);
          angleC.add_pi(angleC_bits);
          target.assign(loope->f);
        }
        else
        { //when squaring f to get target, are we adding a turn or not? We do if f.re<0
          int way=0;
          if (target.isl0()) way|=0x08; //negative() needed over im<0 during forward iteration
          r.assign(target);
          target.sub(mandel?z:c);
          if (target.isl0()) way|=0x04; //.im.isl0 creates thin strips in Julia with c=0
          //r_.recip(tmp);
          //r_.mul(&target, tmp);
          //if (r_.im.isl0()) way|=0x02;
          if (r.ccw_tmp(target)->isl0()) way|=0x02;
          target.sqrt();
          if (target.mulreT_tmp(loope->f)->isl0())
            target.chs();
          if (target.isl0()) way|=0x01; //.im<0 is not enough to tell -3 from +3
          /*static const int table[16]={
             0, //++++ a
            -2, //+++- a
             0, //++-+ a
            -2, //++-- a
             0, //+-++ c,d
            -2, //+-+- c,d
             2, //+--+ b,e
             0, //+--- b,e
             0, //-+++ c,d
            -2, //-++- c,d
             2, //-+-+ b,e
             0, //-+-- b,e
             2, //--++ a
             0, //--+- a
             2, //---+ a
             0  //---- a
          };
          double step=table[way];*/
          double step=0;
          switch (way)
          {
            case  0: step=0; break; //++++ a
            case  1: step=-2; break; //+++- a
            case  2: step=0; break; //++-+ a
            case  3: step=-2; break; //++-- a
            case  4: step=0; break; //+-++ c,d
            case  5: step=-2; break; //+-+- c,d
            case  6: step=2; break; //+--+ b,e
            case  7: step=0; break; //+--- b,e
            case  8: step=0; break; //-+++ c,d
            case  9: step=-2; break; //-++- c,d
            case 10: step=2; break; //-+-+ b,e
            case 11: step=0; break; //-+-- b,e
            case 12: step=2; break; //--++ a
            case 13: step=0; break; //--+- a
            case 14: step=2; break; //---+ a
            case 15: step=0; break; //---- a
          }
          if (step!=0)
          {
            angle.add_double(step);
          };
          angleC.lshift(-1);
          angle.lshift(-1);
        }
        double scale=10/std::sqrt(target.getMag_double());
        target.mul_double(scale);
      };

      //in log-log space: ln(z)+=(ln(target)-ln(f))/ dlog/dlogz f
      //dlog/dlogz f=z f'/f
      //z*=exp( ln(target/f)/ (z f'/f) ) = exp( ln(target/f)*f/(z f') )
      x.assign(loope->f);
      x.recip();
      x.mul(target);
      x.ln_approx();
      if (x.is0())
        break;
      {
        static double maxerr=1.8;
        double x_im=x.im.toDouble();
        if (x_im<-maxerr)//1.7)
        {
          maxerr=-x_im+0.001;
        }
        if (x_im>maxerr)//1.7)
        {
          maxerr=x_im+0.001;
        }
      }

      //dlog/dlogx f(x): d/du (log(f(exp(u)))), exp(u)=x ->
      //dlog/dlogx f(x) = x f'(x)/f(x) a.k.a. dldlz
      //d^2log/dlogx^2 f(x): d^2/du^2 log(f(exp(u))), exp(u)=x ->
      //d^2log/dlogx^2 f(x) = x^2 f'(x)^2/f(x)^2 (f''(x)*f(x)/f'(x)^2 - 1) + x f'(x)/f(x)
      //                    = ( (f''*f/f'^2 - 1)*dldlz + 1)*dldlz
      //find root of function with this d^2/dx^2, d/dx, = log(target/loope->f)
      //using Ostrowski iteration
      // 1/sqrt(F'^2/F^2-F''/F)=1/sqrt(dldlz^2/logx^2+d2ldlogz2/logx)=
      // 1/sqrt(A^2+((f_cc*f/f_c^2-1)*A*logx+1)*A) where A=z*fc/f/logx
      // 1/sqrt((((f_cc*f/f_c^2-1)*logx+1)*A+1)*A) where A=z*fc/f/logx
      dldlz.assign(loope->f);
      dldlz.mul(x);
      dldlz.recip();
      dldlz.mul(z); // z/f/logx
      dldlz.mul(loope->f_c); // z*f'/f/logx=A
#if 1 //use second derivative
      d2ldlz2.assign(loope->f_cc);
      d2ldlz2.mul(loope->f);
      r.assign(loope->f_c);
      r.sqr();
      r.recip();
      d2ldlz2.mul(r);
      d2ldlz2.re.add_double(-1.0);
#else //newton iteration only
      d2ldlz2.zero(-1, 0);
#endif
      d2ldlz2.mul(x);
      d2ldlz2.re.add_double(1.0);
      d2ldlz2.mul(dldlz);
      d2ldlz2.re.add_double(1.0);
      d2ldlz2.mul(dldlz);    // (((f_cc*f/f_c^2-1)*logx+1)*A+1)*A
      d2ldlz2.sqrt();
      //the classical problem of choosing the right sqrt
      //officially should be the one closer in direction to f/f', i.e. lnx*f/f'/z = 1/A
      //1/sqrt() close to 1/A    means    sqrt() close to A
      if (d2ldlz2.mulreT_tmp(dldlz)->isl0())
        d2ldlz2.chs();
      x.assign(d2ldlz2);
      x.recip();

      //bool brake=false;//(x.getMag_double()<x.re.eps234());
      bool enough=x.getMag_double()<0.01; //0.0001 looks good, 0.01 looks good
      x.exp_approx();
      z.mul(x); //log(z)+=x ... z*=exp(x)
      if (i==0 && dbg.first_guess_valid<3)
      {
        if (dbg.first_guess_valid==0)
          dbg.first_guess_0.assign(z);
        else if (dbg.first_guess_valid==1)
          dbg.first_guess_1_.assign(z);
        else
          dbg.first_guess_2.assign(z);
        dbg.first_guess_valid++;
      };
      if (enough)
        break;
    }
    if (dbg.last_guess_valid<3)
    {
      if (dbg.last_guess_valid==0)
        dbg.last_guess_0.assign(z);
      else if (dbg.last_guess_valid==1)
        dbg.last_guess_1_.assign(z);
      else
        dbg.last_guess_2.assign(z);
      dbg.last_guess_valid++;
    };
    if (owner->workIfEpoch!=owner->busyEpoch)
    {
      result->zero(0);
      return;
    };
  }
  if (veryfirst)
  {
    /*veryfirst=false;
    c->arctan2(&angleC, tmp);
    angleC.lshift(ANGLE_ZOOM);*/
    if (!loope->eval_ext_mandMJ(mandel?1.0:0.0, mandel?z:c, z, ANGLE_ZOOM))
    {
      z.arctan2(&angleC);
      angleC.lshift(ANGLE_ZOOM);
    }
    else
      loope->f.arctan2(&angleC);
    //angleC.zero(SPECIAL_VALUE_DEEP);
  };
  result->assign(angleC);
  result->add_pi(angle.toDouble());
  //result->add(angle); //*2*pi
}

template struct LaguerrePoint<double>;
template struct LaguerrePoint<MandelMath::number_any>;
template struct MandelPoint<double>;
template struct MandelPoint<MandelMath::number_any>;
template struct JuliaPoint<double>;
template struct JuliaPoint<MandelMath::number_any>;
template class LaguerreStep<double>;
template class MandelLoopEvaluator<double>;
template class MandelLoopEvaluator<MandelMath::number_any>;
template class MandelEvaluator<double>;
template class MandelEvaluator<MandelMath::number_any>;
template void ComputeMandelParams<double>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<double>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
#if !NUMBER_DOUBLE_ONLY
template struct LaguerrePoint<__float128>;
template struct MandelPoint<__float128>;
template struct JuliaPoint<__float128>;
template class MandelEvaluator<__float128>;
template void ComputeMandelParams<__float128>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<__float128>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
template struct LaguerrePoint<MandelMath::dd_real>;
template struct MandelPoint<MandelMath::dd_real>;
template struct JuliaPoint<MandelMath::dd_real>;
template class MandelEvaluator<MandelMath::dd_real>;
template void ComputeMandelParams<MandelMath::dd_real>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<MandelMath::dd_real>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
template struct LaguerrePoint<MandelMath::dq_real>;
template struct MandelPoint<MandelMath::dq_real>;
template struct JuliaPoint<MandelMath::dq_real>;
template class MandelEvaluator<MandelMath::dq_real>;
template void ComputeMandelParams<MandelMath::dq_real>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<MandelMath::dq_real>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
template struct LaguerrePoint<MandelMath::real642>;
template struct MandelPoint<MandelMath::real642>;
template struct JuliaPoint<MandelMath::real642>;
template class MandelEvaluator<MandelMath::real642>;
template void ComputeMandelParams<MandelMath::real642>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<MandelMath::real642>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
#endif
template void ComputeMandelParams<MandelMath::number_any>::assign_across<MandelMath::number_any>(const ComputeMandelParams<MandelMath::number_any> &src);
template void ComputeJuliaParams<MandelMath::number_any>::assign_across<MandelMath::number_any>(const ComputeJuliaParams<MandelMath::number_any> &src);
//template class MandelEvaluator<MandelMath::worker_multi_double>;
