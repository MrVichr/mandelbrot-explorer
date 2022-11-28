#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "MandelModel.hpp"
#include "MandelEvaluator.hpp"

#define CURRENT_STORE_DIRECT 0 //either both or none of store and numbers must be direct, can't mix
#define UPDATE_CACHED_MOD 0 //works until precision doesn't allow

MandelModel::MandelModel(): QObject(),
  _extAngleZoom(0), pointStore(nullptr), precisionRecord(nullptr)
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);

   /*__float128 p=1.4;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;
   p=(2/p+p)/2;*/

  //selftest of dd_real.split / x87 rounding
  MandelMath::dd_real dd1, dd2;
  /*dd2.hi=0x1fffffff;
  ((uint64_t &)dd2.hi)&=0xfffffffffc000000ull;
  dd2.lo=0x1fffffff-dd2.hi;
  double dd2hi=0x1fffffff;
  ((uint64_t &)dd2hi)&=0xfffffffffc000000ull;
  dd2hi+=3;
  dd1.split(dd2hi);*/
  dd1.split(0x1fffffff); //-> 0x20000000 - 0x1
  dd1.split(0x20000000); //-> 0x20000000 + 0x0
  dd1.split(0x20000001); //-> 0x20000000 + 0x1
  dd1.split(0x20000002); //-> 0x20000000 + 0x2
  dd1.split(0x20000003); //-> 0x20000000 + 0x3
  dd1.split(0x20000004); //-> 0x20000000 + 0x4
  dd1.split(0x20000005); //-> 0x20000000 + 0x5
  dd1.split(0x20000006); //-> 0x20000000 + 0x6
  dd1.split(0x20000007); //-> 0x20000000 + 0x7
  dd1.split(0x20000008); //-> 0x20000000 + 0x8                       100 0000 0000 002                                     80 0000 0000 0010
  dd1.split(0x20000009); //-> 0x20000010 - 0x7  0x20000009*0x8000001=100 0000 6800 0009 @52b= 100 0000 0000 0000 -20000009=FF FFFF DFFF FFF7 rounds to FF FFFF DFFF FFF0
  dd1.split(0x2000000a); //-> 0x20000010 - 0x6
  dd1.split(0x2000000b); //-> 0x20000010 - 0x5
  dd1.split(0x2000000c); //-> 0x20000010 - 0x4
  dd1.split(0x2000000d); //-> 0x20000010 - 0x3
  dd1.split(0x2000000e); //-> 0x20000010 - 0x2
  dd1.split(0x2000000f); //-> 0x20000010 - 0x1
  dd1.split(0x20000010); //-> 0x20000010 + 0x0
  dd1.split(0x20000011); //-> 0x20000010 + 0x1
  dd1.split(-0.00043541188577694845);

  _selectedPaintStyle=paintStyleCls;//Kind;
  _selectedPrecision=precisionDouble;
  epoch=1;
  imageWidth=0;
  imageHeight=0;
  //pointStore=nullptr;
  nextGivenPointIndex=0;
  effortBonus=0;
  //threadCount=4;
  _threadsWorking=0;
  QObject::connect(this, &MandelModel::selectedPrecisionChange,
                   this, &MandelModel::selectedPrecisionChanged);
  selectedPrecisionChanged();
}

MandelModel::~MandelModel()
{
  epoch=(epoch%2000000000)+1;
  delete precisionRecord;
  //delete storeAllocator;
  //delete storeWorker;
  //storeWorker=nullptr;

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

void MandelModel::startRunning()
{
  for (int t=0; t<precisionRecord->threadCount; t++)
  {
    precisionRecord->threads[t]->startRunning();
  }
}

QString MandelModel::pixelXtoRE_str(int x)
{
  MandelMath::number<MandelMath::number_a *> num(precisionRecord->ntype);
  num.assign(precisionRecord->position.center.re);
  num.add_double((x - imageWidth/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString MandelModel::pixelYtoIM_str(int y)
{
  MandelMath::number<MandelMath::number_a *> num(precisionRecord->ntype);
  num.assign(precisionRecord->position.center.im);
  num.add_double((y - imageHeight/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString MandelModel::getTimes()
{
  QString result;
  /*
  for (int t=0; t<precisionRecord->threadCount; t++)
    result+=QString("%1-%2[%3,%4],").
        arg((precisionRecord->threads[t]->timeOuterTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInnerTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInvokePostTotal)/1000000000.0, 0, 'f', 3).
        arg((precisionRecord->threads[t]->timeInvokeSwitchTotal)/1000000000.0, 0, 'f', 3);
  */
  qint64 outer=0, inner=0, invokepost=0, invokeswitch=0, threaded=0;
  for (int t=0; t<precisionRecord->threadCount; t++)
  {
    outer+=precisionRecord->threads[t]->timeOuterTotal;
    inner+=precisionRecord->threads[t]->timeInnerTotal;
    invokepost+=precisionRecord->threads[t]->timeInvokePostTotal;
    invokeswitch+=precisionRecord->threads[t]->timeInvokeSwitchTotal;
    threaded+=precisionRecord->threads[t]->timeThreadedTotal;
  }
  result=QString("%1-%2,%3-%4 =%5,").
      arg((inner)/1000000000.0, 0, 'f', 3).
      arg((invokeswitch)/1000000000.0, 0, 'f', 3).
      arg((invokepost)/1000000000.0, 0, 'f', 3).
      arg((outer)/1000000000.0, 0, 'f', 3).
      arg((threaded)/1000000000.0, 0, 'f', 3);
  return result;
}

QString MandelModel::getTextXY()
{
  if (precisionRecord==nullptr)
    return "-";
  return precisionRecord->orbit.evaluator.currentParams.first_z.re.toString()+" +i* "+
         precisionRecord->orbit.evaluator.currentParams.first_z.im.toString();
}

static QString mandDoubleToString(double x)
{
  if (x>0)
  {
    if (x>10000)
      return QString("+%1").arg(x, 10, 'g');
    else
      return QString("+%1").arg(x, 10, 'f');
  }
  else
  {
    if (x<-10000)
      return QString("%1").arg(x, 10, 'g');
    else
      return QString("%1").arg(x, 10, 'f');//returns like 50 digits QString::number(x, 10, 'f');
  }
}

QString MandelModel::getTextInfoGen()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_a *> tmp(precisionRecord->ntype);
    reimToPixel(&orbit_x, &orbit_y, &precisionRecord->orbit.evaluator.currentParams.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  MandelPointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;
  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*MandelPoint<MandelMath::number_a *>::LEN);
  MandelPoint<MandelMath::number_a *> *data=&precisionRecord->wtiPoint;

  QString state;
  switch (data_store->rstate)
  {
    case MandelPointStore::ResultState::stUnknown:
      if (data_store->wstate.load()==MandelPointStore::WorkState::stIdle)
        state="Unk";
      else
        state="Working...";
      break;
    case MandelPointStore::ResultState::stOutside:
      state="Out"; break;
    case MandelPointStore::ResultState::stOutAngle:
      state="OutA"; break;
    case MandelPointStore::ResultState::stBoundary:
      state="Bound"; break;
    case MandelPointStore::ResultState::stDiverge:
      state="Diver"; break;
    case MandelPointStore::ResultState::stMisiur:
      state="Misiur"; break;
    case MandelPointStore::ResultState::stPeriod2:
      state="Per2"; break;
    case MandelPointStore::ResultState::stPeriod3:
      state="Per3"; break;
    case MandelPointStore::ResultState::stMaxIter:
      state="Max"; break;
  }

  return state+" iter="+QString::number(data_store->iter)+" nearz="+QString::number(data_store->nearziter_0)+" near0="+QString::number(data_store->near0iter_1)+
      " fc="+mandDoubleToString(data->fc_c.re.toDouble())+
             mandDoubleToString(data->fc_c.im.toDouble())+"i";
}

QString MandelModel::getTextInfoSpec()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_a *> tmp(precisionRecord->ntype);
    reimToPixel(&orbit_x, &orbit_y, &precisionRecord->orbit.evaluator.currentParams.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  MandelPointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;
  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*MandelPoint<MandelMath::number_a *>::LEN);
  //MandelPoint<MandelMath::worker_multi> *data=&precisionRecord->wtiPoint;

  switch (data_store->rstate)
  {
    case MandelPointStore::ResultState::stUnknown:
      if (data_store->wstate.load()==MandelPointStore::WorkState::stIdle)
        return " ";
      else
        return "Working...";
      break;
    case MandelPointStore::ResultState::stOutside:
      return QString("ext=")+QString::number(data_store->exterior_hits)+" sure="+QString::number(data_store->surehand); break;
    case MandelPointStore::ResultState::stOutAngle:
      return QString("ext=")+QString::number(data_store->exterior_hits)+
             " phi="+precisionRecord->wtiPoint.extangle.toString()+
             " sure="+QString::number(data_store->surehand); break;
    case MandelPointStore::ResultState::stBoundary:
      return "sure="+QString::number(data_store->surehand);
    case MandelPointStore::ResultState::stDiverge:
      return "sure="+QString::number(data_store->surehand);
    case MandelPointStore::ResultState::stMisiur:
      return "sure="+QString::number(data_store->surehand);
    case MandelPointStore::ResultState::stPeriod2:
    {
      double p=std::atan2((precisionRecord->orbit.bulb.baseFz_.im.toDouble()),
                           precisionRecord->orbit.bulb.baseFz_.re.toDouble())*precisionRecord->orbit.bulb.foundMult_;
      return QString("per=")+QString::number(data_store->period)+" sure="+QString::number(data_store->surehand)+" int="+QString::number(data_store->interior)   +
          " mult="+QString::number(std::round(p/(2*M_PI)))+"/"+QString::number(precisionRecord->orbit.bulb.foundMult_);
    } break;
    case MandelPointStore::ResultState::stPeriod3:
    {
      double p=std::atan2((precisionRecord->orbit.bulb.baseFz_.im.toDouble()),
                           precisionRecord->orbit.bulb.baseFz_.re.toDouble())*precisionRecord->orbit.bulb.foundMult_;
      return QString("per=")+QString::number(data_store->period)+" sure="+QString::number(data_store->surehand)+" int="+QString::number(data_store->interior)   +
          " mult="+QString::number(std::round(p/(2*M_PI)))+"/"+QString::number(precisionRecord->orbit.bulb.foundMult_);
    } break;
    case MandelPointStore::ResultState::stMaxIter:
      return "sure="+QString::number(data_store->surehand);
  }
  return "-?-?-";
}

ShareableViewInfo MandelModel::getViewInfo()
{
  ShareableViewInfo result(precisionRecord->ntype);
  MandelMath::complex<MandelMath::number_a *>::Scratchpad spad(precisionRecord->ntype);
  //result.worker=orbit.worker;
  result.period=precisionRecord->orbit.evaluator.currentData.store->nearziter_0;//evaluator.currentData.lookper_lastGuess;//orbit.pointData.period;
  if (result.period<1)
    result.period=1;
  result.nth_fz=result.period;/*precisionRecord->orbit.evaluator.currentData.store->period;
  if (result.nth_fz<1)
    result.nth_fz=result.period;*/
  result.scale=precisionRecord->position.step_size;
  if ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod2) ||
      (precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod3))
    result.juliaPeriod=precisionRecord->orbit.evaluator.currentData.store->period;
  else
    result.juliaPeriod=1<<MAX_EFFORT;
  /*orbit.worker->init_(&result.re_, &result.re_p);
  orbit.worker->init_(&result.im, &result.im_p);
  orbit.worker->init_(&result.root_re, &result.rre_p);
  orbit.worker->init_(&result.root_im, &result.rim_p);*/
  result.c.assign_across(&precisionRecord->orbit.evaluator.currentParams.c);
  result.root.assign_across(&precisionRecord->orbit.evaluator.currentData.root);
  precisionRecord->orbit.evaluator.loope.eval_zz(&precisionRecord->orbit.evaluator.tmp, result.nth_fz,
      &precisionRecord->orbit.evaluator.currentParams.c,//result.c,
      &precisionRecord->orbit.evaluator.currentData.root,//result.root,
      false);
  result.nth_fz_limit.assign_across(precisionRecord->orbit.evaluator.loope.f_z.getMag_tmp(&spad));

  //TODO: why here? should be somewhere else
  precisionRecord->lagu_c.assign_across(&precisionRecord->orbit.evaluator.currentParams.c);
  precisionRecord->lagu_r.assign_across(&precisionRecord->orbit.evaluator.currentData.root);

  return result;
}

ShareableViewInfo MandelModel::makeViewInfo(const QVariantMap &params)
{
  ShareableViewInfo result(precisionRecord->ntype);
  result.scale=params.value("viewZoom", 0.0).toDouble();
  result.c.zero(params.value("cRe", 0.0).toDouble(), params.value("cIm", 0.0).toDouble());
  //have to be correct result.root.zero(0, 0);
  result.juliaPeriod=params.value("period", 1<<MAX_EFFORT).toInt();
  result.nth_fz=result.juliaPeriod;
  precisionRecord->lagu_c.zero(params.value("cRe", 0.0).toDouble(), params.value("cIm", 0.0).toDouble());
  precisionRecord->lagu_r.zero(0,0);

  //compute the root
  precisionRecord->orbit.evaluator.currentParams.first_z.assign(&result.c);
  precisionRecord->orbit.evaluator.currentParams.epoch=epoch;
  precisionRecord->orbit.evaluator.workIfEpoch=precisionRecord->orbit.evaluator.busyEpoch;//epoch;
  precisionRecord->orbit.evaluator.currentParams.pixelIndex=0;
  precisionRecord->orbit.evaluator.currentParams.c.assign(&precisionRecord->orbit.evaluator.currentParams.first_z);
  //already 0 precisionRecord->orbit.evaluator.currentParams.nth_fz=0;
  //precisionRecord->orbit.evaluator.currentData.store->rstate=MandelPointStore::ResultState::stUnknown_;
  //precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stIdle;
  precisionRecord->orbit.evaluator.currentData.zero(&precisionRecord->orbit.evaluator.currentParams.first_z);
  precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stWorking;
  precisionRecord->orbit.evaluator.currentParams.breakOnNewNearest=false;
  precisionRecord->orbit.evaluator.currentParams.maxiter=1<<MAX_EFFORT;
  precisionRecord->orbit.evaluator.currentParams.want_extangle=false;//don't need for root
  {
    while ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stUnknown) &&
           (precisionRecord->orbit.evaluator.currentData.store->iter<(1<<MAX_EFFORT)))
    {
      //precisionRecord->orbit.evaluator.currentParams.maxiter=1<<MAX_EFFORT; //dont't know->run fully
      precisionRecord->orbit.evaluator.thread.syncMandel();
    }
  }
  if (precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod2 ||
      precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod3)
    result.root.assign(&precisionRecord->orbit.evaluator.currentData.root);
  else
    result.root.zero(0, 0);

  return result;
}

class PixelPositionTransformer
{
public:
  PixelPositionTransformer(int lshift, int new_step_log);
  void setShift(MandelMath::number<MandelMath::number_a *> *shift, int maxOut);
  int transform(int new_value, int offset);
protected:
  bool invalid;
  int new_step_log;
  int step_scale_n_shift, step_scale_d_shift, step_scale_d_mask;
  int delta_int;
};

PixelPositionTransformer::PixelPositionTransformer(int lshift, int new_step_log):
  new_step_log(new_step_log)
{
  if (lshift>=0)
  {
    step_scale_n_shift=0;
    step_scale_d_shift=lshift;
    step_scale_d_mask=(1<<lshift)-1;
  }
  else
  {
    step_scale_n_shift=-lshift;
    step_scale_d_shift=0;
    step_scale_d_mask=0;
  }
  invalid=lshift>10 || lshift<-10;
  delta_int=0;
}

void PixelPositionTransformer::setShift(MandelMath::number<MandelMath::number_a *> *shift, int maxOut)
{
  shift->lshift(new_step_log+step_scale_n_shift); //new_step_log+step_scale_n_shift = max(old_step_log, new_step_log)
  double test=shift->toDouble();
  if (test>1000000 || test<-1000000) //somewhat related to MAX_INT/2^max_valid_shift
    invalid=true;
  delta_int=shift->toRound();//old_worker->toRound(tmp.ptr);
  delta_int-=(maxOut/2)<<step_scale_n_shift;
}

int PixelPositionTransformer::transform(int new_value, int offset)
{
  if (invalid)
    return -1;
  else if ((step_scale_d_mask==0) || ((((new_value<<step_scale_n_shift)+delta_int)&step_scale_d_mask)==0))
    return (offset) + (((new_value<<step_scale_n_shift)+delta_int)>>step_scale_d_shift);
  else
    return -1;//call reset() in the second loop, we may still need the points  =imageHeight;
}

//template <class WORKER_MULTI>
void MandelModel::transformStore(void *old_points, MandelPointStore *old_store, int old_width, int old_height, const MandelMath::complex<MandelMath::number_a *> *old_c,
                                 void *new_points, MandelPointStore *new_store, int new_width, int new_height, const MandelMath::complex<MandelMath::number_a *> *new_c,
                                 int inlog, int new_step_log)
{
  if (precisionRecord==nullptr)
  {
    dbgPoint();
    return;
  };
  /*int indexOfWtiPoint, wtiIndexLen;
  precisionRecord->wtiPoint.self_allocator._getFirstCapac(indexOfWtiPoint, wtiIndexLen);
  working_assert((wtiIndexLen==MandelPoint<WORKER_MULTI, 0>::LEN));*/
  PixelPositionTransformer ytrans=PixelPositionTransformer(inlog, new_step_log);
  PixelPositionTransformer xtrans=PixelPositionTransformer(inlog, new_step_log);
  {
    MandelMath::number<MandelMath::number_a *> tmp(old_c->im.ntype());//precisionRecord->ntype);
    tmp.assign(old_c->im);
    tmp.sub(new_c->im); //and reversing y at the last minute
    ytrans.setShift(&tmp, new_height);

    tmp.assign(new_c->re);
    tmp.sub(old_c->re);
    xtrans.setShift(&tmp, new_width);
  }

  MandelMath::complex<MandelMath::number_a *> first_z(precisionRecord->ntype);
  for (int newy=0; newy<new_height; newy++)
  {
    int oldy=ytrans.transform(newy, old_height/2);
    //-1 is better than imageHeight, call reset() in the second loop, we may still need the points
    first_z.im.assign(precisionRecord->position.center.im);
    first_z.im.add_double((new_height/2-newy)*precisionRecord->position.step_size);
    for (int newx=0; newx<new_width; newx++)
    {
      int oldx=xtrans.transform(newx, old_width/2);
      //-1 is better than imageWidth, call reset() in the second loop, we may still need the points
      if ((oldy>newy) || ((oldy==newy) && (oldx>newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old point to new place
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].wstate==MandelPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=MandelPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*MandelPoint<MandelMath::number_a *>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*MandelPoint<MandelMath::number_a *>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*MandelPoint<WORKER_MULTI, 0>::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint<WORKER_MULTI, 0>::LEN, MandelPoint<WORKER_MULTI, 0>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          first_z.re.assign(precisionRecord->position.center.re);
          first_z.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&first_z);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*MandelPoint<MandelMath::number_a *>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*MandelPoint<WORKER_MULTI, 0>::LEN, precisionRecord->currentWorker.get(), indexOfWtiPoint, MandelPoint<WORKER_MULTI, 0>::LEN);
        }
      };
    }
  }
  for (int newy=(new_height-1)&0xfffffff; newy>=0; newy--) //avoid Clang warning about newy possibly ~ 2^31
  {
    int oldy=ytrans.transform(newy, old_height/2);
    first_z.im.assign(precisionRecord->position.center.im);
    first_z.im.add_double((new_height/2-newy)*precisionRecord->position.step_size);
    for (int newx=new_width-1; newx>=0; newx--)
    {
      int oldx=xtrans.transform(newx, old_width/2);
      if ((oldy<newy) || ((oldy==newy) && (oldx<=newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old to new
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].wstate==MandelPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=MandelPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*MandelPoint<MandelMath::number_a *>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*MandelPoint<MandelMath::number_a *>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*MandelPoint<WORKER_MULTI, 0>::LEN, old_sworker, (oldy*old_width+oldx)*MandelPoint<WORKER_MULTI, 0>::LEN, MandelPoint<WORKER_MULTI, 0>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          first_z.re.assign(precisionRecord->position.center.re);
          first_z.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(&first_z);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*MandelPoint<MandelMath::number_a *>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*MandelPoint<WORKER_MULTI, 0>::LEN, precisionRecord->currentWorker.get(), indexOfWtiPoint, MandelPoint<WORKER_MULTI, 0>::LEN);
        }
      }
    }
  }
}

void MandelModel::setView_double(double c_re, double c_im, double scale)
{
  MandelMath::complex<MandelMath::number_a *> c(precisionRecord->ntype);
  c.zero(c_re, c_im);
  setView(&c, scale);
}

void MandelModel::setView(const MandelMath::complex<MandelMath::number_a *> *c, double scale)
{
  MandelMath::complex<MandelMath::number_a *> old_c(precisionRecord->ntype);
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.setView(c, scale);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void MandelModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::complex<MandelMath::number_a *> old_c(precisionRecord->ntype);
  old_c.assign(&precisionRecord->position.center);

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int dx=qRound(delta_x);
    int dy=qRound(delta_y);
    precisionRecord->position.move(dx, dy);
    //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    startNewEpoch();
#if 0
  for (int newy=0; newy<imageHeight; newy++)
  {
    for (int newx=0; newx<imageWidth; newx++)
    {
      if (pointStore_[newy*imageWidth+newx].state==MandelPointStore::State::stWorking)
        dbgPoint();
    }
  }
#endif
  }
}

void MandelModel::zoom(double x, double y, int inlog)
{
  MandelMath::complex<MandelMath::number_a *> old_c(precisionRecord->ntype);
  old_c.assign(&precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, &old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void MandelModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int newLength=width*height;
    MandelPointStore *newStore=new MandelPointStore[newLength];
    {
      QString size_as_text=QString::number(newLength*sizeof(MandelPointStore));
      for (int pos=size_as_text.length()-3; pos>0; pos-=3)
        size_as_text.insert(pos, '\'');
      qDebug()<<"pointStore uses"<<size_as_text.toLocal8Bit().constData()<<"B"; //lots of work to skip those quotes... can't skip spaces at all
    }
    void *new_points=nullptr;
    switch (precisionRecord->ntype)
    {
      case MandelMath::NumberType::typeEmpty: goto lolwut;
      case MandelMath::NumberType::typeDouble: lolwut:
        new_points=new double[newLength*MandelPoint<MandelMath::number_a *>::LEN];
        break;
#if !NUMBER_DOUBLE_ONLY
      case MandelMath::NumberType::typeFloat128:
        new_points=new __float128[newLength];
        break;
      case MandelMath::NumberType::typeDDouble:
        new_points=new MandelMath::dd_real[newLength];
        break;
      case MandelMath::NumberType::typeQDouble:
        new_points=new MandelMath::dq_real[newLength];
        break;
      case MandelMath::NumberType::typeReal642:
        new_points=new MandelMath::real642[newLength];
        break;
#endif
    }
    MandelMath::complex<MandelMath::number_a *> old_c(precisionRecord->ntype);
    old_c.assign(&precisionRecord->position.center);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, &old_c,
                   new_points, newStore, width, height, &precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    switch (precisionRecord->ntype)
    {
      case MandelMath::NumberType::typeEmpty: goto ughwut;
      case MandelMath::NumberType::typeDouble: ughwut:
        delete[] (double *)precisionRecord->points;
        break;
#if !NUMBER_DOUBLE_ONLY
      case MandelMath::NumberType::typeFloat128:
        delete[] (__float128 *)precisionRecord->points;
        break;
      case MandelMath::NumberType::typeDDouble:
        delete[] (MandelMath::dd_real *)precisionRecord->points;
        break;
      case MandelMath::NumberType::typeQDouble:
        delete[] (MandelMath::dq_real *)precisionRecord->points;
        break;
      case MandelMath::NumberType::typeReal642:
        delete[] (MandelMath::real642 *)precisionRecord->points;
        break;
#endif
    }
    precisionRecord->points=new_points;
    delete[] pointStore;
    pointStore=newStore;
    imageWidth=width;
    imageHeight=height;

    startNewEpoch();
  }
}

void MandelModel::pause(bool pause)
{
  //MandelMath::complex old_c(precisionRecord->currentWorker.get());
  //old_c.assign(&precisionRecord->position.center);

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store

    //change working back to unknown
    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, &precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    if (!pause)
      startNewEpoch();
  }
}

void MandelModel::startNewEpoch()
{
  epoch=(epoch%2000000000)+1;
  nextGivenPointIndex=0;
  effortBonus=0;
#if 0
  for (int t=0; t<precisionRecord->threadCount; t++)
    if (precisionRecord->threads[t]->currentParams.pixelIndex<0)
      giveWork(precisionRecord->threads[t]);
#else
  for (int t=0; t<precisionRecord->threadCount; t++)
  {
    //giveWorkToThread(precisionRecord->threads[t]);
    precisionRecord->threads[t]->workIfEpoch=epoch;
  }
  _threadsWorking+=precisionRecord->threadCount;
  emit triggerComputeThreaded(epoch); //::invokeMethod cannot pass parameters, but ::connect can
#endif
}

void MandelModel::reimToPixel(int *circ_x, int *circ_y, const MandelMath::complex<MandelMath::number_a *> *c, MandelMath::number<MandelMath::number_a *> *tmp)
{
  double scaled;
  tmp->assign_across(&c->re);
  tmp->sub(precisionRecord->position.center.re);
  scaled=tmp->toDouble()/precisionRecord->position.step_size;
  if (scaled<-10003 || scaled>10003)
    *circ_x=-100;
  else
    *circ_x=qRound(scaled)+imageWidth/2;

  tmp->assign_across(&c->im);
  tmp->sub(precisionRecord->position.center.im);
  scaled=tmp->toDouble()/precisionRecord->position.step_size;
  if (scaled<-10003 || scaled>10003)
    *circ_y=-100;
  else
    *circ_y=imageHeight/2-qRound(scaled);
}

void MandelModel::paintOrbit(ShareableImageWrapper image, int x, int y)
{
  if ((x<0) || (x>=imageWidth) || (y<0) || (y>=imageHeight))
    return;
  QImage newOverlay(imageWidth, imageHeight, QImage::Format::Format_ARGB32);
  QPainter painter(&newOverlay);
  //QRgb what=newOverlay.pixel(0, 0);
  //if (what!=0) //0xcdcdcdcd in MSVC compiler
  {
    painter.setCompositionMode(QPainter::CompositionMode::CompositionMode_Source);
    painter.fillRect(0, 0, imageWidth, imageHeight, Qt::GlobalColor::transparent);
  };

  //MandelMath::number_instance_fix<MandelMath::worker_multi, PrecisionRecord::iiw_tmp> tmp(precisionRecord->currentWorker.get(), &precisionRecord->tmp_place);
/* testing how drawLine rounds in Qt
   result: extremely ugly and no way to fix it
  {
    painter.setPen(QColor(0xff, 0xff, 0xff));
    painter.fillRect(0, 0, 20+6*10, 20, Qt::GlobalColor::black);
    for (int i=0; i<6; i++)
    {
      newOverlay.setPixel(10+10*i, 10-5, 0x00ffffff);
      newOverlay.setPixel(10+10*i-1, 10+5, 0x00ffffff);
      newOverlay.setPixel(10+10*i+1, 10+5, 0x00ffffff);

      switch (i)
      {
        case 0: painter.drawLine(10+10*i-2, 10+3, 10+10*i, 10-3); break;
        case 1: painter.drawLine(10+10*i, 10-3, 10+10*i-2, 10+3); break;
        case 2: painter.drawLine(10+10*i+2, 10+3, 10+10*i, 10-3); break;
        case 3: painter.drawLine(10+10*i, 10-3, 10+10*i+2, 10+3); break;
        case 4:
        {
          QLine l3[3]={{10+10*i+1, 10-1, 10+10*i-1, 10-1},
                       {10+10*i-1, 10-1, 10+10*i-1, 10+1},
                       {10+10*i-1, 10+1, 10+10*i+1, 10+1}};
          painter.drawLines(l3, 3);
        } break;
        case 5:
        {
          QLine l2[2]={{10+10*i-2, 10, 10+10*i+2, 10},
                       {10+10*i, 10-2, 10+10*i, 10+2}};
          painter.drawLines(l2, 2);
        } break;
      }
    }
  }*/

  if (precisionRecord==nullptr)
  {
    dbgPoint(); //doesn't happen but makes compiler happy
    return;
  }
  //else if (precisionRecord->orbit.currentWorker->ntype()!=precisionRecord->currentWorker->ntype())
    //dbgPoint();

  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
  //MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
  //MandelPoint *data=&orbit_->evaluator.currentData;
  MandelMath::number<MandelMath::number_a *> tmp(precisionRecord->ntype);
  MandelPointStore *resultStore=&pointStore[y*imageWidth+x];
  if (!precisionRecord->lagu_c.is0())
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, &precisionRecord->lagu_c, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l3[3]={{circ_x+1, circ_y-1, circ_x-1, circ_y-1},
                   {circ_x-1, circ_y-1, circ_x-1, circ_y+1},
                   {circ_x-1, circ_y+1, circ_x+1, circ_y+1}};
      painter.drawLines(l3, 3);
    };

    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0, 0xff, 0)); //paint root as /\    .
    reimToPixel(&circ_x, &circ_y, &precisionRecord->lagu_r, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };
  }
  switch (resultStore->rstate)
  {
    case MandelPointStore::ResultState::stOutside:
    case MandelPointStore::ResultState::stOutAngle:
    { //exterior distance estimate
      int exterior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0xff, 0xff, 0));
      exterior=qRound(resultStore->exterior_hits/precisionRecord->position.step_size);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
      painter.setPen(QColor(0xc0, 0xc0, 0));
      exterior=qRound(resultStore->exterior_avoids/precisionRecord->position.step_size);
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
    } break;
    case MandelPointStore::ResultState::stPeriod2:
    case MandelPointStore::ResultState::stPeriod3:
    { //interior distance estimate
      int interior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0, 0xff, 0xff));
      interior=qRound(resultStore->interior/precisionRecord->position.step_size);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
      painter.setPen(QColor(0, 0xc0, 0xc0));
      interior=qRound(resultStore->interior/4/precisionRecord->position.step_size);
      painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);
    } break;
    default: ;
  }
  precisionRecord->position.pixelXtoRE(x-imageWidth/2, &precisionRecord->orbit.evaluator.currentParams.first_z.re);
  precisionRecord->position.pixelYtoIM(imageHeight/2-y, &precisionRecord->orbit.evaluator.currentParams.first_z.im);
  precisionRecord->orbit.evaluator.currentParams.epoch=epoch;
  precisionRecord->orbit.evaluator.workIfEpoch=precisionRecord->orbit.evaluator.busyEpoch;//epoch;
  precisionRecord->orbit.evaluator.currentParams.pixelIndex=0;
  precisionRecord->orbit.evaluator.currentParams.c.assign(&precisionRecord->orbit.evaluator.currentParams.first_z);
  //already 0 precisionRecord->orbit.evaluator.currentParams.nth_fz=0;
  //precisionRecord->orbit.evaluator.currentData.store->rstate=MandelPointStore::ResultState::stUnknown_;
  //precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stIdle;
  precisionRecord->orbit.evaluator.currentData.zero(&precisionRecord->orbit.evaluator.currentParams.first_z);
  precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stWorking;
  precisionRecord->orbit.evaluator.currentParams.breakOnNewNearest=true;
  precisionRecord->orbit.evaluator.currentParams.maxiter=1<<MAX_EFFORT;
  precisionRecord->orbit.evaluator.currentParams.want_extangle=(_selectedPaintStyle==paintStyleExterAngle) ||
                                                               (resultStore->rstate==MandelPointStore::ResultState::stOutAngle);
  {
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0xff, 0xff)); //paint path
    int line_sx, line_sy;
    reimToPixel(&line_sx, &line_sy, &precisionRecord->orbit.evaluator.currentData.f, &tmp);
    while ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stUnknown) &&
           (precisionRecord->orbit.evaluator.currentData.store->iter<(1<<MAX_EFFORT)))
    {
      int line_ex, line_ey;

      if ((resultStore->rstate==MandelPointStore::ResultState::stPeriod2 || resultStore->rstate==MandelPointStore::ResultState::stPeriod3) &&
          precisionRecord->orbit.evaluator.currentData.store->iter<resultStore->period)
      { //paint first period fully
        precisionRecord->orbit.evaluator.currentParams.maxiter=precisionRecord->orbit.evaluator.currentData.store->iter+1;
      }
      else if (precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess==0)
        precisionRecord->orbit.evaluator.currentParams.maxiter=1<<MAX_EFFORT; //dont't know->run fully
      else //stop at multiples of lookper, +1
        precisionRecord->orbit.evaluator.currentParams.maxiter=1+(precisionRecord->orbit.evaluator.currentData.store->iter/precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess+1)*precisionRecord->orbit.evaluator.currentData.store->lookper_lastGuess;
      precisionRecord->orbit.evaluator.thread.syncMandel();

      reimToPixel(&line_ex, &line_ey, &precisionRecord->orbit.evaluator.currentData.f, &tmp);
      if (line_ex>=-3 && line_ex<=10003 && line_ey>=-3 && line_ey<=10003)
      {
        if (line_sx>=-3 && line_sx<=10003 && line_sy>=-3 && line_sy<=10003)
          painter.drawLine(line_sx, line_sy, line_ex, line_ey);
        line_sx=line_ex;
        line_sy=line_ey;
      };
    }
  }
  if ((precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod2) ||
      (precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stPeriod3))
  {
    int circ_x, circ_y;
    painter.setPen(QColor(0, 0xff, 0xff)); //paint root
    reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.currentData.root, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };

    //orbit.bulbinfo
    /*MandelMath::complex c(orbit.worker, &orbit.evaluator.currentParams.c_re, &orbit.evaluator.currentParams.c_im, true);
    MandelMath::complex cb(orbit.worker, &orbit.bulb.cb_re, &orbit.bulb.cb_im, true);
    MandelMath::complex rb(orbit.worker, &orbit.bulb.rb_re_, &orbit.bulb.rb_im, true);
    MandelMath::complex xc(orbit.worker, &orbit.bulb.xc_re, &orbit.bulb.xc_im, true);
    MandelMath::complex baseZC(orbit.worker, &orbit.bulb.baseZC_re, &orbit.bulb.baseZC_im, true);
    MandelMath::complex baseCC(orbit.worker, &orbit.bulb.baseCC_re, &orbit.bulb.baseCC_im, true);*/
    //precisionRecord->orbit.bulb.foundMult=0;
    //precisionRecord->orbit.bulb.is_card=false;

    /*precisionRecord->orbit.bulb.valid=precisionRecord->orbit.evaluator.bulb.findBulbBase(precisionRecord->orbit.evaluator.currentData.store->period,
        &precisionRecord->orbit.evaluator.currentParams.c, &precisionRecord->orbit.bulb.cb, &precisionRecord->orbit.bulb.rb,
        &precisionRecord->orbit.bulb.xc, &precisionRecord->orbit.bulb.baseZC, &precisionRecord->orbit.bulb.baseCC,
        &precisionRecord->orbit.bulb.is_card, &precisionRecord->orbit.bulb.foundMult, &precisionRecord->orbit.bulb.baseFz);*/
    //precisionRecord->orbit.bulb.valid=
    precisionRecord->orbit.evaluator.bulb.findBulbBase(&precisionRecord->orbit.evaluator.tmp,
        precisionRecord->orbit.evaluator.currentData.store->period,
        &precisionRecord->orbit.evaluator.currentParams.c);
    precisionRecord->orbit.bulb.baseFz_.assign_across(&precisionRecord->orbit.evaluator.bulb.res_baseFz);
    precisionRecord->orbit.bulb.foundMult_=precisionRecord->orbit.evaluator.bulb.res_foundMult;
    if (precisionRecord->orbit.evaluator.bulb.res_valid)
    {
      painter.setBrush(QBrush(QColor(0, 0xff, 0xff)));
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.dbg_first_cb, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0x80, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.dbg_first_rb, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0x80, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.res_xc, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        if (precisionRecord->orbit.evaluator.bulb.res_card)
          painter.setPen(QColor(0, 0xff, 0xff)); //card center
        else
          painter.setPen(QColor(0x80, 0xc0, 0xc0)); //bulb center
        //painter.setBrush(Qt::BrushStyle::SolidPattern);
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }

      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.res_cb, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0, 0xff, 0)); //bulb base c
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.bulb.res_rb, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
      {
        painter.setPen(QColor(0xff, 0, 0)); //bulb base r
        painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      }
    };
  }

  //beginning of external ray
  if (_selectedPaintStyle==paintStyle::paintStyleExterAngle &&
      precisionRecord->orbit.evaluator.currentData.store->rstate==MandelPointStore::ResultState::stOutAngle)
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::SolidPattern);

    painter.setBrush(QBrush(QColor(0xff, 0xff, 0x80)));
    if (precisionRecord->orbit.evaluator.extangle.dbg.first_guess_valid>0)
    {
      painter.setPen(QColor(0xff, 0xff, 0xff)); //first guess of first step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.first_guess_0, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (precisionRecord->orbit.evaluator.extangle.dbg.last_guess_valid>0)
    {
      painter.setPen(QColor(0xff, 0xff, 0x80)); //last guess of first step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.last_guess_0, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    painter.setBrush(QBrush(QColor(0xff, 0xc0, 0x80)));
    if (precisionRecord->orbit.evaluator.extangle.dbg.first_guess_valid>1)
    {
      painter.setPen(QColor(0xff, 0xc0, 0xff)); //first guess of second step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.first_guess_1_, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (precisionRecord->orbit.evaluator.extangle.dbg.last_guess_valid>1)
    {
      painter.setPen(QColor(0xff, 0xc0, 0x80)); //last guess of second step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.last_guess_1_, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    painter.setBrush(QBrush(QColor(0xff, 0x80, 0x80)));
    if (precisionRecord->orbit.evaluator.extangle.dbg.first_guess_valid>2)
    {
      painter.setPen(QColor(0xff, 0x80, 0xff)); //first guess of third step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.first_guess_2, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (precisionRecord->orbit.evaluator.extangle.dbg.last_guess_valid>2)
    {
      painter.setPen(QColor(0xff, 0x80, 0x80)); //last guess of third step of external ray
      reimToPixel(&circ_x, &circ_y, &precisionRecord->orbit.evaluator.extangle.dbg.last_guess_2, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };
  };
  //donePixel1(&orbit.evaluator);
  //orbit.pointData.cleanup(position.worker);
  image.image->swap(newOverlay);
}

int MandelModel::periodToIndex(int period)
{
  //powers of factors -> color index
  //[1]..1
  //[1,1] .. 2
  //[2] .. 2
  //[1,1,1] .. 2
  //[2,1] .. 3
  //[3] .. 3
  //[1,1,1,1] .. 2
  //[2,1,1] .. 4
  //[2,2] .. 4
  //[3,1] .. 4
  //[4] .. 4
  //[1,1,1,1,1] .. 2
  //[2,1,1,1] .. 5
  //[2,2,1] .. 5
  //[3,1,1] .. 5
  //[3,2] .. 5
  //[4,1] .. 5
  //[5] .. 5
  if (period<=1)
    return 0;
  if (periodToIndexCache.length()<=period)
    periodToIndexCache.resize(period+1);
  if (periodToIndexCache[period]!=0)
    return periodToIndexCache[period];
  int c=1;
  for (int i=1; i<=period/2; i++)
  {
    if (period % i ==0)
    {
      int c2=periodToIndex(i);
      if (c2>=c)
        c=c2+1;
    };
  }
  periodToIndexCache[period]=c;
  return c;
}

int MandelModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return -1;
  timerWriteToImage.start();
  //int indexOfWtiPoint, _discard_;
  {
    qint64 totalNewtons=0;
    for (int t=0; t<precisionRecord->threadCount; t++)
      totalNewtons+=precisionRecord->threads[t]->totalNewtonIterations;
    //indexOfWtiPoint=totalNewtons; //zoom0: 47167, zoom1: 188490, zoom2: 754210  with breaker
    (void)totalNewtons;           //       47238         188929         756264  without breaker
  }
  //precisionRecord->wtiPoint.self_allocator._getFirstCapac(indexOfWtiPoint, _discard_);
  double extAngleZoom=1<<this->_extAngleZoom;
  MandelPointStore *wtiStore;
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      //MandelPoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiStore=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleKind:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              if (wtiStore->wstate==MandelPointStore::WorkState::stIdle)
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xffff00ff);
                //image.image->setPixel(x, y, 0xffffffff);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x9f+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b<<0));
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff8000ff);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffc000);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              /*int r;
              switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              if (wtiStore->period>wtiStore->nearziter_0)
                image.image->setPixel(x, y, 0xffff00ff); //seems to only happen by mistake, not in reality
              else
              {
                int index=periodToIndex(wtiStore->period);
                //reverse bottom 7 bits:
                int rh=0x73516240>>((index&7)<<2); //reverse bits 0..2
                int rl=0x73516240>>((index&0x70)>>2); //reverse bits 4..6
                rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
                image.image->setPixel(x, y, 0xff000000+(rh<<16));
              }
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleCls:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00906090);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              #if 0 //smooth color madness
              int r=128+floor(127*cos(data->iter/10.0*2*3.1415926535));
              int g=128+floor(127*cos((data->iter/10.0+0.333)*2*3.1415926535));
              int b=128+floor(127*cos((data->iter/10.0+0.667)*2*3.1415926535));
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                result.setPixel(x, y, 0xffffffff);
              else
                result.setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));*/
              #endif

              #if 0 //by data->iter only
              int iter=data->iter;
              switch (iter % 12)
              {
                case  0: image.image->setPixel(x, y, 0xff0000ff); break;
                case  1: image.image->setPixel(x, y, 0xff0080ff); break;
                case  2: image.image->setPixel(x, y, 0xff00ffff); break;
                case  3: image.image->setPixel(x, y, 0xff00ff80); break;
                case  4: image.image->setPixel(x, y, 0xff00ff00); break;
                case  5: image.image->setPixel(x, y, 0xff80ff00); break;
                case  6: image.image->setPixel(x, y, 0xffffff00); break;
                case  7: image.image->setPixel(x, y, 0xffff8000); break;
                case  8: image.image->setPixel(x, y, 0xffff0000); break;
                case  9: image.image->setPixel(x, y, 0xffff0080); break;
                case 10: image.image->setPixel(x, y, 0xffff00ff); break;
                case 11: image.image->setPixel(x, y, 0xff8000ff); break;
                default: image.image->setPixel(x, y, 0xffffffff);
              }
              #endif

              #if 1 //smooth by iter
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
              double re=precisionRecord->wtiPoint.f.re.toDouble();
              double im=precisionRecord->wtiPoint.f.im.toDouble();
              double iter=wtiStore->iter+6-log2(log2(re*re+im*im)); //+6 to match integer coloring

              //iter=sqrt(iter);
              iter=10*(log(1+iter/10));

              iter=iter/12;
              iter=(iter-floor(iter))*6;
              int iter_phase=iter;
              int iter_256=(iter-iter_phase)*256;
              int r=0xff, g=0xff, b=0xff;
              switch (iter_phase)
              {
                case  0: r=0x00; g=iter_256; b=0xff; break;
                case  1: r=0x00; g=0xff; b=0xff-iter_256; break;
                case  2: r=iter_256; g=0xff; b=0x00; break;
                case  3: r=0xff; g=0xff-iter_256; b=0x00; break;
                case  4: r=0xff; g=0x00; b=iter_256; break;
                case  5: r=0xff-iter_256; g=0x00; b=0xff; break;
              }
              image.image->setPixel(x, y, 0xff000000|(r<<16)|(g<<8)|b);
              #endif

            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              int index=periodToIndex(wtiStore->period);
              /*if (wtiStore->surehand==0 || wtiStore->surehand%wtiStore->period!=0)
              {
                image.image->setPixel(x, y, 0xffff0080);
              }
              else*/
              {
                //reverse bottom 7 bits:
                int rh=0x73516240>>((index&7)<<2);
                int rl=0x73516240>>((index&0x70)>>2);
                rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
                /*switch (data->period)
                {
                  case 1: r=0xc0; break;
                  case 2: r=0xff; break;
                  case 3: r=0x80; break;
                  default: r=0xe0;
                }*/
                //image.image->setPixel(x, y, 0xffffc0c0);
                image.image->setPixel(x, y, 0xff000000+rh*0x010101);
              }
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              //image.image->setPixel(x, y, 0xff808080);
              image.image->setPixel(x, y, 0xff000000);
            } break;
          }
        } break;
        case paintStyle::paintStyleExter:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              double tf;
              if ((wtiStore->exterior_avoids>10000) || (wtiStore->exterior_avoids<=0))
                image.image->setPixel(x, y, 0xff9f9f9f);
              else if (wtiStore->exterior_avoids>=1) //10000..1 -> -9999..0
              {
                tf=(1-wtiStore->exterior_avoids)*1;
                int g=0x9f+qRound(0x60*sin(tf*10)); //green fastest
                int b=0x9f+qRound(0x60*sin(tf/10)); //blue slowest
                int r=0x9f+qRound(0x60*sin(tf)); //red middle
                if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                  image.image->setPixel(x, y, 0xffffffff);
                else
                  image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
              else //1..0 -> 0..inf
              {
                tf=-log(wtiStore->exterior_avoids);//sqrt(1-log(wtiStore->exterior_avoids))*2-2;
                //tf=10*(log(1+tf/10));
                tf=(log(1+tf));
                int g=0x9f+qRound(0x60*sin(tf*10));
                int b=0x9f+qRound(0x60*sin(tf/10));
                //tf=(log(1+tf));
                int r=0x9f+qRound(0x60*sin(tf));
                //tf=(log(1+tf));
                //int b=0x9f+qRound(0x60*sin(tf));
                if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                  image.image->setPixel(x, y, 0xffffffff);
                else
                  image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
              double re=precisionRecord->wtiPoint.fz_r.re.toDouble();
              double im=precisionRecord->wtiPoint.fz_r.im.toDouble();
              //double angle=std::atan2(im, re);
              double mag=sqrt(re*re+im*im);
              int r, g, b;
              if (mag<=0)
               { r=0xff; g=0x00; b=0x00; }
              else
              {
                //re/=mag;
                //im/=mag;
                if (mag>=1)
                  mag=0.99;
                //mag=-log(1-mag)/log(2);
                int expo;
                mag=frexp(1-mag, &expo);
                r=0x40+(int)floor(0xc0*mag);
                g=r;//0xc0+qRound(0x3f*re);
                b=r;//0xc0+qRound(0x3f*im);
              }
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+b);
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleExterAngle:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case MandelPointStore::ResultState::stOutside:
              image.image->setPixel(x, y, 0x00606060);
              break;
            case MandelPointStore::ResultState::stOutAngle:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
              double tf=precisionRecord->wtiPoint.extangle.toDouble();
              if (tf==precisionRecord->orbit.evaluator.extangle.SPECIAL_VALUE_DEEP)
              {
                image.image->setPixel(x, y, 0xff808080);
              }
              else if (tf==precisionRecord->orbit.evaluator.extangle.SPECIAL_VALUE_EDGE)
              {
                image.image->setPixel(x, y, 0xff000000+(0xff<<16));
              }
              else
              {
                int r=(tf<-3.1415926536 || tf>3.1415926535)?0xff:0;//(tf<0)?0x80:0;
                //int r=0;//(tf<-3.1415926535 || tf>3.1415926535)?0x80:0;//qRound(log(wtiStore->iter)*100)%256;
                //int r=qRound(log(wtiStore->iter)*100)%128; //to much red would obscure the angle
                int g=0x9f+qRound(0x60*cos(tf*extAngleZoom));
                int b=0x9f+qRound(0x60*sin(tf*extAngleZoom)); //blue slowest

                /*
                int r=0x9f+qRound(0x60*sin(tf*extAngleZoom*2.828)); //red middle
                int g=0x9f+qRound(0x60*sin(tf*extAngleZoom*6.928)); //green fastest
                int b=0x9f+qRound(0x60*sin(tf*extAngleZoom)); //blue slowest
                */
                if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                  image.image->setPixel(x, y, 0xffffffff);
                else
                  image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
              double re=precisionRecord->wtiPoint.fz_r.re.toDouble();
              double im=precisionRecord->wtiPoint.fz_r.im.toDouble();
              //double angle=std::atan2(im, re);
              double mag=sqrt(re*re+im*im);
              int r, g, b;
              if (mag<=0)
               { r=0xff; g=0x00; b=0x00; }
              else
              {
                //re/=mag;
                //im/=mag;
                if (mag>=1)
                  mag=0.99;
                //mag=-log(1-mag)/log(2);
                int expo;
                mag=frexp(1-mag, &expo);
                r=0x40+(int)floor(0xc0*mag);
                g=r;//0xc0+qRound(0x3f*re);
                b=r;//0xc0+qRound(0x3f*im);
              }
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+b);
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleInter:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              //int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              int b=qRound(log(wtiStore->iter)*100)%256;
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              int ti=30;
              if ((wtiStore->interior>1) || (wtiStore->interior<=0))
                ti=0;
              else
                ti=(qRound(-log(wtiStore->interior/4)*300)+12*0xc0) % (6*0xc0);
              int r, g, b;
              if (ti<0xC0)
              { r=0x3f+ti; g=0xff; b=0x3f; }                           // + H L
              else if (ti<2*0xC0)                                      //
              { r=0xff; g=0xff-(ti-0xc0); b=0x3f; }                    // H - L
              else if (ti<3*0xC0)                                      //
              { r=0xff; g=0x3f; b=0x3f+(ti-2*0xc0); }                  // H L +
              else if (ti<4*0xC0)                                      //
              { r=0xff-(ti-3*0xc0); g=0x3f; b=0xff; }                  // - L H
              else if (ti<5*0xC0)                                      //
              { r=0x3f; g=0x3f+(ti-4*0xc0); b=0xff; }                  // L + H
              else                                                     //
              { r=0x3f; g=0xff; b=0xff-(ti-5*0xc0); }                  // L H -
              if ((r<0) || (r>255) || (g<0) || (g>255) || (b<0) || (b>255))
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleNear:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0xff000000);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int ti=wtiStore->nearziter_0;
              if (ti>=0)
              {
                int tj=0;
                while ((ti%2)==0) { tj+=128*2/5; ti/=2; }
                while ((ti%3)==0) { tj+=128*2/3; ti/=3; }
                while ((ti%5)==0) { tj+=30; ti/=5; }
                while ((ti%7)==0) { tj+=40; ti/=7; }
                while ((ti%11)==0) { tj+=50; ti/=11; }
                while ((ti%13)==0) { tj+=60; ti/=13; }
                while ((ti%17)==0) { tj+=70; ti/=17; }
                int b=0x80+(tj%0x80);
                image.image->setPixel(x, y, 0xff000000+(b<<0));
              }
              else
                image.image->setPixel(x, y, 0xff000080);
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              /*int ti=data->newton_iter; //very noisy, maybe show <=10, >10, >30, 49
              int r;
              if (ti<=10)
                r=0x60;
              else if (ti<30)
                r=0x90;
              else
                r=0xff;
              image.image->setPixel(x, y, 0xff000000+(r<<16));*/

              int ti=wtiStore->nearziter_0;
              if (ti>0)
              {
                int tj=0;
                while ((ti%2)==0) { tj+=128*2/5; ti/=2; }
                while ((ti%3)==0) { tj+=128*2/3; ti/=3; }
                while ((ti%5)==0) { tj+=30; ti/=5; }
                while ((ti%7)==0) { tj+=40; ti/=7; }
                while ((ti%11)==0) { tj+=50; ti/=11; }
                while ((ti%13)==0) { tj+=60; ti/=13; }
                while ((ti%17)==0) { tj+=70; ti/=17; }
                int r=0x80+(tj%0x80);
                image.image->setPixel(x, y, 0xff000000+(r<<16));
              }
              else
                image.image->setPixel(x, y, 0xff800000);
              /* we need func(2)!=func(3) here
              int index=periodToIndex(data->near0iter);
              //reverse bottom 7 bits:
              int rh=0x73516240>>((index&7)<<2);
              int rl=0x73516240>>((index&0x70)>>2);
              rh=0x80 | ((rh&0x07)<<4) | (index&0x08) | (rl&0x07);
              image.image->setPixel(x, y, 0xff000000+(rh<<16));*/
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleFZ:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
              double re=precisionRecord->wtiPoint.fz_r.re.toDouble();
              double im=precisionRecord->wtiPoint.fz_r.im.toDouble();
              double mag=sqrt(MandelMath::sqr_double(re)+MandelMath::sqr_double(im));
              if (mag<0) mag=0;
              else if (mag>1) mag=1;
              /*double mag2=mag*127.49+128;
              double phi=std::atan2(position.worker->toDouble(&data->fz_r_im), position.worker->toDouble(&data->fz_r_re))/(2*M_PI);
              if (phi<0) phi+=1;
              int g=qRound(phi*mag2);
              phi+=0.25;
              if (phi>=1) phi-=1;
              int b=qRound(phi*mag2);*/
              int g=mag==0?128:qRound(re/mag*127.49+128);
              //int b=mag==0?128:qRound(im/mag*127.49+128);
              int b=mag==0?128:qRound(im/mag*127.49+256)%256;
              int r=0;
              mag=-log(1-mag)/log(2);
              mag=mag-floor(mag);
              if (mag<0.2)
                r=128;
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleFC:
        {
          switch (wtiStore->rstate)
          {
            case MandelPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case MandelPointStore::ResultState::stOutside:
            case MandelPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case MandelPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case MandelPointStore::ResultState::stPeriod2:
            case MandelPointStore::ResultState::stPeriod3:
            {
              if (!wtiStore->has_fc_r)
                image.image->setPixel(x, y, 0xffc0c0c0);
              else
              {
                precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_a *>::LEN);
                double re=precisionRecord->wtiPoint.fc_c.re.toDouble();
                double im=precisionRecord->wtiPoint.fc_c.im.toDouble();
                double mag=std::hypot(re, im);//sqrt(MandelMath::sqr_double(re)+MandelMath::sqr_double(im));
                if (mag<0) mag=0;
                else if (mag>1) mag=1;
                //int magi=qRound(mag*127.49)+128;
                double mag2=mag*127.49+128;
                double phi=std::atan2(im, re)/(2*M_PI);
                if (phi<0) phi+=1;
                int g=qRound(phi*mag2);
                phi+=0.25;
                if (phi>=1) phi-=1;
                int b=qRound(phi*mag2);
                int r=0;
                mag=log(mag)/log(2);
                mag=mag-floor(mag);
                //if (mag<0.2)
                  //r=128;
                r=255*mag;
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case MandelPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
      }
    }
  //return result;
  qint64 elapsed=timerWriteToImage.elapsed();
  return elapsed;
}

void MandelModel::doneWorkInThread(MandelEvaluatorThread *)
{
  _threadsWorking--;
}

template <typename BASE>
int MandelModel::giveWorkThreaded(MandelEvaluator<BASE> *me)
{
  QReadLocker locker(&threading_mutex);
  me->timeInvokeSwitchTotal+=me->timeInvoke.nsecsElapsed();
  if (epoch!=me->busyEpoch)
    return 3;
  int retryEffortFrom=0;
  int nextEffortBonus=effortBonus; //don't jump to max instantly once pixels<threads
  //int intoEvaluator, _discard_;
  //me->currentData.self_allocator._getFirstCapac(intoEvaluator, _discard_);
  while (retryEffortFrom>=0)
  {
    retryEffortFrom=-1;
    for (int pi=0; pi<imageWidth*imageHeight; pi++)
    {
      int pointIndex=(nextGivenPointIndex+pi)%(imageWidth*imageHeight);
      //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
        //dbgPoint();
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), pointIndex*MandelPoint::LEN, MandelPoint::LEN, nullptr);
      //MandelPoint pointData_(&pointStore_[pointIndex], &allo);
      MandelPointStore *storeAtIndex=&pointStore[pointIndex];
      MandelPointStore::WorkState state_expected=MandelPointStore::WorkState::stIdle;
      int extra_effort=0;
      if (!storeAtIndex->wstate.compare_exchange_strong(state_expected, MandelPointStore::WorkState::stWorking))
      {
        if ((_selectedPaintStyle==paintStyleFC) &&
            (!storeAtIndex->has_fc_r) &&
            (storeAtIndex->rstate==MandelPointStore::ResultState::stPeriod2 ||
             storeAtIndex->rstate==MandelPointStore::ResultState::stPeriod3))
        {
          state_expected=MandelPointStore::WorkState::stDone;
          if (!storeAtIndex->wstate.compare_exchange_strong(state_expected, MandelPointStore::WorkState::stWorking))
            continue;
            //should check if it still needs but it should be quite rare, evaluate tests it anyway
          extra_effort=1;
        }
        else if (_selectedPaintStyle==paintStyleExterAngle &&
                 storeAtIndex->rstate==MandelPointStore::ResultState::stOutside)
        {
          state_expected=MandelPointStore::WorkState::stDone;
          if (!storeAtIndex->wstate.compare_exchange_strong(state_expected, MandelPointStore::WorkState::stWorking))
            continue;
            //should check if it still needs but it should be quite rare, evaluate tests it anyway
          extra_effort=1;
        }
        else
          continue;
      }
      {
        if (me->currentParams.pixelIndex!=-1)
          dbgPoint();
        working_assert(me->currentParams.pixelIndex==-1);
        {
          int phasex=(pointIndex%imageWidth-imageWidth/2+precisionRecord->position.cached_center_re_mod+32768)%32768;
          int phasey=(pointIndex/imageWidth-imageHeight/2+precisionRecord->position.cached_center_im_mod+32768)%32768;
          //int effort=ctz16(phasex)+ctz16(phasey);
          int effort=MandelMath::ctz16(phasex | phasey);
          if (effort>8)
            effort=8;
          effort+=nextEffortBonus;
          if (effort>=MAX_EFFORT)
            effort=MAX_EFFORT;
          if (effort<2)
            effort=2; //quick run ... at least 4 iterations
          me->currentParams.maxiter=1<<effort;
          if (storeAtIndex->iter >= me->currentParams.maxiter+extra_effort)
          {
            if (effort>=MAX_EFFORT)
            {
              storeAtIndex->rstate=MandelPointStore::ResultState::stMaxIter;
              storeAtIndex->wstate=MandelPointStore::WorkState::stDone;
            }
            else
            {
              if (retryEffortFrom<0)
                retryEffortFrom=pointIndex;
              storeAtIndex->wstate=MandelPointStore::WorkState::stIdle;
            }
          }
          else
          {
            //if (me->currentWorker->ntype()!=precisionRecord->ntype)
              //dbgPoint();
            //evaluator->switchType(position.worker);
            precisionRecord->position.pixelXtoRE<BASE>(pointIndex%imageWidth - imageWidth/2, &me->currentParams.first_z.re);
            precisionRecord->position.pixelYtoIM<BASE>(imageHeight/2-pointIndex/imageWidth, &me->currentParams.first_z.im);
            me->currentParams.c.assign(&me->currentParams.first_z);
            me->currentParams.epoch=me->busyEpoch;
            //storeAtIndex->state=MandelPointStore::State::stWorking;
            me->currentParams.pixelIndex=pointIndex;
            me->currentParams.want_fc_r=(_selectedPaintStyle==paintStyleFC);
            me->currentParams.want_extangle=(_selectedPaintStyle==paintStyleExterAngle);
#if CURRENT_STORE_DIRECT
            me->currentData.store=storeAtIndex;
#else
            me->currentDataStore.assign(storeAtIndex);
#endif
            me->currentData.readFrom(precisionRecord->points, pointIndex*MandelPoint<MandelMath::number_a *>::LEN);
            nextGivenPointIndex=(pointIndex+1)%(imageWidth*imageHeight);
            effortBonus=nextEffortBonus;
            return 0;
          }
        }
      }
        //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
    }
    if ((retryEffortFrom>=0) && (nextEffortBonus<MAX_EFFORT))
    {
      nextEffortBonus++;
      nextGivenPointIndex=retryEffortFrom;
    }
    else
      retryEffortFrom=-1;
  }
  return 2;
}

template <typename BASE>
int MandelModel::doneWorkThreaded(MandelEvaluator<BASE> *me, bool giveWork)
{
  QReadLocker locker(&threading_mutex);
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    MandelPointStore *dstStore=&pointStore[me->currentParams.pixelIndex];
#if CURRENT_STORE_DIRECT
#else
    if (dstStore->wstate!=MandelPointStore::WorkState::stWorking)
      dbgPoint(); //leftovers should be from different epoch
#endif
    {
      //int first, capac;
      me->currentData.writeTo(precisionRecord->points, me->currentParams.pixelIndex*MandelPoint<MandelMath::number_a *>::LEN);
#if CURRENT_STORE_DIRECT
#else
      dstStore->assign(me->currentData.store);
#endif
    }
    if (dstStore->wstate.load()==MandelPointStore::WorkState::stIdle)
      dbgPoint();
    if (dstStore->wstate==MandelPointStore::WorkState::stWorking)
    {
       if (dstStore->iter>=(1<<MAX_EFFORT))
       {
         dstStore->rstate=MandelPointStore::ResultState::stMaxIter;
         dstStore->wstate=MandelPointStore::WorkState::stDone;
       }
       else if (dstStore->rstate!=MandelPointStore::ResultState::stUnknown)
         dstStore->wstate=MandelPointStore::WorkState::stDone;
       else
         dstStore->wstate=MandelPointStore::WorkState::stIdle;
    }
    else
      dbgPoint();
  }
  else if (me->currentParams.epoch!=epoch)
  {
    //qDebug()<<"Old pixel finished";
    /*problem fixed if ((me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
    {
      MandelPointStore *dstStore=&pointStore_[me->currentParams.pixelIndex];
      if (dstStore->state==MandelPointStore::State::stWorking)
        dbgPoint();
    }*/
  }
  else
    qWarning()<<"Invalid pixel finished";
  me->currentParams.pixelIndex=-1;
  if (me->currentParams.epoch==epoch)
  {
    if (giveWork)
      return giveWorkThreaded(me);
    else
      return 0;
  }
  return 1;
}

void MandelModel::selectedPrecisionChanged()
{
  int pointCount;
  if (imageWidth<=0 || imageHeight<=0) //Qt begins with width=0, height=-13
    pointCount=0;
  else
    pointCount=imageWidth*imageHeight;
  void *new_points;
  MandelMath::NumberType newPrecision;
  switch (_selectedPrecision)
  {
    using Type=MandelMath::NumberType;
    case precisionDouble: lolwut:
      newPrecision=Type::typeDouble;
      new_points=MandelMath::number<double>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*MandelPoint<double>::LEN);
      break;
#if !NUMBER_DOUBLE_ONLY
    case precisionFloat128:
      newPrecision=Type::typeFloat128;
      new_points=MandelMath::number<__float128>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*MandelPoint<__float128>::LEN);
      break;
    case precisionDDouble:
      newPrecision=Type::typeDDouble;
      new_points=MandelMath::number<MandelMath::dd_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*MandelPoint<MandelMath::dd_real>::LEN);
      break;
    case precisionQDouble:
      newPrecision=Type::typeQDouble;
      new_points=MandelMath::number<MandelMath::dq_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*MandelPoint<MandelMath::dq_real>::LEN);
      break;
    case precisionReal642:
      newPrecision=Type::typeReal642;
      new_points=MandelMath::number<MandelMath::real642>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*MandelPoint<MandelMath::real642>::LEN);
      break;
#endif
    default: goto lolwut;
  }

  {
    PrecisionRecord *newPrecisionRecord=new PrecisionRecord(newPrecision, precisionRecord, this);
    newPrecisionRecord->points=new_points;
    delete precisionRecord;
    precisionRecord=newPrecisionRecord;
  }

  //pointStore stays

  {
    QWriteLocker locker(&threading_mutex);
    startNewEpoch();
  }
}


MandelModel::Position::Position(MandelMath::NumberType ntype, const Position *source):
  center(ntype)
{
  if (source)
  {
    center.assign_across(&source->center);
    step_log=source->step_log;
    step_size=source->step_size;
  }
  else
  {
    center.zero(-0.5, 0.0);
    step_log=7;
    step_size=1.0/128;
  }
  updateCachedDepth();
}

MandelModel::Position::~Position()
{
}

void MandelModel::Position::setView(const MandelMath::complex<MandelMath::number_a *> *c, double scale)
{
  step_log=-ilogb(scale);
  step_size=ldexp(1.0, -step_log);
  center.assign(c);
  center.lshift(step_log);
  center.re.round();
  center.im.round();
  center.lshift(-step_log);
  updateCachedDepth();
}

void MandelModel::Position::move(int delta_x, int delta_y)
{
  //qDebug()<<"move ("<<delta_x<<","<<delta_y<<")";
  center.re.add_double(-delta_x*step_size);
  center.im.add_double(+delta_y*step_size);
  cached_center_re_mod=(cached_center_re_mod-delta_x+32768)%32768;
  cached_center_im_mod=(cached_center_im_mod+delta_y+32768)%32768;
#if UPDATE_CACHED_MOD
  int ccrm=cached_center_re_mod;
  int ccim=cached_center_im_mod;
  updateCachedDepth();
  if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
    dbgPoint();
#endif
}

void MandelModel::Position::scale(int inlog, int center_x, int center_y)
{
  //qDebug()<<"scale "<<inlog<<" @"<<center_x<<","<<center_y<<"";
  /*
  center_re+center_x*step_size = new_center_re+center_x*new_step_size
  new_step_size=step_size/2^inlog

  center_re+center_x*(step_size-new_step_size) = new_center_re
  */
  if (inlog==0)
    return;
  if (inlog>0)
  {
    int max_zoom=qRound(2-log(center.re.eps2())/log(4));
    if (step_log+inlog>max_zoom)//MAX_ZOOM_IN_DOUBLE)
      return;
    double old_step_size=step_size;
    step_log+=inlog;
    for (int i=0; i<inlog; i++)
    {
      step_size/=2;
    }
    center.re.add_double(center_x*(old_step_size-step_size));
    center.im.add_double(center_y*(old_step_size-step_size));
    //center_re/step_size=center_re/old_step_size*old_step_size/step_size+center_x*(old_step_size-step_size)/step_size;
    cached_center_re_mod=cached_center_re_mod*(1<<inlog)+center_x*((1<<inlog)-1);
    cached_center_re_mod&=0x7fff;//=(cached_center_re_mod+32768)%32768;
    cached_center_im_mod=cached_center_im_mod*(1<<inlog)+center_y*((1<<inlog)-1);
    cached_center_im_mod&=0x7fff;//=(cached_center_im_mod+32768)%32768;
#if UPDATE_CACHED_MOD
    int ccrm=cached_center_re_mod;
    int ccim=cached_center_im_mod;
    updateCachedDepth();
    if ((ccrm!=cached_center_re_mod) || (ccim!=cached_center_im_mod))
      dbgPoint();
#endif
  }
  else
  {
    double old_step_size=step_size;
    int adjust_x=(cached_center_re_mod+center_x)&((1<<-inlog)-1);
    if (adjust_x&(1<<(1-inlog)))
      adjust_x-=(1<<-inlog);
    int adjust_y=(cached_center_im_mod+center_y)&((1<<-inlog)-1);
    if (adjust_y&(1<<(1-inlog)))
      adjust_y-=(1<<-inlog);
    step_log+=inlog;
    for (int i=0; i<-inlog; i++)
    {
      step_size*=2;
    }
    center.re.add_double((center_x-adjust_x)*(old_step_size-step_size)); //(old_step_size-step_size)=(1-(1<<-inlog))*old_step_size
    center.im.add_double((center_y-adjust_y)*(old_step_size-step_size));
    //need to roll in high bits anyway
    /*cached_center_re_mod=cached_center_re_mod*(old_step_size/step_size)+(center_x-adjust_y)*(old_step_size/step_size-1);
    cached_center_re_mod%=32768;
    cached_center_im_mod=cached_center_im_mod*(old_step_size/step_size)+(center_y-adjust_y)*(old_step_size/step_size-1);
    cached_center_im_mod%=32768;*/
    updateCachedDepth();
  }
}

void MandelModel::Position::updateCachedDepth()
{
  MandelMath::complex<MandelMath::number_a *> d(center.re.ntype());

  d.assign(&center);
  d.lshift(step_log-15);
  d.re.mod1();
  d.im.mod1();
  d.lshift(15);
  cached_center_re_mod=d.re.toRound();
  cached_center_im_mod=d.im.toRound();
}

template <typename BASE>
void MandelModel::Position::pixelXtoRE(int x, MandelMath::number<BASE> *result)
{
  result->assign_across(&center.re);
  result->add_double(x*step_size);
}

template <typename BASE>
void MandelModel::Position::pixelYtoIM(int y, MandelMath::number<BASE> *result)
{
  result->assign_across(&center.im);
  result->add_double(y*step_size);
}



MandelModel::Orbit::Orbit(MandelMath::NumberType ntype):
  evaluator(ntype, true), bulb(ntype)
{
}

MandelModel::Orbit::~Orbit()
{
  evaluator.workIfEpoch=-1;
  evaluator.thread.quit();
  evaluator.thread.wait(1000);
}

MandelModel::Orbit::Bulb::Bulb(MandelMath::NumberType ntype):
  cb_(ntype), rb_(ntype), xc_(ntype), baseZC_(ntype), baseCC_(ntype), baseFz_(ntype)
{
}

MandelModel::Orbit::Bulb::~Bulb()
{
}

MandelModel::PrecisionRecord::PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, MandelModel *doneReceiver):
  ntype(ntype), wtiPoint(nullptr, ntype),
  position(ntype, source?&source->position:nullptr),
  orbit(ntype),//, source?&source->orbit:nullptr),
  lagu_c(ntype), lagu_r(ntype), tmp_place(ntype),
  threadCount(source?source->threadCount:0), threads(nullptr)
{
  if (source)
  {
  }
  else
  {
    lagu_c.zero(0, 0);
    lagu_r.zero(0, 0);
  }
  if (threadCount<=0)
  { //first init
    threadCount=QThread::idealThreadCount()-1;
    //threadCount=1;
    if (threadCount<1)
      threadCount=1;
  };
  switch (ntype)
  {
    case MandelMath::NumberType::typeEmpty:
    case MandelMath::NumberType::typeDouble:
      threads=(MandelEvaluator<MandelMath::number_a *> **)new MandelEvaluator<double> *[threadCount];
      for (int t=0; t<threadCount; t++)
      {
        MandelEvaluator<double> *thread=new MandelEvaluator<double>(ntype, source==nullptr);
        threads[t]=(MandelEvaluator<MandelMath::number_a *> *)thread;
        //Qt not allowing parameters in invokeMethod makes it really easy... pfft
        thread->threaded.give=[doneReceiver](MandelEvaluator<double> *me)
            {
              //return doneReceiver->giveWorkThreaded((MandelEvaluator<MandelMath::worker_multi_double> *)me);
              return doneReceiver->giveWorkThreaded<double>(me);
            };
        thread->threaded.doneMandel=[doneReceiver](MandelEvaluator<double> *me, bool giveWork)
            {
              return doneReceiver->doneWorkThreaded<double>(me, giveWork);
            };
        QObject::connect(&thread->thread, &MandelEvaluatorThread::doneMandelThreaded,
                         doneReceiver, &MandelModel::doneWorkInThread,
                         Qt::ConnectionType::QueuedConnection);
        QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                         &thread->thread, &MandelEvaluatorThread::doMandelThreaded,
                         Qt::ConnectionType::QueuedConnection);
      }
      break;
#if !NUMBER_DOUBLE_ONLY
    case MandelMath::NumberType::typeFloat128:
      threads=(MandelEvaluator<MandelMath::number_a *> **)new MandelEvaluator<__float128> *[threadCount];
      for (int t=0; t<threadCount; t++)
      {
        MandelEvaluator<__float128> *thread=new MandelEvaluator<__float128>(ntype, source==nullptr);
        threads[t]=(MandelEvaluator<MandelMath::number_a *> *)thread;
        //Qt not allowing parameters in invokeMethod makes it really easy... pfft
        thread->threaded.give=[doneReceiver](MandelEvaluator<__float128> *me)
            {
              //return doneReceiver->giveWorkThreaded((MandelEvaluator<MandelMath::worker_multi_double> *)me);
              return doneReceiver->giveWorkThreaded<__float128>(me);
            };
        thread->threaded.doneMandel=[doneReceiver](MandelEvaluator<__float128> *me, bool giveWork)
            {
              return doneReceiver->doneWorkThreaded<__float128>(me, giveWork);
            };
        QObject::connect(&thread->thread, &MandelEvaluatorThread::doneMandelThreaded,
                         doneReceiver, &MandelModel::doneWorkInThread,
                         Qt::ConnectionType::QueuedConnection);
        QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                         &thread->thread, &MandelEvaluatorThread::doMandelThreaded,
                         Qt::ConnectionType::QueuedConnection);
      }
      break;
    case MandelMath::NumberType::typeDDouble:
      threads=(MandelEvaluator<MandelMath::number_a *> **)new MandelEvaluator<MandelMath::dd_real> *[threadCount];
      for (int t=0; t<threadCount; t++)
      {
        MandelEvaluator<MandelMath::dd_real> *thread=new MandelEvaluator<MandelMath::dd_real>(ntype, source==nullptr);
        threads[t]=(MandelEvaluator<MandelMath::number_a *> *)thread;
        //Qt not allowing parameters in invokeMethod makes it really easy... pfft
        thread->threaded.give=[doneReceiver](MandelEvaluator<MandelMath::dd_real> *me)
            {
              //return doneReceiver->giveWorkThreaded((MandelEvaluator<MandelMath::worker_multi_double> *)me);
              return doneReceiver->giveWorkThreaded<MandelMath::dd_real>(me);
            };
        thread->threaded.doneMandel=[doneReceiver](MandelEvaluator<MandelMath::dd_real> *me, bool giveWork)
            {
              return doneReceiver->doneWorkThreaded<MandelMath::dd_real>(me, giveWork);
            };
        QObject::connect(&thread->thread, &MandelEvaluatorThread::doneMandelThreaded,
                         doneReceiver, &MandelModel::doneWorkInThread,
                         Qt::ConnectionType::QueuedConnection);
        QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                         &thread->thread, &MandelEvaluatorThread::doMandelThreaded,
                         Qt::ConnectionType::QueuedConnection);
      }
      break;
    case MandelMath::NumberType::typeQDouble:
      threads=(MandelEvaluator<MandelMath::number_a *> **)new MandelEvaluator<MandelMath::dq_real> *[threadCount];
      for (int t=0; t<threadCount; t++)
      {
        MandelEvaluator<MandelMath::dq_real> *thread=new MandelEvaluator<MandelMath::dq_real>(ntype, source==nullptr);
        threads[t]=(MandelEvaluator<MandelMath::number_a *> *)thread;
        //Qt not allowing parameters in invokeMethod makes it really easy... pfft
        thread->threaded.give=[doneReceiver](MandelEvaluator<MandelMath::dq_real> *me)
            {
              //return doneReceiver->giveWorkThreaded((MandelEvaluator<MandelMath::worker_multi_double> *)me);
              return doneReceiver->giveWorkThreaded<MandelMath::dq_real>(me);
            };
        thread->threaded.doneMandel=[doneReceiver](MandelEvaluator<MandelMath::dq_real> *me, bool giveWork)
            {
              return doneReceiver->doneWorkThreaded<MandelMath::dq_real>(me, giveWork);
            };
        QObject::connect(&thread->thread, &MandelEvaluatorThread::doneMandelThreaded,
                         doneReceiver, &MandelModel::doneWorkInThread,
                         Qt::ConnectionType::QueuedConnection);
        QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                         &thread->thread, &MandelEvaluatorThread::doMandelThreaded,
                         Qt::ConnectionType::QueuedConnection);
      }
      break;
    case MandelMath::NumberType::typeReal642:
      threads=(MandelEvaluator<MandelMath::number_a *> **)new MandelEvaluator<MandelMath::real642> *[threadCount];
      for (int t=0; t<threadCount; t++)
      {
        MandelEvaluator<MandelMath::real642> *thread=new MandelEvaluator<MandelMath::real642>(ntype, source==nullptr);
        threads[t]=(MandelEvaluator<MandelMath::number_a *> *)thread;
        //Qt not allowing parameters in invokeMethod makes it really easy... pfft
        thread->threaded.give=[doneReceiver](MandelEvaluator<MandelMath::real642> *me)
            {
              //return doneReceiver->giveWorkThreaded((MandelEvaluator<MandelMath::worker_multi_double> *)me);
              return doneReceiver->giveWorkThreaded<MandelMath::real642>(me);
            };
        thread->threaded.doneMandel=[doneReceiver](MandelEvaluator<MandelMath::real642> *me, bool giveWork)
            {
              return doneReceiver->doneWorkThreaded<MandelMath::real642>(me, giveWork);
            };
        QObject::connect(&thread->thread, &MandelEvaluatorThread::doneMandelThreaded,
                         doneReceiver, &MandelModel::doneWorkInThread,
                         Qt::ConnectionType::QueuedConnection);
        QObject::connect(doneReceiver, &MandelModel::triggerComputeThreaded,
                         &thread->thread, &MandelEvaluatorThread::doMandelThreaded,
                         Qt::ConnectionType::QueuedConnection);
      }
      break;
#endif
  }
}

MandelModel::PrecisionRecord::~PrecisionRecord()
{
  for (int t=threadCount-1; t>=0; t--)
  {
    threads[t]->workIfEpoch=-1;
    threads[t]->thread.quit();
  }
  for (int t=threadCount-1; t>=0; t--)
  {
    threads[t]->thread.wait(1000);
    delete threads[t];
  }
  delete[] threads;
  threadCount=0;
  threads=nullptr;
}
