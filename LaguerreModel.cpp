#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "LaguerreModel.hpp"
#include "MandelEvaluator.hpp"

#define working_assert(x) { if (!(x)) dbgPoint(); }

#define UPDATE_CACHED_MOD 0 //works until precision doesn't allow

LaguerreModel::LaguerreModel(): QObject(),
  pointStore(nullptr), precisionRecord(nullptr)
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);
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
  QObject::connect(this, &LaguerreModel::selectedPrecisionChange,
                   this, &LaguerreModel::selectedPrecisionChanged);
  selectedPrecisionChanged();
}

LaguerreModel::~LaguerreModel()
{
  epoch=(epoch%2000000000)+1;
  delete precisionRecord;

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

void LaguerreModel::startRunning()
{
  MandelMath::visit([](auto &threads, int threadCount){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
        for (int t=0; t<threadCount; t++)
        {
          threads[t]->startRunning();
        }
      };
  }, precisionRecord->threads, precisionRecord->threadCount);
}

QString LaguerreModel::pixelXtoRE_str(int x)
{
  MandelMath::number<MandelMath::number_any> num(precisionRecord->position.center.re);
  num.add_double((x - imageWidth/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString LaguerreModel::pixelYtoIM_str(int y)
{
  MandelMath::number<MandelMath::number_any> num(precisionRecord->position.center.im);
  num.add_double((y - imageHeight/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString LaguerreModel::getTimes()
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
  result="X-X";
  return result;
}

QString LaguerreModel::getTextXY()
{
  if (precisionRecord==nullptr)
    return "-";
  return precisionRecord->orbit.evaluator.currentParams.mandel.first_z.re.toString()+" +i* "+
         precisionRecord->orbit.evaluator.currentParams.mandel.first_z.im.toString();
}

/*static QString mandDoubleToString(double x)
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
}*/

QString LaguerreModel::getTextInfoGen()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_any> tmp(&precisionRecord->orbit.evaluator.tmp);
    reimToPixel(&orbit_x, &orbit_y, precisionRecord->orbit.evaluator.currentParams.mandel.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
  //LaguerrePoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  LaguerrePointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;
  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*LaguerrePoint<MandelMath::number_any>::LEN);
  LaguerrePoint<MandelMath::number_any> *data=&precisionRecord->wtiPoint;

  QString state;
  switch (data_store->state.load())
  {
    case LaguerrePointStore::State::stUnknown:
      state="Unk"; break;
    case LaguerrePointStore::State::stWorking:
      state="Wor"; break;
    case LaguerrePointStore::State::stResolved:
      state="OK"; break;
    case LaguerrePointStore::State::stFail:
      state="Err"; break;
  }

  return "Per="+QString::number(precisionRecord->params.period)+" "+state+" iter="+QString::number(data->store->iter)+
        " mu="+QString::number(precisionRecord->orbit.first_mu_re, 'f', 3)+","+QString::number(precisionRecord->orbit.first_mu_im, 'f', 3)+
        " mum="+QString::number(precisionRecord->orbit.first_mum_re, 'f', 3)+","+QString::number(precisionRecord->orbit.first_mum_im, 'f', 3);
}

QString LaguerreModel::getTextInfoSpec()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_any> tmp(&precisionRecord->orbit.evaluator.tmp);
    reimToPixel(&orbit_x, &orbit_y, precisionRecord->orbit.evaluator.currentParams.mandel.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
  //LaguerrePoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  LaguerrePointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;

  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*LaguerrePoint<MandelMath::number_any>::LEN);
  //LaguerrePoint<MandelMath::number_any> *data=&precisionRecord->wtiPoint;

  switch (data_store->state.load())
  {
    case LaguerrePointStore::State::stUnknown:
      return " ";
      break;
    case LaguerrePointStore::State::stWorking:
      return "Working...";
      break;
    case LaguerrePointStore::State::stResolved:
    {
      //MandelMath::complex fz(orbit.worker, &orbit.pointData.fz_r_re, &orbit.pointData.fz_r_im, true);
      QString naiveChoice;
      switch (precisionRecord->orbit.pointData.store->naiveChoice)
      {
        case NewtonNaiveChoice::nc03: naiveChoice="0.3"; break;
        case NewtonNaiveChoice::nc05: naiveChoice="0.5"; break;
        case NewtonNaiveChoice::nc08: naiveChoice="0.8"; break;
        case NewtonNaiveChoice::ncWide: naiveChoice="Wide"; break;
        case NewtonNaiveChoice::nc100: naiveChoice="100"; break;
        case NewtonNaiveChoice::nc90: naiveChoice="90"; break;
        case NewtonNaiveChoice::nc80: naiveChoice="80"; break;
        case NewtonNaiveChoice::nc60: naiveChoice="60"; break;
        case NewtonNaiveChoice::ncClose: naiveChoice="Close"; break;
      }
      switch (_selectedPaintStyle)
      {
        case paintStyleCls:
          return QString("attr=")+QString::number(precisionRecord->orbit.pointData.fz_r.getMag_double())+
                 QString(" firstM=")+QString::number(precisionRecord->orbit.pointData.store->firstM)+
                 QString(" NChoice=")+naiveChoice;
        case paintStyleNthFz:
          return QString::number(precisionRecord->params.nth_fz)+
                 QString("-th f'=")+precisionRecord->orbit.pointData.nth_fz.toString();
        case paintStyleNthFz1:
          return QString::number(precisionRecord->params.nth_fz)+
                 QString("-th f'=")+precisionRecord->orbit.pointData.nth_fz.toString();
        case paintStyleRoots:
          return QString("attr=")+QString::number(precisionRecord->orbit.pointData.fz_r.getMag_double());
      }

    } break;
    case LaguerrePointStore::State::stFail:
      return " ";
  }
  return "-?-?-";
}

void LaguerreModel::recomputeRoot(int max_effort)
{
  MandelEvaluator<MandelMath::number_any> &evaluator=precisionRecord->orbit.evaluator;
  evaluator.currentParams.mandel.first_z.assign(precisionRecord->params.c);
  evaluator.currentParams.epoch=epoch;
  evaluator.workIfEpoch=evaluator.busyEpoch;//epoch;
  evaluator.currentParams.pixelIndex=0;
  evaluator.currentParams.mandel.c.assign(evaluator.currentParams.mandel.first_z);
  //already 0 precisionRecord->orbit.evaluator.currentParams.nth_fz=0;
  //precisionRecord->orbit.evaluator.currentData.store->rstate=MandelPointStore::ResultState::stUnknown_;
  //precisionRecord->orbit.evaluator.currentData.store->wstate=MandelPointStore::WorkState::stIdle;
  evaluator.mandelData.zero(evaluator.currentParams.mandel.first_z);
  evaluator.mandelData.store->wstate=MandelPointStore::WorkState::stWorking;
  evaluator.currentParams.breakOnNewNearest=false;
  evaluator.currentParams.maxiter=1<<max_effort;
  evaluator.currentParams.want_extangle=false;//don't need for root
  evaluator.thread.syncMandel();
  if (evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stPeriod2 ||
      evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stPeriod3)
  {
    precisionRecord->params.period=evaluator.mandelData.store->near0iter_1;//period;
    precisionRecord->params.root.assign(evaluator.mandelData.root);
  }
  else if (evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stMaxIter)
  {
    precisionRecord->params.period=evaluator.mandelData.store->near0iter_1;
    precisionRecord->params.root.assign(evaluator.mandelData.f);
  }
  else
  {
    precisionRecord->params.period=1;
    precisionRecord->params.root.zero(0, 0);
  }

  MandelMath::visit([](auto &threads, int threadCount, MandelEvaluator<MandelMath::number_any> &evaluator){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
        //using Evaluator=std::remove_pointer_t<std::remove_pointer_t<std::remove_reference_t<decltype(threads)>>>;
        for (int t=0; t<threadCount; t++)
        {
          threads[t]->currentParams.mandel.assign_across(evaluator.currentParams.mandel);
        }
      };
  }, precisionRecord->threads, precisionRecord->threadCount, evaluator);
}

void LaguerreModel::setParams(ShareableViewInfo viewInfo)
{
  //position.setNumberType(viewInfo.worker->ntype());
  if (precisionRecord==nullptr)
  {
    dbgPoint();
    return; //doesn't happen but makes compiler happy
  }
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->params.nth_fz=viewInfo.nth_fz;
    precisionRecord->params.nth_fz_limit.assign_across(viewInfo.nth_fz_limit);

    precisionRecord->orbit.evaluator.currentParams.mandel.c.assign_across(viewInfo.c);
    precisionRecord->params.c.assign(precisionRecord->orbit.evaluator.currentParams.mandel.c);
    recomputeRoot(viewInfo.max_root_effort);

    MandelMath::complex<MandelMath::number_any> old_center(&precisionRecord->orbit.evaluator.tmp);
    old_center.assign(precisionRecord->position.center);
    int old_step_log=precisionRecord->position.step_log;

    MandelMath::complex<MandelMath::number_any> view_here(&precisionRecord->orbit.evaluator.tmp);
    view_here.assign_across(viewInfo.view);
    precisionRecord->position.setView(view_here, viewInfo.scale);

    transformStore(precisionRecord->points, pointStore, 0, 0, old_center,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

class PixelPositionTransformer
{
public:
  PixelPositionTransformer(int lshift, int new_step_log);
  void setShift(MandelMath::number<MandelMath::number_any> *shift, int maxOut);
  int transform(int new_value, int offset);
protected:
  bool invalid;
  int new_step_log;
  int step_scale_n_shift, step_scale_d_shift, step_scale_d_mask;
  int delta_int;
};

/*PixelPositionTransformer::PixelPositionTransformer(int lshift, int new_step_log):
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

void PixelPositionTransformer::setShift(MandelMath::number<MandelMath::number_any> *shift, int maxOut)
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
}*/

void LaguerreModel::transformStore(void *old_points, LaguerrePointStore *old_store, int old_width, int old_height, MandelMath::complex<MandelMath::number_any> const &old_c,
                                   void *new_points, LaguerrePointStore *new_store, int new_width, int new_height, MandelMath::complex<MandelMath::number_any> const &new_c,
                                   int inlog, int new_step_log)
{
  if (precisionRecord==nullptr)
  {
    dbgPoint();
    return;
  };
  /*int indexOfWtiPoint, wtiIndexLen;
  precisionRecord->wtiPoint.self_allocator._getFirstCapac(indexOfWtiPoint, wtiIndexLen);
  working_assert(wtiIndexLen==LaguerrePoint<WORKER_MULTI>::LEN);*/
  //(oldx-width/2)*old_step+old_cre = (newx-width/2)*new_step+new_cre
  //oldx = (width/2) + (newx-width/2+(new_cre-old_cre)/new_step)*new_step/old_step
  PixelPositionTransformer ytrans(inlog, new_step_log);
  PixelPositionTransformer xtrans(inlog, new_step_log);
  {
    MandelMath::number<MandelMath::number_any> tmp(old_c.im);
    //tmp.assign(old_c.im);
    tmp.sub(new_c.im); //and reversing y at the last minute
    ytrans.setShift(&tmp, new_height);

    tmp.assign(new_c.re);
    tmp.sub(old_c.re);
    xtrans.setShift(&tmp, new_width);
  }

  MandelMath::complex<MandelMath::number_any> c(&precisionRecord->orbit.evaluator.tmp);
  for (int newy=0; newy<new_height; newy++)
  {
    int oldy=ytrans.transform(newy, old_height/2);
    //-1 is better than imageHeight, call reset() in the second loop, we may still need the points
    c.im.assign(precisionRecord->position.center.im);
    c.im.add_double((new_height/2-newy)*precisionRecord->position.step_size);
    for (int newx=0; newx<new_width; newx++)
    {
      int oldx=xtrans.transform(newx, old_width/2);
      //-1 is better than imageWidth, call reset() in the second loop, we may still need the points
      if ((oldy>newy) || ((oldy==newy) && (oldx>newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old point to new place
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].state==LaguerrePointStore::State::stWorking)
            new_store[newy*new_width+newx].state=LaguerrePointStore::State::stUnknown; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*LaguerrePoint<MandelMath::number_any>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*LaguerrePoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint<WORKER_MULTI>::LEN, old_sworker, (oldy*old_width+oldx)*LaguerrePoint<WORKER_MULTI>::LEN, LaguerrePoint<WORKER_MULTI>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          c.re.assign(precisionRecord->position.center.re);
          c.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(c);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*LaguerrePoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint<WORKER_MULTI>::LEN, precisionRecord->currentWorker, indexOfWtiPoint, LaguerrePoint<WORKER_MULTI>::LEN);
        }
      };
    }
  }
  for (int newy=(new_height-1)&0xfffffff; newy>=0; newy--) //avoid Clang warning about newy possibly ~ 2^31
  {
    int oldy=ytrans.transform(newy, old_height/2);
    c.im.assign(precisionRecord->position.center.im);
    c.im.add_double((new_height/2-newy)*precisionRecord->position.step_size);
    for (int newx=new_width-1; newx>=0; newx--)
    {
      int oldx=xtrans.transform(newx, old_width/2);
      if ((oldy<newy) || ((oldy==newy) && (oldx<=newx)))
      {
        if ((oldy>=0) && (oldy<old_height) && (oldx>=0) && (oldx<old_width))
        { //copy old to new
          new_store[newy*new_width+newx].assign(&old_store[oldy*old_width+oldx]);
          if (new_store[newy*new_width+newx].state==LaguerrePointStore::State::stWorking)
            new_store[newy*new_width+newx].state=LaguerrePointStore::State::stUnknown; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*LaguerrePoint<MandelMath::number_any>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*LaguerrePoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint<WORKER_MULTI>::LEN, old_sworker, (oldy*old_width+oldx)*LaguerrePoint<WORKER_MULTI>::LEN, LaguerrePoint<WORKER_MULTI>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          c.re.assign(precisionRecord->position.center.re);
          c.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(c);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*LaguerrePoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*LaguerrePoint<WORKER_MULTI>::LEN, precisionRecord->currentWorker, indexOfWtiPoint, LaguerrePoint<WORKER_MULTI>::LEN);
        }
      }
    }
  }
}

void LaguerreModel::setView(MandelMath::complex<MandelMath::number_any> const &c, double scale)
{
  MandelMath::complex<MandelMath::number_any> old_c(&precisionRecord->orbit.evaluator.tmp);
  old_c.assign(precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.setView(c, scale);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void LaguerreModel::drag(double delta_x, double delta_y)
{
  //qDebug()<<"drag ("<<delta_x<<","<<delta_y<<")";
  MandelMath::complex<MandelMath::number_any> old_c(&precisionRecord->orbit.evaluator.tmp);
  old_c.assign(precisionRecord->position.center);

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int dx=qRound(delta_x);
    int dy=qRound(delta_y);
    precisionRecord->position.move(dx, dy);
    //qDebug()<<"new c: re="<<position.worker->toString(&position.center_re_s)<<",im="<<position.worker->toString(&position.center_im_s);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void LaguerreModel::zoom(double x, double y, int inlog)
{
  MandelMath::complex<MandelMath::number_any> old_c(&precisionRecord->orbit.evaluator.tmp);

  old_c.assign(precisionRecord->position.center);
  int old_step_log=precisionRecord->position.step_log;

  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    precisionRecord->position.scale(inlog, qRound(x)-imageWidth/2, imageHeight/2-qRound(y));

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, old_c,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   precisionRecord->position.step_log-old_step_log, precisionRecord->position.step_log);

    startNewEpoch();
  }
}

void LaguerreModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int newLength=width*height;
    LaguerrePointStore *newStore=new LaguerrePointStore[newLength];
    {
      QString size_as_text=QString::number(newLength*sizeof(LaguerrePointStore));
      for (int pos=size_as_text.length()-3; pos>0; pos-=3)
        size_as_text.insert(pos, '\'');
      qDebug()<<"laguerreStore uses"<<size_as_text.toLocal8Bit().constData()<<"B";
    }
    void *new_points=nullptr;
    switch (precisionRecord->ntype)
    {
      case MandelMath::NumberType::typeEmpty: goto lolwut;
      case MandelMath::NumberType::typeDouble: lolwut:
        new_points=new double[newLength*LaguerrePoint<MandelMath::number_any>::LEN];
        break;
#if !NUMBER_DOUBLE_ONLY
      case MandelMath::NumberType::typeFloat128:
        new_points=new __float128[newLength*LaguerrePoint<MandelMath::number_any>::LEN];
        break;
      case MandelMath::NumberType::typeDDouble:
        new_points=new MandelMath::dd_real[newLength*LaguerrePoint<MandelMath::number_any>::LEN];
        break;
      case MandelMath::NumberType::typeQDouble:
        new_points=new MandelMath::dq_real[newLength*LaguerrePoint<MandelMath::number_any>::LEN];
        break;
      case MandelMath::NumberType::typeReal642:
        new_points=new MandelMath::real642[newLength*LaguerrePoint<MandelMath::number_any>::LEN];
        break;
#endif
    }
    MandelMath::complex<MandelMath::number_any> old_c(&precisionRecord->orbit.evaluator.tmp);
    old_c.assign(precisionRecord->position.center);

    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, old_c,
                   new_points, newStore, width, height, precisionRecord->position.center,
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

void LaguerreModel::pause(bool pause)
{
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    transformStore(precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   precisionRecord->points, pointStore, imageWidth, imageHeight, precisionRecord->position.center,
                   0, precisionRecord->position.step_log);

    if (!pause)
      startNewEpoch();
  }
}

void LaguerreModel::startNewEpoch()
{
  epoch=(epoch%2000000000)+1;
  nextGivenPointIndex=0;
  effortBonus=0;
#if 0
  for (int t=0; t<precisionRecord->threadCount; t++)
    if (precisionRecord->threads[t]->currentParams.pixelIndex<0)
      giveWork(precisionRecord->threads[t]);
#else
  MandelMath::visit([](auto &threads, int threadCount, int epoch){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
        for (int t=0; t<threadCount; t++)
        {
          //giveWorkToThread(precisionRecord->threads[t]);
          threads[t]->workIfEpoch=epoch;
        }
      };
  }, precisionRecord->threads, precisionRecord->threadCount, epoch);
  _threadsWorking+=precisionRecord->threadCount;
  emit triggerLaguerreThreaded(epoch, precisionRecord->params.period); //::invokeMethod cannot pass parameters, but ::connect can
#endif
}

void LaguerreModel::reimToPixel(int *circ_x, int *circ_y, MandelMath::complex<MandelMath::number_any> const &z, MandelMath::number<MandelMath::number_any> *tmp)
{
  double scaled;
  tmp->assign_across(z.re);
  tmp->sub(precisionRecord->position.center.re);
  scaled=tmp->toDouble()/precisionRecord->position.step_size;
  if (scaled<-10003 || scaled>10003)
    *circ_x=-100;
  else
    *circ_x=qRound(scaled)+imageWidth/2;

  tmp->assign_across(z.im);
  tmp->sub(precisionRecord->position.center.im);
  scaled=tmp->toDouble()/precisionRecord->position.step_size;
  if (scaled<-10003 || scaled>10003)
    *circ_y=-100;
  else
    *circ_y=imageHeight/2-qRound(scaled);
}

void LaguerreModel::paintOrbit(ShareableImageWrapper image, int x, int y)
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

  if (precisionRecord==nullptr)
  {
    dbgPoint(); //doesn't happen but makes compiler happy
    return;
  }
  /*else if (precisionRecord->currentWorker==nullptr)
    dbgPoint();
  else if (precisionRecord->orbit.currentWorker->ntype()!=precisionRecord->currentWorker->ntype())
    dbgPoint();*/
  if (precisionRecord->params.period<=0)
    return;

  //LaguerrePointStore *resultStore=&pointStore[y*imageWidth+x];
  //LaguerrePoint *data=&pointStore[y*imageWidth+x];
  /*switch (resultStore->state)
  {
    case LaguerrePoint::State::stUnknown:
    case LaguerrePoint::State::stResolved:
    case LaguerrePoint::State::stFail:
    default: ;
  }*/

  MandelMath::number<MandelMath::number_any> tmp(&precisionRecord->orbit.evaluator.tmp);
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint c
    reimToPixel(&circ_x, &circ_y, precisionRecord->params.c, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l3[3]={{circ_x+1, circ_y-1, circ_x-1, circ_y-1},
                   {circ_x-1, circ_y-1, circ_x-1, circ_y+1},
                   {circ_x-1, circ_y+1, circ_x+1, circ_y+1}};
      painter.drawLines(l3, 3);
    };

    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0, 0)); //paint root as +
    reimToPixel(&circ_x, &circ_y, precisionRecord->params.root, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };
  }


  //evaluator.params.c .. saved mouse coordinates
  //evaluator.data.root .. found root (temporary)
  //pointData.f .. not used
  //pointData.fz .. result f'
  precisionRecord->position.pixelXtoRE(x-imageWidth/2, &precisionRecord->orbit.evaluator.currentParams.mandel.first_z.re);
  precisionRecord->position.pixelYtoIM(imageHeight/2-y, &precisionRecord->orbit.evaluator.currentParams.mandel.first_z.im);
  precisionRecord->orbit.evaluator.currentParams.epoch=epoch;
  precisionRecord->orbit.evaluator.workIfEpoch=precisionRecord->orbit.evaluator.busyEpoch;//epoch;
  precisionRecord->orbit.evaluator.currentParams.pixelIndex=0;
  precisionRecord->orbit.evaluator.currentParams.nth_fz=precisionRecord->params.nth_fz;
  precisionRecord->orbit.evaluator.laguerreData.zero(precisionRecord->orbit.evaluator.currentParams.mandel.first_z);
  //precisionRecord->orbit.evaluator.currentData.store->wstate=LaguerrePointStore::State::stWorking;
  /*MandelMath::complex base(orbit.worker, &params.base_re_s_, &params.base_im_s, true);
  currentWorker->assign(&orbit.evaluator.currentData.root_re, &orbit.evaluator.currentParams.c_re);
  orbit.worker->assign(&orbit.evaluator.currentData.root_im, &orbit.evaluator.currentParams.c_im);
  MandelMath::complex root(orbit.worker, &orbit.evaluator.currentData.root_re, &orbit.evaluator.currentData.root_im, true);*/
  precisionRecord->orbit.evaluator.laguerreData.r.assign(precisionRecord->orbit.evaluator.currentParams.mandel.first_z);
  precisionRecord->orbit.evaluator.currentParams.mandel.c.assign_across(precisionRecord->params.c);
  {
    precisionRecord->orbit.evaluator.loope.eval_zz(1,
                                                   precisionRecord->orbit.evaluator.currentParams.mandel.c,
                                                   precisionRecord->orbit.evaluator.laguerreData.r, false, true);
    //precisionRecord->orbit.evaluator.bulb.bulbe.eval_zz(1, &precisionRecord->params.base, &precisionRecord->orbit.evaluator.currentData.root, false, true);
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0xff, 0xff)); //white X for next iteration
    reimToPixel(&circ_x, &circ_y, precisionRecord->orbit.evaluator.loope.f, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      //painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y-2, circ_x+2, circ_y+2},
                   {circ_x-2, circ_y+2, circ_x+2, circ_y-2}};
      painter.drawLines(l2, 2);
    };
    precisionRecord->orbit.evaluator.loope.eval_zz(precisionRecord->params.period-1,
                                                   precisionRecord->orbit.evaluator.currentParams.mandel.c,
                                                   precisionRecord->orbit.evaluator.laguerreData.r, false, false);
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0xff, 0xff)); //white (X) for next period
    reimToPixel(&circ_x, &circ_y, precisionRecord->orbit.evaluator.loope.f, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y-2, circ_x+2, circ_y+2},
                   {circ_x-2, circ_y+2, circ_x+2, circ_y-2}};
      painter.drawLines(l2, 2);
    };
  }
  precisionRecord->orbit.evaluator.thread.syncLaguerre(precisionRecord->params.period, &precisionRecord->orbit.evaluator.laguerreData.r, true);
  //orbit.pointData.assign(orbit.worker, orbit.evaluator.currentData);
  precisionRecord->orbit.pointDataStore.iter=precisionRecord->orbit.evaluator.newtres.cyclesNeeded;
  precisionRecord->orbit.pointDataStore.firstM=precisionRecord->orbit.evaluator.newtres.firstMum_re;
  precisionRecord->orbit.pointDataStore.firstStep_re=
      precisionRecord->orbit.evaluator.newtres.first_guess_lagu.re.toDouble()-
      precisionRecord->orbit.evaluator.currentParams.mandel.first_z.re.toDouble();
  precisionRecord->orbit.pointDataStore.firstStep_im=
      precisionRecord->orbit.evaluator.newtres.first_guess_lagu.im.toDouble()-
      precisionRecord->orbit.evaluator.currentParams.mandel.first_z.im.toDouble();
  precisionRecord->orbit.pointData.r.assign_across(precisionRecord->orbit.evaluator.laguerreData.r); //not used
  precisionRecord->orbit.pointData.fz_r.assign_across(precisionRecord->orbit.evaluator.newtres.fz_r);
  precisionRecord->orbit.pointData.fz_r.re.add_double(1);
  precisionRecord->orbit.pointData.nth_fz.assign_across(precisionRecord->orbit.evaluator.newtres.nth_fz);
  precisionRecord->orbit.pointDataStore.naiveChoice=precisionRecord->orbit.evaluator.newtres.naiveChoice;
  precisionRecord->orbit.first_mu_re=precisionRecord->orbit.evaluator.newtres.firstMu_re;
  precisionRecord->orbit.first_mu_im=precisionRecord->orbit.evaluator.newtres.firstMu_im;
  precisionRecord->orbit.first_mum_re=precisionRecord->orbit.evaluator.newtres.firstMum_re;
  precisionRecord->orbit.first_mum_im=precisionRecord->orbit.evaluator.newtres.firstMum_im;

  int circ_x, circ_y;
  double tmp_re, tmp_im;
  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0xff, 0));
  reimToPixel(&circ_x, &circ_y, precisionRecord->orbit.evaluator.laguerreData.r, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    double eps2=precisionRecord->position.center.re.eps2();
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3); //final=yellow /
    painter.drawLine(circ_x-2, circ_y+2, circ_x+2, circ_y-2);
    //size of r that maps to 1 epsilon: about the width of transition from >0 to <0
    int circ_r=ldexp(sqrt(precisionRecord->orbit.evaluator.newtres.accy_tostop*eps2), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
      painter.drawEllipse(circ_x-circ_r, circ_y-circ_r, 2*circ_r, 2*circ_r);
    //size of r that maps to a few epsilon: the accuracy of root and also size of the dead pool
    circ_r=ldexp(sqrt((precisionRecord->orbit.evaluator.newtres.accy_tostop*precisionRecord->orbit.evaluator.newtres.accy_multiplier)*eps2), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
      painter.drawEllipse(circ_x-circ_r, circ_y-circ_r, 2*circ_r, 2*circ_r);
    //new way to guess noise after iteration
    painter.setPen(QPen(QBrush(QColor(0xff, 0xff, 0)), 1, Qt::PenStyle::DashLine));
    circ_r=ldexp(sqrt((precisionRecord->orbit.evaluator.newtres.accy_noise)*eps2), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
      painter.drawEllipse(circ_x-circ_r, circ_y-circ_r, 2*circ_r, 2*circ_r);
    painter.setPen(QColor(0xff, 0xff, 0));
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0xff, 0, 0xff)); //Newton=purple \      .
  reimToPixel(&circ_x, &circ_y, precisionRecord->orbit.evaluator.newtres.first_guess_newt, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-2, circ_y-2, circ_x+2, circ_y+2);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //white = is better visible than cyan | so paint it under
  painter.setPen(QColor(0xff, 0xff, 0xff)); //Fejer=white =
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_fejer_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_fejer_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-3, circ_y-1, circ_x+3, circ_y-1);
    painter.drawLine(circ_x-3, circ_y+1, circ_x+3, circ_y+1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush); //
  painter.setPen(QColor(0xff, 0x00, 0x00)); //Naive=red o
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_naive1_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_naive1_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
  };
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_naive2_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_naive2_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_naive_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_naive_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre=cyan ||
  reimToPixel(&circ_x, &circ_y, precisionRecord->orbit.evaluator.newtres.first_guess_lagu, &tmp);
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x-1, circ_y-3, circ_x-1, circ_y+3);
    painter.drawLine(circ_x+1, circ_y-3, circ_x+1, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1=cyan |
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_lagu1_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_lagu1_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawLine(circ_x, circ_y-3, circ_x, circ_y+3);
  };

  painter.setBrush(Qt::BrushStyle::NoBrush);
  painter.setPen(QColor(0x80, 0xff, 0xff)); //Laguerre1other=cyan o
  tmp_re=precisionRecord->orbit.evaluator.newtres.first_lagu1o_re-precisionRecord->position.center.re.toDouble();
  tmp_re=ldexp(tmp_re, precisionRecord->position.step_log);
  circ_x=tmp_re+imageWidth/2;
  tmp_im=precisionRecord->orbit.evaluator.newtres.first_lagu1o_im-precisionRecord->position.center.im.toDouble();
  tmp_im=ldexp(tmp_im, precisionRecord->position.step_log);
  circ_y=imageHeight/2-tmp_im;
  if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
  {
    painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
    painter.drawEllipse(circ_x-1, circ_y-1, 2*1, 2*1);
  };


  painter.setBrush(Qt::BrushStyle::NoBrush);
#if 1 //Neumaier bounds
  if (precisionRecord->orbit.evaluator.newtres.first_neumaier1_im!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Neumaier1_im
    int circ_r=ldexp(abs(precisionRecord->orbit.evaluator.newtres.first_neumaier1_im), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  };
  if ((precisionRecord->orbit.evaluator.newtres.first_neumaier1_im!=0) || (precisionRecord->orbit.evaluator.newtres.first_neumaier1_re<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Neumaier1_mag
    double mag=sqrt(precisionRecord->orbit.evaluator.newtres.first_neumaier1_re*precisionRecord->orbit.evaluator.newtres.first_neumaier1_re+
                    precisionRecord->orbit.evaluator.newtres.first_neumaier1_im*precisionRecord->orbit.evaluator.newtres.first_neumaier1_im);
    int circ_r=ldexp(mag, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Neumaier1
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres.first_neumaier1_re, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }

  painter.setBrush(Qt::BrushStyle::NoBrush);
  if (precisionRecord->orbit.evaluator.newtres.first_neumaier2_im!=0)
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Neumaier2_im
    int circ_r=ldexp(abs(precisionRecord->orbit.evaluator.newtres.first_neumaier2_im), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  };
  if ((precisionRecord->orbit.evaluator.newtres.first_neumaier2_im!=0) || (precisionRecord->orbit.evaluator.newtres.first_neumaier2_re<0))
  {
    painter.setPen(QColor(0xff, 0x00, 0xff)); //Neumaier2_mag
    double mag=sqrt(precisionRecord->orbit.evaluator.newtres.first_neumaier2_re*precisionRecord->orbit.evaluator.newtres.first_neumaier2_re+
                    precisionRecord->orbit.evaluator.newtres.first_neumaier2_im*precisionRecord->orbit.evaluator.newtres.first_neumaier2_im);
    int circ_r=ldexp(mag, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
  else
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Neumaier2
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres.first_neumaier2_re, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
#else //Ostrowski bounds, don't look as good
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Ostrowski1c
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres_.ostrowski_r1c, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
    };
  }
  {
    painter.setPen(QColor(0xff, 0x00, 0x00)); //Ostrowski1x
    int circ_r=ldexp(abs(precisionRecord->orbit.evaluator.newtres_.ostrowski_r1x), precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      //painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
  {
    painter.setPen(QColor(0xff, 0xff, 0xff)); //Ostrowski2c
    int circ_r=ldexp(precisionRecord->orbit.evaluator.newtres_.ostrowski_r2c, precisionRecord->position.step_log);
    if ((circ_r>=1) && (circ_r<=10003))
    {
      painter.drawEllipse(x-circ_r, y-circ_r, 2*circ_r, 2*circ_r);
      painter.drawEllipse(x-circ_r-2, y-circ_r-2, 2*circ_r+4, 2*circ_r+4);
    };
  }
#endif


  image.image->swap(newOverlay);
}

int LaguerreModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return -1;
  timerWriteToImage.start();
  double params_nth_fz_limit=precisionRecord->params.nth_fz_limit.toDouble();
  LaguerrePointStore *wtiStore;
  static constexpr unsigned int FAIL_COLOR=0xff404040;
  for (int y=0; y<imageHeight; y++)
    for (int x=0; x<imageWidth; x++)
    {
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
      //LaguerrePoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiStore=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleCls:
        {
          switch (wtiStore->state.load())
          {
            case LaguerrePointStore::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              break;
            case LaguerrePointStore::State::stWorking:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xffff00ff);
              break;
            case LaguerrePointStore::State::stResolved:
            {
              if (wtiStore->iter<=1)
                image.image->setPixel(x, y, 0xff606060); //the dead pool
              else
              {
                int r;
                {
                  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*LaguerrePoint<MandelMath::number_any>::LEN);
                  double tr=precisionRecord->wtiPoint.fz_r.re.toDouble();
                  double ti=precisionRecord->wtiPoint.fz_r.im.toDouble();
                  if (tr*tr+ti*ti>1)
                    r=0xc0;
                  else
                    r=0x00;
                }
                int g;
                /*if (data->firstM>=3)
                  g=0xff;
                else if (data->firstM>=2)
                  g=0xbf+0x40*(data->firstM-2);
                else if (data->firstM>=1)
                  g=0x7f+0x40*(data->firstM-1);
                else if (data->firstM>=0)
                  g=0x3f+0x40*(data->firstM-0);
                else if (data->firstM>=-1)
                  g=0x00+0x40*(data->firstM+1);
                else
                  g=0xdf;*/
                if (wtiStore->firstM>=3)
                  g=0x7f;
                else if (wtiStore->firstM>=2)
                  g=0xbf-0x40*(wtiStore->firstM-2);
                else if (wtiStore->firstM>=1)
                  g=0xff-0x40*(wtiStore->firstM-1);
                else if (wtiStore->firstM>=0)
                  g=0x00+0x40*(wtiStore->firstM-0);
                else if (wtiStore->firstM>=-1)
                  g=0x40+0x40*(wtiStore->firstM+1);
                else
                  g=0x80;
                int b;
                switch (wtiStore->iter % 5)
                {
                  case  0: b=0x00; break;
                  case  1: b=0xff; break;
                  case  2: b=0xc0; break;
                  case  3: b=0x80; break;
                  case  4: b=0x40; break;
                  default: b=0xff;
                }
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case LaguerrePointStore::State::stFail:
            {
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, FAIL_COLOR);
            } break;
          }
        } break;
        case paintStyle::paintStyleNthFz:
        {
          switch (wtiStore->state.load())
          {
            case LaguerrePointStore::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              break;
            case LaguerrePointStore::State::stWorking:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xffff00ff);
              break;
            case LaguerrePointStore::State::stResolved:
            {
              if (wtiStore->iter<=1)
                image.image->setPixel(x, y, 0xff606060); //the dead pool
              else
              {
                precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*LaguerrePoint<MandelMath::number_any>::LEN);
                double fz_re=precisionRecord->wtiPoint.nth_fz.re.toDouble();
                double fz_im=precisionRecord->wtiPoint.nth_fz.im.toDouble();
                double fz_mag=fz_re*fz_re+fz_im*fz_im;
                double fz_abs=std::sqrt(fz_mag);
                int r;
                if (fz_mag>1)
                  r=params_nth_fz_limit>=1? 0xff : 0xc0;
                else if (fz_mag<params_nth_fz_limit)
                  r=0x60;
                else
                  r=0x00;
                int g=fz_abs==0?128:qRound(fz_re/fz_abs*127.49+128);
                int b=fz_abs==0?128:qRound(fz_im/fz_abs*127.49+256)%256;
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case LaguerrePointStore::State::stFail:
            {
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, FAIL_COLOR);
            } break;
          }
        } break;
        case paintStyle::paintStyleNthFz1:
        {
          switch (wtiStore->state.load())
          {
            case LaguerrePointStore::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              break;
            case LaguerrePointStore::State::stWorking:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xffff00ff);
              break;
            case LaguerrePointStore::State::stResolved:
            {
              if (wtiStore->iter<=1)
                image.image->setPixel(x, y, 0xff606060); //the dead pool
              else
              {
                precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*LaguerrePoint<MandelMath::number_any>::LEN);
                double fz_re=precisionRecord->wtiPoint.nth_fz.re.toDouble()-1;
                double fz_im=precisionRecord->wtiPoint.nth_fz.im.toDouble();
                double fz_mag=fz_re*fz_re+fz_im*fz_im;
                double fz_abs=std::sqrt(fz_mag);
                int r;
                if (fz_mag>1)
                  r=0xc0;//params_nth_fz_limit>=1? 0xff : 0xc0;
                else if (fz_mag<0.5)//params_nth_fz_limit)
                  r=0x60;
                else
                  r=0x00;
                int g=fz_abs==0?128:qRound(fz_re/fz_abs*127.49+128);
                int b=fz_abs==0?128:qRound(fz_im/fz_abs*127.49+256)%256;
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case LaguerrePointStore::State::stFail:
            {
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, FAIL_COLOR);
            } break;
          }
        } break;
        case paintStyle::paintStyleRoots:
        {
          switch (wtiStore->state.load())
          {
            case LaguerrePointStore::State::stUnknown:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xff000000);
              break;
            case LaguerrePointStore::State::stWorking:
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, 0xffff00ff);
              break;
            case LaguerrePointStore::State::stResolved:
            {
              if (wtiStore->iter<=1)
                image.image->setPixel(x, y, 0xff606060); //the dead pool
              else
              {
                precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*LaguerrePoint<MandelMath::number_any>::LEN);
                int r;
                {
                  double fz_re=precisionRecord->wtiPoint.fz_r.re.toDouble();
                  double fz_im=precisionRecord->wtiPoint.fz_r.im.toDouble();
                  double fz_mag=fz_re*fz_re+fz_im*fz_im;
                  //double fz_abs=std::sqrt(fz_mag);
                  if (fz_mag<1)
                    r=0x00;
                  else if (fz_mag>100)
                    r=0xff;
                  else
                    r=qRound(std::log(fz_mag)/std::log(100)*0x9f)+0x60;
                }
                int g, b;
                {
                  double step_re=wtiStore->firstStep_re;
                  double step_im=wtiStore->firstStep_im;
                  double step_abs=std::sqrt(step_re*step_re+step_im*step_im);
                  step_abs=step_abs==0?0:1.0/step_abs;
                  g=qRound(step_re*step_abs*127.49+128);
                  b=qRound(step_im*step_abs*127.49+256)%256;
                }
                image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
              }
            } break;
            case LaguerrePointStore::State::stFail:
            {
              //image.image->setPixel(x, y, 0xffffffff);
              image.image->setPixel(x, y, FAIL_COLOR);
            } break;
          }
        } break;
      }
    }
  //return result;
  qint64 elapsed=timerWriteToImage.elapsed();
  return elapsed;
}

void LaguerreModel::doneWorkInThread(MandelEvaluatorThread *)
{
  _threadsWorking--;
}

template <typename BASE>
int LaguerreModel::giveWorkThreaded(MandelEvaluator<BASE> *me)
{
  QReadLocker locker(&threading_mutex);
  me->timeInvokeSwitchTotal+=me->timeInvoke.nsecsElapsed();
  if (epoch!=me->busyEpoch)
    return 3;
  int retryEffortFrom=0;
  int nextEffortBonus=effortBonus; //don't jump to max instantly once pixels<threads
  while (retryEffortFrom>=0)
  {
    retryEffortFrom=-1;
    //MandelMath::complex<MandelMath::number_any> &tmpc=precisionRecord->params.c;//will assign to currentParams.c (position.worker, &params.base_re_s_, &params.base_im_s, true);
    MandelMath::complex<BASE> &root=me->laguerreData.r;
    for (int pi=0; pi<imageWidth*imageHeight; pi++)
    {
      int pointIndex=(nextGivenPointIndex+pi)%(imageWidth*imageHeight);
      //if ((lastGivenPointIndex_!=0) && (pointIndex==0))
        //dbgPoint();
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), pointIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
      LaguerrePointStore *storeAtIndex=&pointStore[pointIndex];
      if (storeAtIndex->state.load(std::memory_order_relaxed)!=LaguerrePointStore::State::stUnknown) //speedup/precheck
        continue;
      LaguerrePointStore::State state_prev=storeAtIndex->state.exchange(LaguerrePointStore::State::stWorking, std::memory_order_acquire);
      if (state_prev!=LaguerrePointStore::State::stUnknown)
      {
        //was working->do not touch, could be done by now
        //was done/fail->restore old state
        if (state_prev!=LaguerrePointStore::State::stWorking)
          storeAtIndex->state.store(state_prev, std::memory_order_relaxed);
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
          if (effort<8)
          {
            if (retryEffortFrom<0)
              retryEffortFrom=pointIndex;
            storeAtIndex->state.store(LaguerrePointStore::State::stUnknown, std::memory_order_release);
          }
          else
          {
            //if (me->currentWorker->ntype()!=precisionRecord->ntype)
              //dbgPoint();
            //evaluator->switchType(position.worker);
            //set up me->currentData.f
            precisionRecord->position.pixelXtoRE<BASE>(pointIndex%imageWidth - imageWidth/2, &me->currentParams.mandel.first_z.re);
            precisionRecord->position.pixelYtoIM<BASE>(imageHeight/2-pointIndex/imageWidth, &me->currentParams.mandel.first_z.im);
            root.assign(me->currentParams.mandel.first_z);
            //no can do me->currentData.f.assign(&root); //set both root and f: root will change, f needed in doneWork
            storeAtIndex->firstStep_re=root.re.toDouble();
            storeAtIndex->firstStep_im=root.im.toDouble();
            me->currentParams.epoch=me->busyEpoch;
            me->currentParams.pixelIndex=pointIndex;
            me->currentParams.nth_fz=precisionRecord->params.nth_fz;
            //me->startNewton(precisionRecord->params.period, &tmpc); //evaluator->currentData.f is additional hidden parameter
            //me->currentParams.store->period=precisionRecord->params.period;

            //pre-assigned in setParams() me->currentParams.c.assign_across(&tmpc);
            //me->currentWorker->assign_across(me->currentParams.c.re, precisionRecord->currentWorker, tmpc.re);
            //me->currentWorker->assign_across(me->currentParams.c.im, precisionRecord->currentWorker, tmpc.im);

            nextGivenPointIndex=(pointIndex+1)%(imageWidth*imageHeight);
            effortBonus=nextEffortBonus;
            return 0;
          }
        }
      }
      //MandelEvaluator::simple(cr, ci, pointStore[y*imageWidth+x]);
    }
    if ((retryEffortFrom>=0))// && (effortBonus<MAX_EFFORT))
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
int LaguerreModel::doneWorkThreaded(MandelEvaluator<BASE> *me, int result, bool giveWork)
{
  QReadLocker locker(&threading_mutex);
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), me->currentParams.pixelIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
    //LaguerrePoint point_(&pointStore_[me->currentParams.pixelIndex], &allo);
    /*
    if (point->state!=LaguerrePoint::State::stUnknown)
      qDebug()<<"Finished pixel finished again";
    else*/
    {
      /*if (me->currentWorker==nullptr)
        dbgPoint();
      else*/
      {
        LaguerrePointStore *dstStore=&pointStore[me->currentParams.pixelIndex];
        if (dstStore->state.load(std::memory_order_relaxed)==LaguerrePointStore::State::stUnknown)
          dbgPoint();
        if (dstStore->state.load(std::memory_order_relaxed)==LaguerrePointStore::State::stWorking)
        {
          if (result>0)
          {
            dstStore->state.store(LaguerrePointStore::State::stResolved, std::memory_order_release);
          }
          else
          {
            me->newtres.cyclesNeeded=-1;
            dstStore->state.store(LaguerrePointStore::State::stFail, std::memory_order_release);
          };
          //point->assign(position.worker, me->currentData);
          //in theory, LaguerrePoint can be assigned from NewtRes; in practice, they are quite different
          //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), me->currentParams.pixelIndex*LaguerrePoint::LEN, LaguerrePoint::LEN, nullptr);
          //LaguerrePoint point_(dstStore, &allo, nullptr); mhmm
          //me->laguerreData.r.assign(&me->laguerreData.r); //root
          me->laguerreData.fz_r.assign(me->newtres.fz_r);
          me->laguerreData.fz_r.re.add_double(1);
          me->laguerreData.nth_fz.assign(me->newtres.nth_fz);
          me->laguerreData.writeTo(precisionRecord->points, me->currentParams.pixelIndex*LaguerrePoint<MandelMath::number_any>::LEN);
          dstStore->iter=me->newtres.cyclesNeeded;
          dstStore->naiveChoice=me->newtres.naiveChoice;
          dstStore->firstM=me->newtres.firstMum_re;
          /* first lagu step - it's easier to just show direction to the root itself, not first step
          dstStore->firstStep_re=me->currentWorker->toDouble(me->newtres.first_guess_lagu.re)-
                                 me->currentWorker->toDouble(me->currentData.f.re);
          dstStore->firstStep_im=me->currentWorker->toDouble(me->newtres.first_guess_lagu.im)-
                                 me->currentWorker->toDouble(me->currentData.f.im);*/
          /* direction to the root - confusing boundaries because of jumps in multi
          dstStore->firstStep_re=me->currentWorker->toDouble(me->currentData.root.re)-
                                 me->currentWorker->toDouble(me->currentData.f.re);
          dstStore->firstStep_im=me->currentWorker->toDouble(me->currentData.root.im)-
                                 me->currentWorker->toDouble(me->currentData.f.im);*/
          dstStore->firstStep_re=me->newtres.first_lagu1_re-
                                 dstStore->firstStep_re;//me->currentData.f.re.toDouble();
          dstStore->firstStep_im=me->newtres.first_lagu1_im-
                                 dstStore->firstStep_im;//me->currentData.f.im.toDouble();
        }
        else
          dbgPoint();
      }
    }
  }
  else if (me->currentParams.epoch!=epoch)
  { }//qDebug()<<"Old pixel finished";
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

void LaguerreModel::selectedPrecisionChanged()
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
      new_points=MandelMath::number<double>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*LaguerrePoint<double>::LEN);
      break;
#if !NUMBER_DOUBLE_ONLY
    case precisionFloat128:
      newPrecision=Type::typeFloat128;
      new_points=MandelMath::number<__float128>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*LaguerrePoint<__float128>::LEN);
      break;
    case precisionDDouble:
      newPrecision=Type::typeDDouble;
      new_points=MandelMath::number<MandelMath::dd_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*LaguerrePoint<MandelMath::dd_real>::LEN);
      break;
    case precisionQDouble:
      newPrecision=Type::typeQDouble;
      new_points=MandelMath::number<MandelMath::dq_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*LaguerrePoint<MandelMath::dq_real>::LEN);
      break;
    case precisionReal642:
      newPrecision=Type::typeReal642;
      new_points=MandelMath::number<MandelMath::real642>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*LaguerrePoint<MandelMath::real642>::LEN);
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

LaguerreModel::Params::Params(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Params *source):
  period(source?source->period:1), nth_fz(source?source->nth_fz:1), c(spad), root(spad), nth_fz_limit(spad)
{
  if (source)
  {
    c.assign_across(source->c);
    root.assign_across(source->root);
    nth_fz_limit.assign_across(source->nth_fz_limit);
  }
  else
  {
    c.zero(0, 0);
    root.zero(0, 0);
    nth_fz_limit.zero(1);
  }
}


LaguerreModel::Position::Position(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Position *source):
  center(spad)
{
  if (source)
  {
    center.assign_across(source->center);
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

LaguerreModel::Position::~Position()
{
}

void LaguerreModel::Position::setView(MandelMath::complex<MandelMath::number_any> const &c, double scale)
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

void LaguerreModel::Position::move(int delta_x, int delta_y)
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

void LaguerreModel::Position::scale(int inlog, int center_x, int center_y)
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
    center.re.add_double((center_x-adjust_x)*(old_step_size-step_size)); //(old_step_size-step_size__)=(1-(1<<-inlog))*old_step_size
    center.im.add_double((center_y-adjust_y)*(old_step_size-step_size));
    //need to roll in high bits anyway
    /*cached_center_re_mod=cached_center_re_mod*(old_step_size/step_size)+(center_x-adjust_y)*(old_step_size/step_size-1);
    cached_center_re_mod%=32768;
    cached_center_im_mod=cached_center_im_mod*(old_step_size/step_size)+(center_y-adjust_y)*(old_step_size/step_size-1);
    cached_center_im_mod%=32768;*/
    updateCachedDepth();
  }
}

void LaguerreModel::Position::updateCachedDepth()
{
  MandelMath::complex<MandelMath::number_any> d(center);

  //d.assign(center);
  d.lshift(step_log-15);
  d.re.mod1();
  d.im.mod1();
  d.lshift(15);
  cached_center_re_mod=d.re.toRound();
  cached_center_im_mod=d.im.toRound();
}

template <typename BASE>
void LaguerreModel::Position::pixelXtoRE(int x, MandelMath::number<BASE> *result)
{
  result->assign_across(center.re);
  result->add_double(x*step_size);
}

template <typename BASE>
void LaguerreModel::Position::pixelYtoIM(int y, MandelMath::number<BASE> *result)
{
  result->assign_across(center.im);
  result->add_double(y*step_size);
}



LaguerreModel::Orbit::Orbit(MandelMath::NumberType ntype):
  evaluator(ntype, true),
  pointDataStore(), pointData(&pointDataStore, &evaluator.tmp),
  first_mu_re(0), first_mu_im(0), first_mum_re(0), first_mum_im(0)
{
}

LaguerreModel::Orbit::~Orbit()
{
  evaluator.workIfEpoch=-1;
  evaluator.thread.quit();
  evaluator.thread.wait(1000);
}


LaguerreModel::PrecisionRecord::PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, LaguerreModel *doneReceiver):
  ntype(ntype), orbit(ntype), wtiPoint(nullptr, &orbit.evaluator.tmp),
  params(&orbit.evaluator.tmp, source?&source->params:nullptr),
  position(&orbit.evaluator.tmp, source?&source->position:nullptr),
  tmp_place(&orbit.evaluator.tmp),
  threadCount(source?source->threadCount:0), threads()
{
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
    threads=new MandelEvaluator<double> *[threadCount];
    break;
//somehow...
#if !NUMBER_DOUBLE_ONLY
  case MandelMath::NumberType::typeFloat128:
    threads=new MandelEvaluator<__float128> *[threadCount];
    break;
  case MandelMath::NumberType::typeDDouble:
    threads=new MandelEvaluator<MandelMath::dd_real> *[threadCount];
    break;
  case MandelMath::NumberType::typeQDouble:
    threads=new MandelEvaluator<MandelMath::dq_real> *[threadCount];
    break;
  case MandelMath::NumberType::typeReal642:
    threads=new MandelEvaluator<MandelMath::real642> *[threadCount];
    break;
#endif
  }
  MandelMath::visit([](auto &threads, PrecisionRecord &self, LaguerreModel *doneReceiver, PrecisionRecord *source){
      for (int t=0; t<self.threadCount; t++)
      {
        if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
        {
          using Evaluator=std::remove_pointer_t<std::remove_pointer_t<std::remove_reference_t<decltype(threads)>>>;
          Evaluator *thread=new Evaluator(self.ntype, source==nullptr);
          threads[t]=thread;
          thread->threaded.give=[doneReceiver](Evaluator *me)
          {
            return doneReceiver->giveWorkThreaded(me);
          };
          thread->threaded.doneLaguerre=[doneReceiver](Evaluator *me, int result, bool giveWork)
          {
            return doneReceiver->doneWorkThreaded(me, result, giveWork);
          };
          QObject::connect(&thread->thread, &MandelEvaluatorThread::doneLaguerreThreaded,
                           doneReceiver, &LaguerreModel::doneWorkInThread,
                           Qt::ConnectionType::QueuedConnection);
          QObject::connect(doneReceiver, &LaguerreModel::triggerLaguerreThreaded,
                           &thread->thread, &MandelEvaluatorThread::doLaguerreThreaded,
                           Qt::ConnectionType::QueuedConnection);
        }
      }
  }, threads, *this, doneReceiver, source);
}

LaguerreModel::PrecisionRecord::~PrecisionRecord()
{
  MandelMath::visit([](auto &threads, int threadCount){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
          for (int t=threadCount-1; t>=0; t--)
          {
              threads[t]->workIfEpoch=-1;
              threads[t]->thread.quit();
          }
      };
  }, threads, threadCount);
  MandelMath::visit([](auto &threads, int threadCount){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
          for (int t=threadCount-1; t>=0; t--)
          {
              threads[t]->thread.wait(1000);
              delete threads[t];
          }
          delete[] threads;
      };
  }, threads, threadCount);
  threadCount=0;
  threads=nullptr;
}
