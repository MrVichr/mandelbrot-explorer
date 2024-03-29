#include <math.h>
#include <QObject>
#include <QDebug>
#include <QPainter>

#include "JuliaModel.hpp"
#include "MandelEvaluator.hpp"

#define working_assert(x) { if (!(x)) dbgPoint(); }

#define UPDATE_CACHED_MOD 0 //works until precision doesn't allow

JuliaModel::JuliaModel(): QObject(),
  _extAngleZoom(0), _selectedEdePatchAlgo(ext_de_patch_algo::EedepaBlend), _selectedExteriorColoring(exterior_coloring::EextcolRainbow),
  pointStore(nullptr), precisionRecord(nullptr)
{
  unsigned int oldcw; //524319 = 0x8001F = mask all interrupts, 80bit precision
  MandelMath::fpu_fix_start(&oldcw);
  _selectedPaintStyle=paintStyleCls;//Kind;
  _selectedPrecision=EprecisionDouble;
  epoch=1;
  imageWidth=0;
  imageHeight=0;
  //pointStore=nullptr;
  nextGivenPointIndex=0;
  effortBonus=0;
  //threadCount=4;
  _threadsWorking=0;
  QObject::connect(this, &JuliaModel::selectedPrecisionChange,
                   this, &JuliaModel::selectedPrecisionChanged);
  selectedPrecisionChanged();
  QObject::connect(this, &JuliaModel::selectedEdePatchAlgoChange,
                   this, &JuliaModel::selectedEdePatchAlgoChanged);
  selectedEdePatchAlgoChanged();
  QObject::connect(this, &JuliaModel::selectedExteriorColoringChange,
                   this, &JuliaModel::selectedExteriorColoringChanged);
  selectedExteriorColoringChanged();
}

JuliaModel::~JuliaModel()
{
  epoch=(epoch%2000000000)+1;
  delete precisionRecord;

  delete[] pointStore;
  pointStore=nullptr;
  imageWidth=0;
  imageHeight=0;
}

void JuliaModel::startRunning()
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

QString JuliaModel::pixelXtoRE_str(int x)
{
  MandelMath::number<MandelMath::number_any> num(precisionRecord->position.center.re);
  num.add_double((x - imageWidth/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString JuliaModel::pixelYtoIM_str(int y)
{
  MandelMath::number<MandelMath::number_any> num(precisionRecord->position.center.im);
  num.add_double((y - imageHeight/2)*precisionRecord->position.step_size);
  QString result=num.toString();
  return result;
}

QString JuliaModel::getTimes()
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

QString JuliaModel::getTextXY()
{
  if (precisionRecord==nullptr)
    return "-";
  return precisionRecord->orbit.evaluator.currentParams.julia.first_z.re.toString()+" +i* "+
         precisionRecord->orbit.evaluator.currentParams.julia.first_z.im.toString();
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

QString JuliaModel::getTextInfoGen()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_any> tmp(&precisionRecord->orbit.evaluator.tmp);
    reimToPixel(&orbit_x, &orbit_y, precisionRecord->orbit.evaluator.currentParams.julia.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*JuliaPoint::LEN, JuliaPoint::LEN, nullptr);
  //JuliaPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  JuliaPointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;
  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*JuliaPoint<MandelMath::number_any>::LEN);
  JuliaPoint<MandelMath::number_any> *data=&precisionRecord->wtiPoint;

  QString state;
  QString fc_fz;
  switch (data_store->rstate)
  {
    case JuliaPointStore::ResultState::stUnknown:
      if (data_store->wstate.load()==JuliaPointStore::WorkState::stIdle)
        state="Unk";
      else
        state="Working...";
      fc_fz=" fz="+data->fz_z.toString();
      break;
    case JuliaPointStore::ResultState::stOutside:
      state="Out"; fc_fz=" fz="+data->fz_z.toString(); break;
    case JuliaPointStore::ResultState::stOutAngle:
      state="OutA"; fc_fz=" fz="+data->fz_z.toString(); break;
    case JuliaPointStore::ResultState::stBoundary:
      state="Bound"; fc_fz=" fz="+data->fz_z.toString(); break;
    case JuliaPointStore::ResultState::stDiverge:
      state="Diver"; fc_fz=" fz="+data->fz_z.toString(); break;
    case JuliaPointStore::ResultState::stMisiur:
      state="Misiur"; fc_fz=" fz="+data->fz_z.toString(); break;
    case JuliaPointStore::ResultState::stPeriod2:
      state="Per2";
      switch (_selectedPaintStyle)
      {
        case paintStylePhi: fc_fz=" fi="+QString::number(data_store->interior.phi_re)+"+i*"+QString::number(data_store->interior.phi_im); break;
        case paintStylePhi1: fc_fz=" fi'="+QString::number(data_store->interior.phi1_re)+"+i*"+QString::number(data_store->interior.phi1_im); break;
        case paintStylePhi2: fc_fz=" fi''="+QString::number(data_store->interior.phi2_re)+"+i*"+QString::number(data_store->interior.phi2_im); break;
        case paintStyleF1: fc_fz=" f'="+QString::number(data_store->interior.f1_re)+"+i*"+QString::number(data_store->interior.f1_im); break;
        case paintStyleFirstUnder1: fc_fz=" first="+QString::number(data_store->interior.first_under_1); break;
        default: fc_fz=" fz="+data->fz_z.toString();
      }
      break;
    case JuliaPointStore::ResultState::stPeriod3:
      state="Per3";
      switch (_selectedPaintStyle)
      {
        case paintStylePhi: fc_fz=" fi="+QString::number(data_store->interior.phi_re)+"+i*"+QString::number(data_store->interior.phi_im); break;
        case paintStylePhi1: fc_fz=" fi'="+QString::number(data_store->interior.phi1_re)+"+i*"+QString::number(data_store->interior.phi1_im); break;
        case paintStylePhi2: fc_fz=" fi''="+QString::number(data_store->interior.phi2_re)+"+i*"+QString::number(data_store->interior.phi2_im); break;
        case paintStyleF1: fc_fz=" f'="+QString::number(data_store->interior.f1_re)+"+i*"+QString::number(data_store->interior.f1_im); break;
        case paintStyleFirstUnder1: fc_fz=" first="+QString::number(data_store->interior.first_under_1); break;
        default: fc_fz=" fz="+data->fz_z.toString();
      }
      break;
    case JuliaPointStore::ResultState::stMaxIter:
      state="Max"; fc_fz=" fz="+data->fz_z.toString(); break;
  }

  if (precisionRecord->params.period<=0)
    return state+" iter="+QString::number(data_store->iter)+
                 //" nearz="+QString::number(data_store->nearziter_0)+"("+QString::number(data_store->nearziter_0)+")"+
                 //" bigfz="+QString::number(data_store->bigfzm_iter)+"("+QString::number(data_store->bigfzm_iter)+")"+
                 " bigfzfzz="+QString::number(data_store->bigfzfzziter)+"("+QString::number(data_store->bigfzfzziter)+")"+
                 " near0="+QString::number(data_store->near0iter_1)+"("+QString::number(data_store->near0iter_1)+")"+
                 " ("+QString::number((data_store->near0iter_1-data_store->nearziter_0))+")"+
           fc_fz;
  else
    return state+" iter="+QString::number(data_store->iter)+
                 //" nearz="+QString::number(data_store->nearziter_0)+"("+QString::number(data_store->nearziter_0%precisionRecord->params.period)+")"+
                 //" bigfz="+QString::number(data_store->bigfzm_iter)+"("+QString::number(data_store->bigfzm_iter%precisionRecord->params.period)+")"+
                 " bigfzfzz="+QString::number(data_store->bigfzfzziter)+"("+QString::number(data_store->bigfzfzziter%precisionRecord->params.period)+")"+
                 " near0="+QString::number(data_store->near0iter_1)+"("+QString::number(data_store->near0iter_1%precisionRecord->params.period)+")"+
                 " ("+QString::number((data_store->near0iter_1-data_store->nearziter_0))+")"+
           fc_fz;
}

QString JuliaModel::getTextInfoSpec()
{
  if (precisionRecord==nullptr)
    return "-";
  int orbit_x, orbit_y;
  {
    MandelMath::number<MandelMath::number_any> tmp(&precisionRecord->orbit.evaluator.tmp);
    reimToPixel(&orbit_x, &orbit_y, precisionRecord->orbit.evaluator.currentParams.julia.first_z, &tmp);
  }
  if ((orbit_x<0) || (orbit_x>=imageWidth) || (orbit_y<0) | (orbit_y>=imageHeight))
    return "? +i* ?";
  //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (orbit_x+imageWidth*orbit_y)*JuliaPoint::LEN, JuliaPoint::LEN, nullptr);
  //JuliaPoint data_(&pointStore_[orbit_x+imageWidth*orbit_y], &allo);
  JuliaPointStore *data_store=&pointStore[orbit_x+imageWidth*orbit_y];
  precisionRecord->wtiPoint.store=data_store;

  precisionRecord->wtiPoint.readFrom(precisionRecord->points, (orbit_x+imageWidth*orbit_y)*JuliaPoint<MandelMath::number_any>::LEN);
  //JuliaPoint<MandelMath::number_any> *data=&precisionRecord->wtiPoint;

  double exterior;
  switch (_selectedEdePatchAlgo)
  {
    case EedepaAlways: exterior=data_store->exterior.always_mixed; break;
    case EedepaConditional: exterior=data_store->exterior.condi_mixed; break;
    case EedepaBlend: exterior=data_store->exterior.blend_mixed; break;
    case EedepaVerifyH: exterior=data_store->exterior.verify_mixed/data_store->exterior.none_mixed-1; break;
    case EedepaVerifyA: exterior=data_store->exterior.verify_empty/data_store->exterior.none_empty-1; break;
    default: exterior=data_store->exterior.none_mixed; break;
  }
  switch (data_store->rstate)
  {
    case JuliaPointStore::ResultState::stUnknown:
      if (data_store->wstate.load()==JuliaPointStore::WorkState::stIdle)
        return QString("per=")+QString::number(precisionRecord->params.period)+" ";
      else
        return "Working...";
      break;
    case JuliaPointStore::ResultState::stOutside:
    {
      return QString("per=")+QString::number(precisionRecord->params.period)+
             QString(" ext=")+QString::number(exterior)+
             " near-r="+QString::number(data_store->nearmriter)+(data_store->nearmr_difficult?" !!":"");
    } break;
    case JuliaPointStore::ResultState::stOutAngle:
    {
      return QString("per=")+QString::number(precisionRecord->params.period)+
             QString(" ext=")+QString::number(exterior)+
             " near-r="+QString::number(data_store->nearmriter)+(data_store->nearmr_difficult?" !!":"")+
             " phi="+precisionRecord->wtiPoint.extangle.toString();
    } break;
    case JuliaPointStore::ResultState::stBoundary:
      return QString("per=")+QString::number(precisionRecord->params.period);
    case JuliaPointStore::ResultState::stDiverge:
      return QString("per=")+QString::number(precisionRecord->params.period);
    case JuliaPointStore::ResultState::stMisiur:
      return QString("per=")+QString::number(precisionRecord->params.period);
    case JuliaPointStore::ResultState::stPeriod2:
    {
      //double p=std::atan2((precisionRecord->orbit.bulb.baseFz_.im.toDouble()),
      //                     precisionRecord->orbit.bulb.baseFz_.re.toDouble())*precisionRecord->orbit.bulb.foundMult_;
      return QString("per=")+QString::number(precisionRecord->params.period)+
                    " int="+QString::number(data_store->interior.mixed)+
                    " near-r="+QString::number(data_store->nearmriter)+(data_store->nearmr_difficult?" !!":"");
        //  +" mult="+QString::number(std::round(p/(2*M_PI)))+"/"+QString::number(precisionRecord->orbit.bulb.foundMult_);
    } break;
    case JuliaPointStore::ResultState::stPeriod3:
    {
      //double p=std::atan2((precisionRecord->orbit.bulb.baseFz_.im.toDouble()),
      //                     precisionRecord->orbit.bulb.baseFz_.re.toDouble())*precisionRecord->orbit.bulb.foundMult_;
      return QString("per=")+QString::number(precisionRecord->params.period)+
                     " int="+QString::number(data_store->interior.mixed)+
                     " near-r="+QString::number(data_store->nearmriter)+(data_store->nearmr_difficult?" !!":"");
        //  +" mult="+QString::number(std::round(p/(2*M_PI)))+"/"+QString::number(precisionRecord->orbit.bulb.foundMult_);
    } break;
    case JuliaPointStore::ResultState::stMaxIter:
      return QString("per=")+QString::number(precisionRecord->params.period);
  }
  return "-?-?-";
}

void JuliaModel::recomputeRoot(int max_effort)
{
  MandelEvaluator<MandelMath::number_any> &evaluator=precisionRecord->orbit.evaluator;
  evaluator.currentParams.epoch=epoch;
  evaluator.workIfEpoch=evaluator.busyEpoch;//epoch;
  evaluator.currentParams.pixelIndex=0;
  evaluator.currentParams.mandel.c.assign(precisionRecord->params.c);//evaluator.currentParams.julia.c);
  evaluator.currentParams.mandel.first_z.assign(precisionRecord->params.c);//&evaluator.currentParams.julia.c);
  evaluator.mandelData.zero(evaluator.currentParams.mandel.first_z);
  evaluator.mandelData.store->wstate=MandelPointStore::WorkState::stWorking;
  evaluator.currentParams.breakOnNewNearest=false;
  evaluator.currentParams.maxiter=1<<max_effort;
  evaluator.currentParams.want_extangle=false;

  evaluator.thread.syncMandel();
  if (evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stPeriod2 ||
      evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stPeriod3)
  {
    precisionRecord->params.period=evaluator.mandelData.store->period;
    evaluator.currentParams.julia.period=evaluator.mandelData.store->period;
    evaluator.currentParams.julia.root.assign(evaluator.mandelData.root);
    precisionRecord->params.root.assign(evaluator.mandelData.root);
    evaluator.loope.eval_zz(precisionRecord->params.period-1, evaluator.currentParams.julia.c, evaluator.currentParams.julia.root, false, true);
    evaluator.currentParams.julia.root0.assign(evaluator.loope.f);
    evaluator.currentParams.julia.alphaShort.assign(evaluator.loope.f_z);
    evaluator.loope.eval_zz(1, evaluator.currentParams.julia.c, evaluator.currentParams.julia.root, false, false);
    evaluator.currentParams.julia.alpha.assign(evaluator.loope.f_z);
    //also jas*2*r0=alpha -> jas=alpha/(2*r0)=alpha/(2*sqrt(r-c))
    /* looks the same, up to precision, except has division by 0 and has ambiguity of +-sqrt
    evaluator.bulb.t1.assign(&evaluator.currentParams.juliaRoot);
    evaluator.bulb.t1.sub(&evaluator.currentParams.c);
    evaluator.bulb.t1.sqrt(&evaluator.tmp);
    evaluator.bulb.t1.lshift(1);
    evaluator.bulb.t1.recip(&evaluator.tmp);
    evaluator.bulb.t1.mul(&evaluator.currentParams.juliaAlpha, &evaluator.tmp);
    */
  }
  else if (precisionRecord->orbit.evaluator.mandelData.store->rstate==MandelPointStore::ResultState::stMaxIter)
  {
    precisionRecord->params.period=evaluator.mandelData.store->near0iter_1;
    evaluator.currentParams.julia.period=evaluator.mandelData.store->near0iter_1;
    evaluator.currentParams.julia.root.assign(evaluator.mandelData.f);
    precisionRecord->params.root.assign(evaluator.mandelData.f);
    evaluator.loope.eval_zz(precisionRecord->params.period-1, evaluator.currentParams.julia.c, evaluator.currentParams.julia.root, false, true);
    evaluator.currentParams.julia.root0.assign(evaluator.mandelData.f);
    evaluator.currentParams.julia.alphaShort.assign(evaluator.loope.f_z);
    evaluator.loope.eval_zz(1, evaluator.currentParams.julia.c, evaluator.currentParams.julia.root, false, false);
    evaluator.currentParams.julia.alpha.assign(evaluator.loope.f_z);
    //precisionRecord->orbit.evaluator.currentParams.juliaAlpha.recip(&precisionRecord->orbit.evaluator.tmp);
  }
  else
  {
    precisionRecord->params.period=0;
    evaluator.currentParams.julia.period=0;
    precisionRecord->params.root.zero(0, 0);
    evaluator.currentParams.julia.root.zero(0, 0);
    evaluator.currentParams.julia.root0.zero(0, 0);
    evaluator.currentParams.julia.alphaShort.zero(2, 0);
    evaluator.currentParams.julia.alpha.zero(2, 0);
  }

  //find distance estimate at c (interested in exterior for now)
  evaluator.currentParams.julia.patchSizeExterior=0; //already used by syncJulia()
  evaluator.currentParams.julia.first_z.assign(evaluator.currentParams.julia.c);
  evaluator.juliaData.zero(evaluator.currentParams.julia.first_z);
  evaluator.juliaData.store->wstate=JuliaPointStore::WorkState::stWorking;
  evaluator.thread.syncJulia();
  if (MandelMath::enum_is_one_of<JuliaPointStore::ResultState,
      JuliaPointStore::ResultState::stOutside,
      JuliaPointStore::ResultState::stOutAngle>(evaluator.juliaData.store->rstate))
  {
    //will be comparing .near0 to sqrt(d.e.) is same as .near0m to d.e.
    if (evaluator.juliaData.store->exterior.none_empty>=0)
    {
      // for what x does r=d/(2*x) differ from r^2+2xr-d=0 by 10% ?
      // r1=d/(x+sqrt(x^2+d))  r2=d/(2x)    r1/r2=d/(x+sqrt(x^2+d))/d*2*x=2/(1+sqrt(1+d/x^2))
      // k=2/(1+sqrt(1+d/x^2))   x=sqrt(d)/sqrt((2/k-1)^2-1))
      // k=0.828 <- r=sqrt(d); k=0.9 -> 1.423*sqrt(d); k=0.95->2.124*sqrt(d)   k=0.899 <- sqrt(2*d)
      evaluator.currentParams.julia.patchSizeExterior=evaluator.juliaData.store->exterior.none_empty   *2;
    }
  };

  // z^2+c=z*e, for which |z| is |e|>1 ?
  // for which z is |e|=1?
  // z^2-e*z+c=0
  // (e+-sqrt(e^2-4*c))/2=z1[e],z2[e]
  // |(1+-sqrt(1-4*c/e^2))|/2=|z1|[e],|z2|[e]
  // worst case |z|=(1+sqrt(1+4*|c|))/2 so everything above (and some below) flies even farther
  //evaluator.currentParams.juliaBailout=MandelMath::sqr_double(1+std::sqrt(1+4*std::sqrt(evaluator.currentParams.c.getMag_double())))/4;
  evaluator.currentParams.julia.bailout.assign(*evaluator.currentParams.julia.c.getMag_tmp());
  evaluator.currentParams.julia.bailout.sqrt();
  evaluator.currentParams.julia.bailout.lshift(2);
  evaluator.currentParams.julia.bailout.add_double(1);
  evaluator.currentParams.julia.bailout.sqrt();
  evaluator.currentParams.julia.bailout.add_double(1);
  evaluator.currentParams.julia.bailout.lshift(-1);
  evaluator.currentParams.julia.bailout.sqr(); // ( (1+sqrt(1+4*|c|))/2 )^2



  precisionRecord->params.root.assign(evaluator.currentParams.julia.root);

  MandelMath::visit([](auto &threads, int threadCount, MandelEvaluator<MandelMath::number_any> &evaluator){
      if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
      {
        //using Evaluator=std::remove_pointer_t<std::remove_pointer_t<std::remove_reference_t<decltype(threads)>>>;
        for (int t=0; t<threadCount; t++)
        {
          threads[t]->currentParams.julia.assign_across(evaluator.currentParams.julia);
        }
      };
  }, precisionRecord->threads, precisionRecord->threadCount, evaluator);
}

void JuliaModel::setParams(ShareableViewInfo viewInfo)
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

    precisionRecord->orbit.evaluator.currentParams.julia.c.assign_across(viewInfo.c);
    precisionRecord->params.c.assign(precisionRecord->orbit.evaluator.currentParams.julia.c);
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

void JuliaModel::transformStore(void *old_points, JuliaPointStore *old_store, int old_width, int old_height, MandelMath::complex<MandelMath::number_any> const &old_c,
                                void *new_points, JuliaPointStore *new_store, int new_width, int new_height, MandelMath::complex<MandelMath::number_any> const &new_c,
                                int inlog, int new_step_log)
{
  if (precisionRecord==nullptr)
  {
    dbgPoint();
    return;
  };
  /*int indexOfWtiPoint, wtiIndexLen;
  precisionRecord->wtiPoint.self_allocator._getFirstCapac(indexOfWtiPoint, wtiIndexLen);
  working_assert(wtiIndexLen==JuliaPoint<WORKER_MULTI>::LEN);*/
  //(oldx-width/2)*old_step+old_cre = (newx-width/2)*new_step+new_cre
  //oldx = (width/2) + (newx-width/2+(new_cre-old_cre)/new_step)*new_step/old_step
  PixelPositionTransformer ytrans(inlog, new_step_log);
  PixelPositionTransformer xtrans(inlog, new_step_log);
  {
    MandelMath::number<MandelMath::number_any> tmp(old_c.im);
    tmp.sub(new_c.im); //and reversing y at the last minute
    ytrans.setShift(&tmp, new_height);

    tmp.assign(new_c.re);
    tmp.sub(old_c.re);
    xtrans.setShift(&tmp, new_width);
  }

  MandelMath::complex<MandelMath::number_any> first_z(&precisionRecord->orbit.evaluator.tmp);
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
          if (new_store[newy*new_width+newx].wstate==JuliaPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=JuliaPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*JuliaPoint<MandelMath::number_any>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*JuliaPoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*JuliaPoint<WORKER_MULTI>::LEN, old_sworker, (oldy*old_width+oldx)*JuliaPoint<WORKER_MULTI>::LEN, JuliaPoint<WORKER_MULTI>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          first_z.re.assign(precisionRecord->position.center.re);
          first_z.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(first_z);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*JuliaPoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*JuliaPoint<WORKER_MULTI>::LEN, precisionRecord->currentWorker, indexOfWtiPoint, JuliaPoint<WORKER_MULTI>::LEN);
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
          if (new_store[newy*new_width+newx].wstate==JuliaPointStore::WorkState::stWorking)
            new_store[newy*new_width+newx].wstate=JuliaPointStore::WorkState::stIdle; //work will be cancelled because of new epoch
          precisionRecord->wtiPoint.readFrom(old_points, (oldy*old_width+oldx)*JuliaPoint<MandelMath::number_any>::LEN);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*JuliaPoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*JuliaPoint<WORKER_MULTI>::LEN, old_sworker, (oldy*old_width+oldx)*JuliaPoint<WORKER_MULTI>::LEN, JuliaPoint<WORKER_MULTI>::LEN);
        }
        else
        { //make fresh point in wti and copy to new
          first_z.re.assign(precisionRecord->position.center.re);
          first_z.re.add_double((newx - new_width/2)*precisionRecord->position.step_size);
          precisionRecord->wtiPoint.store=&new_store[newy*new_width+newx];
          precisionRecord->wtiPoint.zero(first_z);
          precisionRecord->wtiPoint.writeTo(new_points, (newy*new_width+newx)*JuliaPoint<MandelMath::number_any>::LEN);
          //new_sworker->assign_block((newy*new_width+newx)*JuliaPoint<WORKER_MULTI>::LEN, precisionRecord->currentWorker, indexOfWtiPoint, JuliaPoint<WORKER_MULTI>::LEN);
        }
      }
    }
  }
}

void JuliaModel::setView(MandelMath::complex<MandelMath::number_any> const &c, double scale)
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

void JuliaModel::drag(double delta_x, double delta_y)
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

void JuliaModel::zoom(double x, double y, int inlog)
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

void JuliaModel::setImageSize(int width, int height)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  if ((width==imageWidth) && (height==imageHeight))
    return;
  {
    QWriteLocker locker(&threading_mutex);
    epoch=(epoch%2000000000)+1; //invalidate threads while transforming store
    int newLength=width*height;
    JuliaPointStore *newStore=new JuliaPointStore[newLength];
    {
      QString size_as_text=QString::number(newLength*sizeof(JuliaPointStore));
      for (int pos=size_as_text.length()-3; pos>0; pos-=3)
        size_as_text.insert(pos, '\'');
      qDebug()<<"laguerreStore uses"<<size_as_text.toLocal8Bit().constData()<<"B";
    }
    void *new_points=nullptr;
    switch (precisionRecord->ntype)
    {
      case MandelMath::NumberType::typeEmpty: goto lolwut;
      case MandelMath::NumberType::typeDouble: lolwut:
        new_points=MandelMath::number<double>::convert_block(MandelMath::NumberType::typeEmpty, nullptr, newLength*JuliaPoint<MandelMath::number_any>::LEN);
        break;
#if !NUMBER_DOUBLE_ONLY
      case MandelMath::NumberType::typeFloat128:
        new_points=MandelMath::number<__float128>::convert_block(MandelMath::NumberType::typeEmpty, nullptr, newLength*JuliaPoint<MandelMath::number_any>::LEN);
        break;
      case MandelMath::NumberType::typeDDouble:
        new_points=MandelMath::number<MandelMath::dd_real>::convert_block(MandelMath::NumberType::typeEmpty, nullptr, newLength*JuliaPoint<MandelMath::number_any>::LEN);
        break;
      case MandelMath::NumberType::typeQDouble:
        new_points=MandelMath::number<MandelMath::dq_real>::convert_block(MandelMath::NumberType::typeEmpty, nullptr, newLength*JuliaPoint<MandelMath::number_any>::LEN);
        break;
      case MandelMath::NumberType::typeReal642:
        new_points=MandelMath::number<MandelMath::real642>::convert_block(MandelMath::NumberType::typeEmpty, nullptr, newLength*JuliaPoint<MandelMath::number_any>::LEN);
        break;
#endif
    }
    MandelMath::complex<MandelMath::number_any> old_c(precisionRecord->position.center);
    //old_c.assign(precisionRecord->position.center);

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

    image_dirty.setBase(0, imageWidth, 0, imageHeight);
    startNewEpoch();
  }
}

void JuliaModel::pause(bool pause)
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

void JuliaModel::startNewEpoch()
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

  invalidateMainImage();

  emit triggerJuliaThreaded(epoch, precisionRecord->params.period); //::invokeMethod cannot pass parameters, but ::connect can
#endif
}

void JuliaModel::reimToPixel(int *circ_x, int *circ_y, MandelMath::complex<MandelMath::number_any> const &z, MandelMath::number<MandelMath::number_any> *tmp)
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

void JuliaModel::paintOrbit(ShareableImageWrapper image, int x, int y)
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
  if (precisionRecord->params.period<0)
    return;

  //JuliaPointStore *resultStore=&pointStore[y*imageWidth+x];
  //JuliaPoint *data=&pointStore[y*imageWidth+x];
  /*switch (resultStore->state)
  {
    case JuliaPoint::State::stUnknown:
    case JuliaPoint::State::stResolved:
    case JuliaPoint::State::stFail:
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

    if (precisionRecord->params.period>0)
    {
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
    };
  }
  JuliaPointStore *resultStoreHover=&pointStore[y*imageWidth+x];
  if (!_orbit_frozen)
    precisionRecord->orbit.pointDataStore.assign(resultStoreHover);
  JuliaPointStore *resultStoreOrbit=&precisionRecord->orbit.pointDataStore;
  switch (resultStoreHover->rstate) //TODO: for debugging it might be better to read from evaluator.juliaData.store (after it's computed)
  {
    case JuliaPointStore::ResultState::stOutside:
    case JuliaPointStore::ResultState::stOutAngle:
    { //exterior distance estimate
      int exterior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      painter.setPen(QColor(0xff, 0xff, 0));
      switch (_selectedEdePatchAlgo)
      {
        case EedepaAlways: exterior=qRound(resultStoreHover->exterior.always_mixed/precisionRecord->position.step_size); break;
        case EedepaConditional: exterior=qRound(resultStoreHover->exterior.condi_mixed/precisionRecord->position.step_size); break;
        case EedepaBlend: exterior=qRound(resultStoreHover->exterior.blend_mixed/precisionRecord->position.step_size); break;
        case EedepaVerifyH: exterior=qRound(resultStoreHover->exterior.verify_mixed/precisionRecord->position.step_size); break;
        case EedepaVerifyA: exterior=qRound(resultStoreHover->exterior.verify_mixed/precisionRecord->position.step_size); break;
        default: exterior=qRound(resultStoreHover->exterior.none_mixed/precisionRecord->position.step_size); break;
      }
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
      painter.setPen(QColor(0xc0, 0xc0, 0));
      switch (_selectedEdePatchAlgo)
      {
        case EedepaAlways: exterior=qRound(resultStoreHover->exterior.always_empty/precisionRecord->position.step_size); break;
        case EedepaConditional: exterior=qRound(resultStoreHover->exterior.condi_empty/precisionRecord->position.step_size); break;
        case EedepaBlend: exterior=qRound(resultStoreHover->exterior.blend_empty/precisionRecord->position.step_size); break;
        case EedepaVerifyH: exterior=qRound(resultStoreHover->exterior.verify_empty/precisionRecord->position.step_size); break;
        case EedepaVerifyA: exterior=qRound(resultStoreHover->exterior.verify_empty/precisionRecord->position.step_size); break;
        default: exterior=qRound(resultStoreHover->exterior.none_empty/precisionRecord->position.step_size); break;
      }
      painter.drawEllipse(x-exterior, y-exterior, 2*exterior, 2*exterior);
    } break;
    case JuliaPointStore::ResultState::stPeriod2:
    case JuliaPointStore::ResultState::stPeriod3:
    { //interior distance estimate
      int interior;
      painter.setBrush(Qt::BrushStyle::NoBrush);
      if (_selectedPaintStyle==paintStylePhi2)
      {
        painter.setPen(QColor(0, 0xc0, 0xc0));
        interior=qRound(std::hypot(resultStoreHover->interior.phi2_re, resultStoreHover->interior.phi2_im)/precisionRecord->position.step_size);
        painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);

        //painter.setPen(QColor(0, 0xff, 0xff));
        //interior=qRound(resultStore->interior.mixed/precisionRecord->position.step_size);
        //painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);

        painter.setPen(QColor(0x00, 0x00, 0x00));
        int inte_x=qRound(resultStoreHover->interior.phi2_re/precisionRecord->position.step_size);
        int inte_y=qRound(resultStoreHover->interior.phi2_im/precisionRecord->position.step_size);
        painter.drawEllipse(x+inte_x-1, y-inte_y-1, 3, 3);
        painter.setPen(QColor(0xff, 0xff, 0xff));
        painter.drawEllipse(+x+inte_x-2, y-inte_y-2, 5, 5);
      }
      else if (_selectedPaintStyle==paintStylePhi)
      {
        double phi_mag=MandelMath::sqr_double(resultStoreHover->interior.phi_re)+MandelMath::sqr_double(resultStoreHover->interior.phi_im);
        painter.setPen(QColor(0, 0xc0, 0xc0));
        interior=qRound(1.0/std::sqrt(phi_mag)/precisionRecord->position.step_size);
        painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);

        /*painter.setPen(QColor(0, 0xff, 0xff));
        interior=qRound(resultStoreHover->interior.mixed/precisionRecord->position.step_size);
        painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);*/

        painter.setPen(QColor(0x00, 0x00, 0x00)); //black dot at -1/phi
        int inte_x=qRound(-resultStoreHover->interior.phi_re/phi_mag/precisionRecord->position.step_size);
        int inte_y=qRound( resultStoreHover->interior.phi_im/phi_mag/precisionRecord->position.step_size);
        painter.drawEllipse(x+inte_x-1, y-inte_y-1, 3, 3);
        painter.setPen(QColor(0xff, 0xff, 0xff));
        painter.drawEllipse(+x+inte_x-2, y-inte_y-2, 5, 5);
      }
      else
      {
        painter.setPen(QColor(0, 0xc0, 0xc0));
        interior=qRound(resultStoreHover->interior.full/precisionRecord->position.step_size);
        painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);

        painter.setPen(QColor(0, 0xff, 0xff));
        interior=qRound(resultStoreHover->interior.mixed/precisionRecord->position.step_size);
        painter.drawEllipse(x-interior, y-interior, 2*interior, 2*interior);

        painter.setPen(QColor(0x00, 0x00, 0x00));
        int inte_x=qRound(resultStoreHover->interior.mixed_re/precisionRecord->position.step_size);
        int inte_y=qRound(resultStoreHover->interior.mixed_im/precisionRecord->position.step_size);
        painter.drawEllipse(x+inte_x-1, y-inte_y-1, 3, 3);
        painter.setPen(QColor(0xff, 0xff, 0xff));
        painter.drawEllipse(+x+inte_x-2, y-inte_y-2, 5, 5);
      }
    } break;
    default: ;
  }


  //evaluator.params.c .. saved mouse coordinates
  //evaluator.data.root .. found root (temporary)
  //pointData.f .. not used
  //pointData.fz .. result f'

  MandelEvaluator<MandelMath::number_any> &evaluator=precisionRecord->orbit.evaluator;
  //TODO: a bit backwards dance with currentParams.c and currentData.f, needs to be straightened
  if (!_orbit_frozen)
  {
    precisionRecord->position.pixelXtoRE(x-imageWidth/2, &evaluator.currentParams.julia.first_z.re);
    precisionRecord->position.pixelYtoIM(imageHeight/2-y, &evaluator.currentParams.julia.first_z.im);
  };
  evaluator.currentParams.epoch=epoch;
  evaluator.workIfEpoch=evaluator.busyEpoch;//epoch;
  evaluator.currentParams.pixelIndex=0;
  //already 0 precisionRecord->orbit.evaluator.currentParams.nth_fz=0;
  evaluator.juliaData.zero(evaluator.currentParams.julia.first_z);
  evaluator.juliaData.store->wstate=JuliaPointStore::WorkState::stWorking;
  evaluator.currentParams.breakOnNewNearest=false;
  evaluator.currentParams.maxiter=1<<MAX_EFFORT;
  evaluator.currentParams.want_extangle=(_selectedPaintStyle==paintStyleExterAngle) ||
                                        (resultStoreOrbit->rstate==JuliaPointStore::ResultState::stOutAngle);
  evaluator.currentParams.julia.c.assign_across(precisionRecord->params.c);

  //precisionRecord->orbit.evaluator.currentData.root.assign(&precisionRecord->orbit.evaluator.currentParams.c);
  //precisionRecord->orbit.evaluator.loope.help_c.assign_across(&precisionRecord->params.base);
  {
    painter.setBrush(Qt::BrushStyle::NoBrush);
    painter.setPen(QColor(0xff, 0xff, 0xff)); //paint path
    int line_sx, line_sy;
    reimToPixel(&line_sx, &line_sy, evaluator.juliaData.f, &tmp);
    while ((evaluator.juliaData.store->rstate==JuliaPointStore::ResultState::stUnknown) &&
           (evaluator.juliaData.store->iter<(1<<MAX_EFFORT)))
    {
      int line_ex, line_ey;

      if (MandelMath::enum_is_one_of<JuliaPointStore::ResultState,
                          JuliaPointStore::ResultState::stPeriod2,
                          JuliaPointStore::ResultState::stPeriod3,
                          JuliaPointStore::ResultState::stUnknown,
                          JuliaPointStore::ResultState::stMaxIter>(resultStoreOrbit->rstate) &&
          evaluator.juliaData.store->iter<precisionRecord->params.period+resultStoreOrbit->near0iter_1)
      { //paint first period fully
        evaluator.currentParams.maxiter=evaluator.juliaData.store->iter+1;
      }
      else if (precisionRecord->params.period<=0)
        evaluator.currentParams.maxiter=evaluator.juliaData.store->iter*2+1;
      else //stop at multiples of lookper, +near0iter_1
        evaluator.currentParams.maxiter=resultStoreOrbit->near0iter_1+((evaluator.juliaData.store->iter-resultStoreOrbit->near0iter_1)/precisionRecord->params.period+1)*precisionRecord->params.period;
      evaluator.thread.syncJulia();

      reimToPixel(&line_ex, &line_ey, evaluator.juliaData.f, &tmp);
      if (line_ex>=-3 && line_ex<=10003 && line_ey>=-3 && line_ey<=10003)
      {
        if (line_sx>=-3 && line_sx<=10003 && line_sy>=-3 && line_sy<=10003)
          painter.drawLine(line_sx, line_sy, line_ex, line_ey);
        line_sx=line_ex;
        line_sy=line_ey;
      };

      if (evaluator.juliaData.store->iter<=precisionRecord->params.period)
      {
        evaluator.bulb.t1.assign(evaluator.currentParams.julia.first_z);
        evaluator.bulb.t1.mul(evaluator.juliaData.fz_z);
        evaluator.bulb.t2.assign(evaluator.juliaData.f);
        //f_z_target=1 evaluator.bulb.t2.mul(f_z_target, tmp);
        evaluator.bulb.t2.rsub(evaluator.bulb.t1); //z*f_z-f*f_z_target
        evaluator.bulb.t1.assign(evaluator.juliaData.fz_z);
        evaluator.bulb.t1.re.add_double(-1); //evaluator.bulb.t1.sub(f_z_target);
        evaluator.bulb.t1.recip();
        evaluator.bulb.t1.mul(evaluator.bulb.t2); //s1=A=(z*f_z-f*1)/(f_z-1)
        reimToPixel(&line_ex, &line_ey, evaluator.bulb.t1, &tmp);
        if (line_ex>=-3 && line_ex<=10003 && line_ey>=-3 && line_ey<=10003)
        {
          painter.drawEllipse(line_ex-1, line_ey-1, 3, 3);
        };
      };
    }
  }

  //circle one iteration forward
  {
    constexpr int straight=10;
    constexpr int diagonal=7;
    constexpr struct {int x; int y;} offsets[]={
      {.x=0,         .y=0},
      {.x=+straight, .y=0},
      {.x=+diagonal, .y=+diagonal},
      {.x=0,         .y=+straight},
      {.x=-diagonal, .y=+diagonal},
      {.x=-straight, .y=0},
      {.x=-diagonal, .y=-diagonal},
      {.x=0,         .y=-straight},
      {.x=+diagonal, .y=-diagonal},
    };
    QPoint aroundBefore[9];
    QPoint aroundAfter[9];
    QPoint aroundDeriv[9];
    double f_re, f_im, fz_re, fz_im, fzz_re, fzz_im;
    for (int around=0; around<=8; around++)
    {
      aroundBefore[around].setX(x+offsets[around].x);
      aroundBefore[around].setY(y-offsets[around].y);
      precisionRecord->position.pixelXtoRE(aroundBefore[around].x()-imageWidth/2, &evaluator.loope.sumA.re);
      precisionRecord->position.pixelYtoIM(imageHeight/2-aroundBefore[around].y(), &evaluator.loope.sumA.im);

      evaluator.loope.eval_zz(std::max(1, precisionRecord->params.period*1),
                              precisionRecord->params.c, evaluator.loope.sumA, false, true);

      int circ_x, circ_y;
      reimToPixel(&circ_x, &circ_y, evaluator.loope.f, &tmp);
      aroundAfter[around].setX(circ_x);
      aroundAfter[around].setY(circ_y);
      if (around==0)
      {
        f_re=evaluator.loope.f.re.toDouble();
        f_im=evaluator.loope.f.im.toDouble();
        fz_re=evaluator.loope.f_z.re.toDouble();
        fz_im=evaluator.loope.f_z.im.toDouble();
        fzz_re=evaluator.loope.f_zz.re.toDouble();
        fzz_im=evaluator.loope.f_zz.im.toDouble();
        aroundDeriv[0].setX(circ_x);
        aroundDeriv[0].setY(circ_y);
      }
      else
      {
        double dre=offsets[around].x*precisionRecord->position.step_size;
        double dim=offsets[around].y*precisionRecord->position.step_size;
        //(dr+i*di)*(dr+i*di)*(zr+i*zi)=(dr*dr-di*di+2*i*di*dr)*(zr+i*zi)=dr*dr*zr+i*dr*dr*zi-di*di*zr-i*di*di*zi+2*i*di*dr*zr-2*di*dr*zi=
        //  =dr*dr*zr-di*di*zr-2*di*dr*zi+i*dr*dr*zi-i*di*di*zi+2*i*di*dr*zr
        double ore=dre*fz_re-dim*fz_im+dre*dre*fzz_re-dim*dim*fzz_re-2*dre*dim*fzz_im;
        double oim=dre*fz_im+dim*fz_re+dre*dre*fzz_im-dim*dim*fzz_im+2*dre*dim*fzz_re;
        evaluator.loope.sumA.zero(f_re+ore, f_im+oim);
        reimToPixel(&circ_x, &circ_y, evaluator.loope.sumA, &tmp);
        aroundDeriv[around].setX(circ_x);
        aroundDeriv[around].setY(circ_y);
      }
    }
    painter.setPen(QColor(0xff, 0x80, 0xff));
    painter.drawPolyline(aroundBefore, 9);
    painter.drawPolyline(aroundAfter, 9);
    painter.setPen(QColor(0xff, 0xc0, 0xff));
    painter.drawPolyline(aroundDeriv, 9);
  }

  if ((evaluator.juliaData.store->rstate==JuliaPointStore::ResultState::stPeriod2) ||
      (evaluator.juliaData.store->rstate==JuliaPointStore::ResultState::stPeriod3))
  {
    int circ_x, circ_y;
    painter.setPen(QColor(0, 0xff, 0xff)); //paint root
    reimToPixel(&circ_x, &circ_y, evaluator.currentParams.julia.root, &tmp);
    if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
    {
      painter.drawEllipse(circ_x-3, circ_y-3, 2*3, 2*3);
      QLine l2[2]={{circ_x-2, circ_y, circ_x+2, circ_y},
                   {circ_x, circ_y-2, circ_x, circ_y+2}};
      painter.drawLines(l2, 2);
    };
    //no bulbs in Julia fractal...yet?
  };
  //beginning of external ray
  if (_selectedPaintStyle==paintStyle::paintStyleExterAngle &&
      evaluator.juliaData.store->rstate==JuliaPointStore::ResultState::stOutAngle)
  {
    int circ_x, circ_y;
    painter.setBrush(Qt::BrushStyle::SolidPattern);

    painter.setBrush(QBrush(QColor(0xff, 0xff, 0x80)));
    if (evaluator.extangle.dbg.first_guess_valid>0)
    {
      painter.setPen(QColor(0xff, 0xff, 0xff)); //first guess of first step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.first_guess_0, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (evaluator.extangle.dbg.last_guess_valid>0)
    {
      painter.setPen(QColor(0xff, 0xff, 0x80)); //last guess of first step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.last_guess_0, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    painter.setBrush(QBrush(QColor(0xff, 0xc0, 0x80)));
    if (evaluator.extangle.dbg.first_guess_valid>1)
    {
      painter.setPen(QColor(0xff, 0xc0, 0xff)); //first guess of second step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.first_guess_1_, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (evaluator.extangle.dbg.last_guess_valid>1)
    {
      painter.setPen(QColor(0xff, 0xc0, 0x80)); //last guess of second step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.last_guess_1_, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    painter.setBrush(QBrush(QColor(0xff, 0x80, 0x80)));
    if (evaluator.extangle.dbg.first_guess_valid>2)
    {
      painter.setPen(QColor(0xff, 0x80, 0xff)); //first guess of third step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.first_guess_2, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };

    if (evaluator.extangle.dbg.last_guess_valid>2)
    {
      painter.setPen(QColor(0xff, 0x80, 0x80)); //last guess of third step of external ray
      reimToPixel(&circ_x, &circ_y, evaluator.extangle.dbg.last_guess_2, &tmp);
      if ((circ_x>=-3) && (circ_x<=10003) && (circ_y>=-3) && (circ_y<=10003))
        painter.drawEllipse(circ_x-2, circ_y-2, 2*2, 2*2);
    };
  };

  image.image->swap(newOverlay);
}

int JuliaModel::periodToIndex(int period)
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

int JuliaModel::writeToImage(ShareableImageWrapper image)
{
  //QImage result(imageWidth, imageHeight, QImage::Format::Format_ARGB32);//setPixel in _Premultiplied affects surrounding pixels?! _Premultiplied);
  if (image.image->isNull() || (image.image->width()!=imageWidth) || (image.image->height()!=imageHeight))
    return -1;
  timerWriteToImage.start();
  //double params_nth_fz_limit=precisionRecord->params.nth_fz_limit.toDouble();
  double extAngleZoom=1<<this->_extAngleZoom;
  JuliaPointStore *wtiStore;
  //static constexpr unsigned int FAIL_COLOR=0xff404040;

  //need some std::atomic expert here...
  image_dirty.fetch();

  //precisionRecord->wtiPoint.self_allocator._getFirstCapac(indexOfWtiPoint, _discard_);
  if (image_dirty.fleft<=image_dirty.fright && image_dirty.ftop<=image_dirty.fbottom)
  for (int y=image_dirty.ftop; y<=image_dirty.fbottom; y++)
    for (int x=image_dirty.fleft; x<=image_dirty.fright; x++)
    {
      //MandelMath::worker_multi::Allocator allo(storeWorker->getAllocator(), (y*imageWidth+x)*JuliaPoint::LEN, JuliaPoint::LEN, nullptr);
      //JuliaPoint data_(&pointStore_[y*imageWidth+x], &allo);
      wtiStore=&pointStore[y*imageWidth+x];
      //if ((x==222) && (y==170))
        //data->state=MandelPoint::State::stOutside;
      switch (_selectedPaintStyle)
      {
        case paintStyle::paintStyleKind:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              if (wtiStore->wstate==JuliaPointStore::WorkState::stIdle)
                image.image->setPixel(x, y, 0xffffffff);
              else
                image.image->setPixel(x, y, 0xffff00ff);
                //image.image->setPixel(x, y, 0xffffffff);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              int b=0x9f+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b<<0));
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff8000ff);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffc000);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              /*int r;
              switch (data->period)
              {
                case 1: r=0xc0; break;
                case 2: r=0xff; break;
                case 3: r=0x80; break;
                default: r=0xe0;
              }*/
              /*if (wtiStore->period!=precisionRecord->params.period)//if (wtiStore->period>wtiStore->nearziter_0)
                image.image->setPixel(x, y, 0xffff00ff); //seems to only happen by mistake, not in reality
              else
              {
                int index=periodToIndex(wtiStore->period);
                //reverse bottom 7 bits:
                int r=0x80 | MandelMath::ReverseBits<7,1>(index);
                image.image->setPixel(x, y, 0xff000000+(r<<16));
              }*/
              double iterx=std::log(wtiStore->iter+1);
              int r=(iterx-floor(iterx))*256;
              image.image->setPixel(x, y, 0xff000000+(r<<16));
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleCls:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00906090);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
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
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
              double re=precisionRecord->wtiPoint.f.re.toDouble();
              double im=precisionRecord->wtiPoint.f.im.toDouble();
              double iter=wtiStore->iter+6-log2(log2(re*re+im*im)); //+6 to match integer coloring
              iter=sqrt(iter);
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
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xffffffff);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              //ok but shouldn't use raw iter int index=periodToIndex(wtiStore->iter % wtiStore->period);
              //ok int index=periodToIndex(wtiStore->nearziter_0 % wtiStore->period); fails with c around -1+0i : all chunks at phase=0
              //int index=periodToIndex(wtiStore->near0iter_1 % wtiStore->period);//also ok?
              int index;
              if (precisionRecord->params.period<=0)
                index=wtiStore->near0iter_1;
              else
                index=wtiStore->near0iter_1 % precisionRecord->params.period;//also ok?
              /*if (wtiStore->surehand==0 || wtiStore->surehand%wtiStore->period!=0)
              {
                image.image->setPixel(x, y, 0xffff0080);
              }
              else*/
              {
                int r=0x80 | MandelMath::ReverseBits<7,1>(index);
                /*switch (data->period)
                {
                  case 1: r=0xc0; break;
                  case 2: r=0xff; break;
                  case 3: r=0x80; break;
                  default: r=0xe0;
                }*/
                //image.image->setPixel(x, y, 0xffffc0c0);
                image.image->setPixel(x, y, 0xff000000+r*0x010101);
              }
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
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
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              double tf, exterior;
              switch (_selectedEdePatchAlgo)
              {
                case EedepaAlways: exterior=wtiStore->exterior.always_mixed; break;
                case EedepaConditional: exterior=wtiStore->exterior.condi_mixed; break;
                case EedepaBlend: exterior=wtiStore->exterior.blend_mixed; break;
                case EedepaVerifyH: exterior=std::abs(wtiStore->exterior.verify_mixed/wtiStore->exterior.none_mixed-1); break;
                case EedepaVerifyA: exterior=std::abs(wtiStore->exterior.verify_empty/wtiStore->exterior.none_empty-1); break;
                default: exterior=wtiStore->exterior.none_mixed; break;
              }
              if ((exterior>10000) || (exterior<=0))
                image.image->setPixel(x, y, 0xff9f9f9f);
              else if (_selectedExteriorColoring==EextcolScaledWaves)
              {
                if (exterior>=1) //10000..1 -> -9999..0
                {
                  tf=(1-exterior)*1;
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
                  tf=-log(exterior);//sqrt(1-log(wtiStore->exterior_avoids))*2-2;
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
              }
              else //(_selectedExteriorColoring==EextcolRainbow)
              {
                int ti;
                if (exterior<=0)
                  ti=0;
                else if (exterior>=4)
                  ti=6*0xc0-1-(qRound(log(exterior/4)*600)) % (6*0xc0);
                else
                  ti=(qRound(-log(exterior/4)*600)+12*0xc0) % (6*0xc0);
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
              }
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              int index;
              if (precisionRecord->params.period<=0)
                index=wtiStore->near0iter_1;
              else
                index=wtiStore->near0iter_1 % precisionRecord->params.period;//also ok?
              int r=0x80 | MandelMath::ReverseBits<7,1>(index);
              image.image->setPixel(x, y, 0xff000000+r*0x010101);
              /*
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
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
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+b);*/
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleExterAngle:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
              image.image->setPixel(x, y, 0x00606060);
              break;
            case JuliaPointStore::ResultState::stOutAngle:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
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
                int r=(tf<-3.1415926536 || tf>3.1415926535)?0xff:(tf<0)?0x80:0;//qRound(log(wtiStore->iter)*100)%256;
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
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              int index;
              if (precisionRecord->params.period<=0)
                index=wtiStore->nearziter_0;
              else
                index=wtiStore->nearziter_0 % precisionRecord->params.period;//also ok?
              int r=0x80 | MandelMath::ReverseBits<7,1>(index);
              image.image->setPixel(x, y, 0xff000000+r*0x010101);
              /*
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
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
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+b);*/
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleInter:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              //int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              int b=qRound(log(wtiStore->iter)*100)%256;
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              int ti=30;
              if (wtiStore->interior.mixed<=0)
                ti=0;
              else if (wtiStore->interior.mixed>4)
                //ti=(qRound((1-wtiStore->interior.hits/4)*300)+12*0xc0) % (6*0xc0);
                ti=6*0xc0-1-(qRound((wtiStore->interior.mixed/4-1)*300)+1)%(6*0xc0);
              else
                ti=(qRound(-log(wtiStore->interior.mixed/4)*300)+12*0xc0) % (6*0xc0);
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
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleNear0:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0xff000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              //int ti=wtiStore->nearziter_0;
              int ti;
              if (precisionRecord->params.period<=0)
                ti=wtiStore->near0iter_1;
              else
                ti=wtiStore->near0iter_1;// % precisionRecord->params.period;
              image.image->setPixel(x, y, 0xff000080+(MandelMath::ReverseBits<7,1>(ti)<<0));
              /*if (ti>=0)
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
                image.image->setPixel(x, y, 0xff000080);*/
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
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

              //int ti=wtiStore->nearziter_0;
              int ti;
              if (precisionRecord->params.period<=0)
                ti=wtiStore->near0iter_1;
              else
                ti=wtiStore->near0iter_1;// % precisionRecord->params.period;
              image.image->setPixel(x, y, qRgb(0x80+MandelMath::ReverseBits<7,1>(ti), 0, 0));
              /*if (ti>0)
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
                image.image->setPixel(x, y, 0xff800000);*/
              /* we need func(2)!=func(3) here
              int index=periodToIndex(data->near0iter);
              int r=0x80 | MandelMath::ReverseBits<7,1>(index);
              image.image->setPixel(x, y, 0xff000000+(r<<16));*/
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleNearZ:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0xff000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              //int ti=wtiStore->nearziter_0;
              int ti;
              if (precisionRecord->params.period<=0)
                ti=wtiStore->nearziter_0;
              else
                ti=wtiStore->nearziter_0 % precisionRecord->params.period;
              image.image->setPixel(x, y, 0xff000080+(MandelMath::ReverseBits<7,1>(ti)<<0));
              /*if (ti>=0)
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
                image.image->setPixel(x, y, 0xff000080);*/
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff00ff00);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff00c000);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff008000);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
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

              //int ti=wtiStore->nearziter_0;
              int ti;
              if (precisionRecord->params.period<=0)
                ti=wtiStore->nearziter_0;
              else
                ti=wtiStore->nearziter_0 % precisionRecord->params.period;
              image.image->setPixel(x, y, 0xff800000+(MandelMath::ReverseBits<7,1>(ti)<<16));
              /*if (ti>0)
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
                image.image->setPixel(x, y, 0xff800000);*/
              /* we need func(2)!=func(3) here
              int index=periodToIndex(data->near0iter);
              int r=0x80 | MandelMath::ReverseBits<7,1>(index);
              image.image->setPixel(x, y, 0xff000000+(r<<16));*/
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStylePhi:
        case paintStyle::paintStylePhi1:
        case paintStyle::paintStylePhi2:
        case paintStyle::paintStyleF1:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
              //double re=precisionRecord->wtiPoint.fz_r.re.toDouble();
              //double im=precisionRecord->wtiPoint.fz_r.im.toDouble(); //fz_r
              double re, im;
              switch (_selectedPaintStyle)
              {
                case paintStylePhi: re=wtiStore->interior.phi_re; im=wtiStore->interior.phi_im; break;
                case paintStylePhi1: re=wtiStore->interior.phi1_re; im=wtiStore->interior.phi1_im; break;
                case paintStylePhi2: re=wtiStore->interior.phi2_re; im=wtiStore->interior.phi2_im; break;
                case paintStyleF1: re=wtiStore->interior.f1_re; im=wtiStore->interior.f1_im; break;
                default: re=0; im=0;
              }
              double mag=sqrt(MandelMath::sqr_double(re)+MandelMath::sqr_double(im));
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
              //mag=log(1+mag)/log(2);
              mag=log(mag)/log(2);
              mag=mag-floor(mag);
              if (mag<0.2)
                r=128;
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleFirstUnder1:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              int b=0x60+floor(0x60*cos((wtiStore->iter/10.0+0)*2*3.1415926535));
              image.image->setPixel(x, y, 0xff000000+(b*0x010101));
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
              //int index=wtiStore->interior.first_under_1;
              int index=wtiStore->bigfzfzziter;
              //reverse bottom 7 bits:
              int g=0x80 | MandelMath::ReverseBits<7,1>(index);
              int b=0;
              int r=0;
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;
        case paintStyle::paintStyleNearMR:
        {
          switch (wtiStore->rstate)
          {
            case JuliaPointStore::ResultState::stUnknown:
              image.image->setPixel(x, y, 0x00000000);
              break;
            case JuliaPointStore::ResultState::stOutside:
            case JuliaPointStore::ResultState::stOutAngle:
            {
              int index=wtiStore->nearmriter;
              //reverse bottom 7 bits:
              int g=0;
              if (precisionRecord->params.period>0)
              {
                if ((wtiStore->near0iter_1-wtiStore->nearmriter) % precisionRecord->params.period!=0)
                  g|=0x60;
              };
              int b=0x80 | MandelMath::ReverseBits<7,1>(index);
              int r=0;
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case JuliaPointStore::ResultState::stBoundary:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stMisiur:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stDiverge:
            {
              image.image->setPixel(x, y, 0xff308030);
            } break;
            case JuliaPointStore::ResultState::stPeriod2:
            case JuliaPointStore::ResultState::stPeriod3:
            {
              //precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*JuliaPoint<MandelMath::number_any>::LEN);
              //int index=wtiStore->interior.first_under_1;
              //would have to copy wtiPoint into precisionRecord->wtiPoint->store but that's nullptr     int index=precisionRecord->wtiPoint.nearmr.get();

              int index=wtiStore->nearmriter;
              //reverse bottom 7 bits:
              int g=0;
              if (precisionRecord->params.period>0)
              {
                if ((wtiStore->near0iter_1-wtiStore->nearmriter) % precisionRecord->params.period!=0)
                  g|=0x60;
              };
              int b=0;
              if (wtiStore->nearmr_difficult)
                b|=0xa0;
              int r=0x80 | MandelMath::ReverseBits<7,1>(index);
              image.image->setPixel(x, y, 0xff000000+(r<<16)+(g<<8)+(b));
            } break;
            case JuliaPointStore::ResultState::stMaxIter:
            {
              image.image->setPixel(x, y, 0xff808080);
            } break;
          }
        } break;

        /*case paintStyle::paintStyleFC:
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
                precisionRecord->wtiPoint.readFrom(precisionRecord->points, (y*imageWidth+x)*MandelPoint<MandelMath::number_any>::LEN);
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
        } break;*/
      }
    }
  //return result;
  qint64 elapsed=timerWriteToImage.elapsed();
  return elapsed;
}

void JuliaModel::doneWorkInThread(MandelEvaluatorThread *)
{
  _threadsWorking--;
}

template <typename BASE>
int JuliaModel::giveWorkThreaded(MandelEvaluator<BASE> *me)
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
      JuliaPointStore *storeAtIndex=&pointStore[pointIndex];
      if (storeAtIndex->wstate.load(std::memory_order_relaxed)==JuliaPointStore::WorkState::stWorking)
        continue; //speedup/precheck
      int extra_effort=0;
      JuliaPointStore::WorkState state_prev=storeAtIndex->wstate.exchange(JuliaPointStore::WorkState::stWorking, std::memory_order_acquire);
      if (state_prev==JuliaPointStore::WorkState::stWorking)
        continue;
      if (state_prev==JuliaPointStore::WorkState::stDone)
      {
        /*if ((_selectedPaintStyle==paintStyleFC) &&
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
        else*/ if (_selectedPaintStyle==paintStyleExterAngle &&
                 storeAtIndex->rstate==JuliaPointStore::ResultState::stOutside)
        {
          extra_effort=1;
        }
        else
        {
          storeAtIndex->wstate.store(JuliaPointStore::WorkState::stDone, std::memory_order_relaxed);
          continue;
        }
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
              storeAtIndex->rstate=JuliaPointStore::ResultState::stMaxIter;
              storeAtIndex->wstate.store(JuliaPointStore::WorkState::stDone, std::memory_order_release);
            }
            else
            {
              if (retryEffortFrom<0)
                retryEffortFrom=pointIndex;
              storeAtIndex->wstate.store(JuliaPointStore::WorkState::stIdle, std::memory_order_release);
            }
          }
          else
          {
            //if (me->currentWorker->ntype()!=precisionRecord->ntype)
              //dbgPoint();
            //evaluator->switchType(position.worker);
            precisionRecord->position.pixelXtoRE<BASE>(pointIndex%imageWidth - imageWidth/2, &me->currentParams.julia.first_z.re);
            precisionRecord->position.pixelYtoIM<BASE>(imageHeight/2-pointIndex/imageWidth, &me->currentParams.julia.first_z.im);
            //pre-assigned in setParams me->currentParams.c_.assign_across(&precisionRecord->params.base);
            me->currentParams.epoch=me->busyEpoch;
            //storeAtIndex->state=MandelPointStore::State::stWorking;
            me->currentParams.pixelIndex=pointIndex;
            me->currentParams.want_fc_r=false;//(_selectedPaintStyle==paintStyleFC);
            me->currentParams.want_extangle=(_selectedPaintStyle==paintStyleExterAngle);
#if CURRENT_STORE_DIRECT
            me->currentData.store=storeAtIndex;
#else
            me->juliaStore.assign(storeAtIndex);
#endif
            me->juliaData.readFrom(precisionRecord->points, pointIndex*JuliaPoint<MandelMath::number_any>::LEN);
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
int JuliaModel::doneWorkThreaded(MandelEvaluator<BASE> *me, bool giveWork)
{
  QReadLocker locker(&threading_mutex);
  if ((me->currentParams.epoch==epoch) && (me->currentParams.pixelIndex>=0) && (me->currentParams.pixelIndex<imageWidth*imageHeight))
  {
    JuliaPointStore *dstStore=&pointStore[me->currentParams.pixelIndex];
    if (dstStore->wstate.load(std::memory_order_relaxed)!=JuliaPointStore::WorkState::stWorking)
      dbgPoint(); //leftovers should be from different epoch
    me->juliaData.writeTo(precisionRecord->points, me->currentParams.pixelIndex*JuliaPoint<MandelMath::number_any>::LEN);
    dstStore->assign(me->juliaData.store);
    if (dstStore->wstate.load(std::memory_order_relaxed)==JuliaPointStore::WorkState::stIdle)
      dbgPoint();
    if (dstStore->wstate.load(std::memory_order_relaxed)==JuliaPointStore::WorkState::stWorking)
    {
       if (dstStore->iter>=(1<<MAX_EFFORT))
       {
         dstStore->rstate=JuliaPointStore::ResultState::stMaxIter;
         dstStore->wstate.store(JuliaPointStore::WorkState::stDone, std::memory_order_release);
       }
       else if (dstStore->rstate!=JuliaPointStore::ResultState::stUnknown)
         dstStore->wstate.store(JuliaPointStore::WorkState::stDone, std::memory_order_release);
       else
         dstStore->wstate.store(JuliaPointStore::WorkState::stIdle, std::memory_order_release);

       int x=me->currentParams.pixelIndex%imageWidth;
       int y=me->currentParams.pixelIndex/imageWidth;
       image_dirty.addPoint(x, y);
    }
    else
      dbgPoint();
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

void JuliaModel::selectedPrecisionChanged()
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
    case EprecisionDouble: lolwut:
      newPrecision=Type::typeDouble;
      new_points=MandelMath::number<double>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*JuliaPoint<double>::LEN);
      break;
#if !NUMBER_DOUBLE_ONLY
    case EprecisionFloat128:
      newPrecision=Type::typeFloat128;
      new_points=MandelMath::number<__float128>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*JuliaPoint<__float128>::LEN);
      break;
    case EprecisionDDouble:
      newPrecision=Type::typeDDouble;
      new_points=MandelMath::number<MandelMath::dd_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*JuliaPoint<MandelMath::dd_real>::LEN);
      break;
    case EprecisionQDouble:
      newPrecision=Type::typeQDouble;
      new_points=MandelMath::number<MandelMath::dq_real>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*JuliaPoint<MandelMath::dq_real>::LEN);
      break;
    case EprecisionReal642:
      newPrecision=Type::typeReal642;
      new_points=MandelMath::number<MandelMath::real642>::convert_block(precisionRecord?precisionRecord->ntype:Type::typeEmpty, precisionRecord?precisionRecord->points:nullptr, pointCount*JuliaPoint<MandelMath::real642>::LEN);
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

void JuliaModel::selectedEdePatchAlgoChanged()
{

}

void JuliaModel::selectedExteriorColoringChanged()
{

}

JuliaModel::Params::Params(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Params *source):
  period(source?source->period:1), c(spad), root(spad)
{
  if (source)
  {
    c.assign_across(source->c);
    root.assign_across(source->root);
  }
  else
  {
    c.zero(0, 0);
    root.zero(0, 0);
  }
}


JuliaModel::Position::Position(MandelMath::number<MandelMath::number_any>::Scratchpad *spad, const Position *source):
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

JuliaModel::Position::~Position()
{
}

void JuliaModel::Position::setView(MandelMath::complex<MandelMath::number_any> const &c, double scale)
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

void JuliaModel::Position::move(int delta_x, int delta_y)
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

void JuliaModel::Position::scale(int inlog, int center_x, int center_y)
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

void JuliaModel::Position::updateCachedDepth()
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
void JuliaModel::Position::pixelXtoRE(int x, MandelMath::number<BASE> *result)
{
  result->assign_across(center.re);
  result->add_double(x*step_size);
}

template <typename BASE>
void JuliaModel::Position::pixelYtoIM(int y, MandelMath::number<BASE> *result)
{
  result->assign_across(center.im);
  result->add_double(y*step_size);
}



JuliaModel::Orbit::Orbit(MandelMath::NumberType ntype, Orbit const *source):
  evaluator(ntype, true),
  pointDataStore(),
  first_mu_re(0), first_mu_im(0), first_mum_re(0), first_mum_im(0)
{
  if (source)
  {
    evaluator.currentParams.julia.c.assign_across(source->evaluator.currentParams.julia.c);
  };
}

JuliaModel::Orbit::~Orbit()
{
  evaluator.workIfEpoch=-1;
  evaluator.thread.quit();
  evaluator.thread.wait(1000);
}


JuliaModel::PrecisionRecord::PrecisionRecord(MandelMath::NumberType ntype, PrecisionRecord *source, JuliaModel *doneReceiver):
  ntype(ntype), orbit(ntype, source?&source->orbit:nullptr),
  wtiPoint(nullptr, &orbit.evaluator.tmp),
  params(&orbit.evaluator.tmp, source?&source->params:nullptr),
  position(&orbit.evaluator.tmp, source?&source->position:nullptr),
  //tmp_place(&orbit.evaluator.tmp),
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
  MandelMath::visit([](auto &threads, PrecisionRecord &self, JuliaModel *doneReceiver, PrecisionRecord *source){
      for (int t=0; t<self.threadCount; t++)
      {
        if constexpr (!std::is_same_v<std::decay_t<decltype(threads)>, std::nullptr_t>)
        {
          using Evaluator=std::remove_pointer_t<std::remove_pointer_t<std::remove_reference_t<decltype(threads)>>>;
          Evaluator *thread=new Evaluator(self.ntype, source==nullptr);
          threads[t]=thread;
          if (source)
            thread->currentParams.julia.assign_across(source->orbit.evaluator.currentParams.julia);
          thread->threaded.give=[doneReceiver](Evaluator *me)
          {
            return doneReceiver->giveWorkThreaded(me);
          };
          thread->threaded.doneJulia=[doneReceiver](Evaluator *me, bool giveWork)
          {
            return doneReceiver->doneWorkThreaded(me, giveWork);
          };
          QObject::connect(&thread->thread, &MandelEvaluatorThread::doneJuliaThreaded,
                           doneReceiver, &JuliaModel::doneWorkInThread,
                           Qt::ConnectionType::QueuedConnection);
          QObject::connect(doneReceiver, &JuliaModel::triggerJuliaThreaded,
                           &thread->thread, &MandelEvaluatorThread::doJuliaThreaded,
                           Qt::ConnectionType::QueuedConnection);
        }
      }
    }, threads, *this, doneReceiver, source);
}

JuliaModel::PrecisionRecord::~PrecisionRecord()
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

  switch (ntype)
  {
  case MandelMath::NumberType::typeEmpty: goto ughwut;
  case MandelMath::NumberType::typeDouble: ughwut:
      delete[] (double *)points;
      break;
#if !NUMBER_DOUBLE_ONLY
  case MandelMath::NumberType::typeFloat128:
      delete[] (__float128 *)points;
      break;
  case MandelMath::NumberType::typeDDouble:
      delete[] (MandelMath::dd_real *)points;
      break;
  case MandelMath::NumberType::typeQDouble:
      delete[] (MandelMath::dq_real *)points;
      break;
  case MandelMath::NumberType::typeReal642:
      delete[] (MandelMath::real642 *)points;
      break;
#endif
  }
}
