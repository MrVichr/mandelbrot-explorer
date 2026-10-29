#include "double_double.hpp"
#include <cctype>
#include <clocale>
#include <cmath>
#include <string>
#include <quadmath.h>

void __nop()
{

}

void dd_dbgPoint()
{
  __nop();
}

#if defined(_WIN32) && !defined(__GNUC__)
#ifdef __BORLANDC__
#else
  /* Win 32 MSVC */
  #include <float.h>
#endif
#endif

namespace MandelMath {

// https://www.davidhbailey.com/dhbsoftware/
// https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz
// https://www.davidhbailey.com/dhbsoftware/dqfun-v01.tar.gz

#define X86

void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#if defined(_WIN32) && !defined(__GNUC__)
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  #ifndef _WIN64 //why can't switch precision? why do they return value that is not in use?
  if ((cw & _MCW_PC) != _PC_53)
    _control87(_PC_53, _MCW_PC);
  #endif
  if ((cw&_MCW_RC)!=_RC_NEAR)
    _control87(_RC_NEAR, _MCW_RC);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  //_FPU_GETCW(cw);
  asm volatile ("fnstcw %0":"=m" (cw));


  new_cw = (cw & ~0x0300) | 0x0200;
  //_FPU_SETCW(new_cw);
  asm volatile ("fldcw %0": :"m" (new_cw));

  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#if defined(_WIN32) && !defined(__GNUC__)

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    //_FPU_SETCW(cw);
    asm volatile ("fldcw %0": :"m" (cw));
  }
#endif
#endif
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

dd_real &dd_real::round()
{
  //if lo==+-0.5 exactly, we should decide based on hi whether to round to +-1 or 0 but who cares
  //or could do this+=2^103-=2^103 or something like that
  if ((lo_<=-0.5) || (lo_>=0.5))
  {
    //hi=round(hi);
    lo_=std::round(lo_);
    int discard;
    double testme=std::frexp(lo_, &discard);
    if (testme==-1 || testme==-0.5 || testme==0.5 || testme==1) //spec says it returns +-0.5 instead of 1 but that's just stupid
      dd_dbgPoint(); //may need renormalization
  }
  else
  {
    hi=std::round(hi);
    lo_=0;
    //shoud be good even at hi= +-2^52
  }
  return *this;
}

dd_real &dd_real::frac()
{
  if ((lo_<=-0.5) || (lo_>=0.5))
  {
    if ((hi<0))// && (lo_>=0.5))
    { //-1000+0.75 -> -0.25
      //-1000-0.75 -> -0.75
      hi=lo_-std::ceil(lo_);
      lo_=0;
    }
    else //if ((hi>0) && (lo_<=-0.5))
    { //1000-0.75 -> 0.25
      //1000+0.75 -> 0.75
      hi=lo_-std::floor(lo_);
      lo_=0;
    }
    /*else if (lo_>0)
    { //1000+0.75
      hi=lo_-floor(lo_);
      lo_=0;
    }
    else
    { //-1000-0.75 -> -0.75
      hi=lo_-ceil(lo_);
      lo_=0;
    }*/
  }
  else
  {
    if (hi<0)
      add_double(-std::ceil(hi));
    else
      add_double(-std::floor(hi));
  }
  return *this;
}

dd_real &dd_real::mod1()
{
  if ((lo_<=-0.5) || (lo_>=0.5))
  {
    //1000+0.75->0.75
    //1000-0.75->0.25
    //-1000+0.75->0.75
    //-1000-0.75->0.25
    hi=lo_-floor(lo_);
    lo_=0;
  }
  else
  {
    add_double(-std::floor(hi));
  }
  return *this;
}

bool dd_real::reduce_angle()
{
  static const dd_real M_PIdd(7074237752028440.0/2251799813685248.0, 2483878800010755.0/20282409603651670423947251286016.0);
  static const dd_real M_PIdd_chs(-M_PIdd.hi, -M_PIdd.lo_);
  static const dd_real M_PIdd_2(2*M_PIdd.hi, 2*M_PIdd.lo_);
  //static const dd_real M_PIdd_2chs(-2*3.141, -2*0.0005926);
  if (M_PIdd_chs.isle(*this))
  {
    add(M_PIdd_2);
    return true;
  }
  else if (isle(M_PIdd))
  {
    sub(M_PIdd_2);
    return true;
  }
  else
    return false;
}

dd_real &dd_real::add_pi(double x)
{
  dd_real M_PIdd(7074237752028440.0/2251799813685248.0, 2483878800010755.0/20282409603651670423947251286016.0);
  M_PIdd.mul_double(x);
  add(M_PIdd);
  return *this;
}

// libquadmath follows LC_NUMERIC, which Qt sets from the environment (e.g. ',')
static char locale_decimal_point()
{
  char const *dp=std::localeconv()->decimal_point;
  return (dp && dp[0] && !dp[1])?dp[0]:'.';
}

bool float128_fromString(char const *str, __float128 &result)
{
  std::string s(str);
  char dp=locale_decimal_point();
  if (dp!='.')
    for (char &c: s)
      if (c=='.')
        c=dp;
  char *end;
  result=strtoflt128(s.c_str(), &end);
  if (end==s.c_str())
    return false;
  while (std::isspace((unsigned char)*end))
    end++;
  return *end==0;
}

QString float128_toString(__float128 x, int prec)
{
  char buf[80];
  quadmath_snprintf(buf, sizeof(buf), "%.*Qg", prec, x);
  char dp=locale_decimal_point();
  if (dp!='.')
    for (char *c=buf; *c; c++)
      if (*c==dp)
        *c='.';
  return QString(buf);
}

QString float128_toString(__float128 x)
{
  QString result;
  for (int prec=33; prec<=36; prec++) // 36 digits always round-trip a binary128
  {
    result=float128_toString(x, prec);
    __float128 back;
    if (float128_fromString(result.toUtf8().constData(), back) && back==x)
      break;
  }
  return result;
}

bool dd_real::fromString(char const *str)
{
  __float128 q;
  if (!float128_fromString(str, q))
    return false;
  hi=(double)q;
  lo_=(double)(q-hi);
  return true;
}

QString dd_real::toString() const
{
  __float128 q=(__float128)hi+lo_; // exact for a normalized dd (<=107 bits)
  QString result;
  for (int prec=31; prec<=36; prec++)
  {
    result=float128_toString(q, prec);
    dd_real back;
    if (back.fromString(result.toUtf8().constData()) && back.hi==hi && back.lo_==lo_)
      break;
  }
  return result;
}


/*
inline doubledouble recip(const doubledouble& y) {
  x86_FIX
  double  hc, tc, hy, ty, C, c, U, u;
  C = 1.0/y.h();
  c = doubledouble::Split*C;
  hc =c-C;
  u = doubledouble::Split*y.h();
  hc = c-hc; tc = C-hc; hy = u-y.h(); U = C*y.h(); hy = u-hy; ty = y.h()-hy;
  u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
  c = ((((1.0-U)-u))-C*y.l())/y.h();
  doubledouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
  END_x86_FIX
  return z;
}

inline doubledouble operator /(const doubledouble& x,const doubledouble& y ) {
  x86_FIX
  double hc, tc, hy, ty, C, c, U, u;
  C = x.hi/y.hi; c = doubledouble::Split*C; hc =c-C;  u = doubledouble::Split*y.hi; hc = c-hc;
  tc = C-hc; hy = u-y.hi; U = C * y.hi; hy = u-hy; ty = y.hi-hy;
  u = (((hc*hy-U)+hc*ty)+tc*hy)+tc*ty;
  c = ((((x.hi-U)-u)+x.lo)-C*y.lo)/y.hi;
  doubledouble z; z.hi = C+c; z.lo = double(C-z.hi)+c;
  END_x86_FIX
  return z;
}
*/

#undef MAKE_DOUBLE

} // namespace MandelMath
