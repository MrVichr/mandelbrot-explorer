#include "double_double.hpp"
#include <cmath>

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

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline void dd_real::quick_two_sum(double a, double b)
{
  hi = a + b;
  lo_ = b + (a - hi);
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline void dd_real::two_sum(double a, double b)
{
  hi = a + b;
  double bb = hi - a;
  lo_ = (a - (hi - bb)) + (b - bb);
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

/* Computes high word and lo word of a */
/*inline*/ void dd_real::split(double a)
{
  double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo_ = a - hi;
    hi *= 268435456.0;          // 2^28
    lo_ *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a; //a=a1*2^27+a2*1  temp=a1*2^54+(a1+a2)*2^27 +(a2)*1
    hi = temp - (temp - a);  //temp-a=a1*2^54+a2*2^27    temp-(temp-a)=a1*2^27
    lo_ = a - hi;
  }
}

/* Computes fl(a*b) and err(a*b). */
inline void dd_real::two_prod(double a, double b)
{
  dd_real aa;
  dd_real bb;
  hi = a * b;
  aa.split(a);
  bb.split(b);
  lo_ = ((aa.hi * bb.hi - hi) + aa.hi * bb.lo_ + aa.lo_ * bb.hi) + aa.lo_ * bb.lo_;
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline void dd_real::two_sqr(double a)
{
  dd_real aa;
  hi = a * a;
  aa.split(a);
  lo_ = ((aa.hi * aa.hi - hi) + 2.0 * aa.hi * aa.lo_) + aa.lo_ * aa.lo_;
}

//hack if you don't want to mess with FPU control word... but is a nightmare
#define MAKE_DOUBLE( tmp, val ) ( tmp=(val), *(volatile double *)&tmp )

double dd_real::eps2() const
{
  return 6.07717e-64; // 2^-(2*(53+52)) rounded up
}

double dd_real::eps234() const
{
  return 3.87e-48; // eps2^(3/4)
}

dd_real &dd_real::lshift(int exp)
{
  //lolwut "On many implementations, std::ldexp is less efficient than multiplication or division
  //        by a power of two using arithmetic operators"
  //should I just adjust exponent in *(int *)&hi or what
  hi=ldexp(hi, exp);
  lo_=ldexp(lo_, exp);
  return *this;
}

dd_real &dd_real::add_double(double h2)
{
  dd_real se;
  se.two_sum(hi, h2);
  se.lo_ += lo_;
  quick_two_sum(se.hi, se.lo_);
  return *this;
}

dd_real &dd_real::mul_double(double h2)
{
  dd_real p;
  p.two_prod(hi, h2);
  p.lo_ += lo_*h2;
  quick_two_sum(p.hi, p.lo_);
  return *this;
}

dd_real &dd_real::add(dd_real const &other)
{
//dqfun.f90#dqadd is the same as dd_inline.h#dd_real::sloppy_add
  dd_real se;

  se.two_sum(hi, other.hi);
  se.lo_ += (lo_ + other.lo_);
  quick_two_sum(se.hi, se.lo_);

//1998 inline.h#operator +(const doubledouble& x is the same as dd_inline.h#dd_real::ieee_add
/*
  dd_real s;
  dd_real t;

  s.two_sum(hi, h2);
  t.two_sum(lo, l2);
  s.lo += t.hi;
  s.quick_two_sum(s.hi, s.lo);
  s.lo += t.lo;
  quick_two_sum(s.hi, s.lo);
*/
  return *this;
}

dd_real &dd_real::sub(dd_real const &other)
{
  //dqfun.f90#dqadd is the same as dd_inline.h#dd_real::sloppy_add
  dd_real se;

  se.two_sum(hi, -other.hi);
  se.lo_ += (lo_ - other.lo_);
  quick_two_sum(se.hi, se.lo_);
  return *this;
}

dd_real &dd_real::mul(dd_real const &other)
{
  dd_real p;

  p.two_prod(hi, other.hi);
  p.lo_ += (hi * other.lo_ + lo_ * other.hi);
  quick_two_sum(p.hi, p.lo_);
/*
  double hx, tx, hy, ty, C, c;
  double th, tl;
  C = MAKE_DOUBLE(th, double_double_split*hi);
  hx = MAKE_DOUBLE(tl, C-hi);
  c = MAKE_DOUBLE(th, double_double_split*h2);
  hx = MAKE_DOUBLE(tl, C-hx);
  tx = MAKE_DOUBLE(th, hi-hx);
  hy = MAKE_DOUBLE(tl, c-h2);
  C = MAKE_DOUBLE(th, hi*h2);
  hy = MAKE_DOUBLE(tl, c-hy);
  ty = MAKE_DOUBLE(th, h2-hy);
  c = MAKE_DOUBLE(tl, ((((hx*hy-C)+hx*ty)+tx*hy)+tx*ty)+(hi*l2+lo*h2));

  hi = MAKE_DOUBLE(th, C+c);
  hx = MAKE_DOUBLE(th, C-hi);
  lo = MAKE_DOUBLE(tl, c+hx);
*/
  return *this;
}

dd_real &dd_real::sqr()
{
  dd_real p;
  dd_real s;
  p.two_sqr(hi);
  p.lo_ += 2.0 * hi * lo_;
  p.lo_ += lo_ * lo_;
  quick_two_sum(p.hi, p.lo_);
  return *this;
}

double dd_real::radixfloor() const
{
  int ilog1=std::ilogb(hi);
  return ldexp(1, ilog1);
}

dd_real &dd_real::recip()
{
#if 1 //sloppy_div
  double s1, s2;
  double q1, q2;
  dd_real r;

  q1 = 1 / hi;  /* approximate quotient */

  /* compute  this - q1 * dd */
  r.hi=hi;
  r.lo_=lo_;
  r.mul_double(q1);
  s1 = two_diff(1, r.hi, s2);
  s2 -= r.lo_;
  s2 += 0;

  /* get next approximation */
  q2 = (s1 + s2) / hi;

  /* renormalize */
  quick_two_sum(q1, q2);
#else //accurate_div
  double q1, q2, q3;
  dd_real r, tmp;

  q1 = 1 / hi;  /* approximate quotient */

  r.hi=hi;
  r.lo=lo;
  r.mul(q1, 0);
  r.chs();
  r.add(1, 0); // [1,0] - q1 * self;

  q2 = r.hi / hi;
  tmp.hi=hi;
  tmp.lo=lo;
  tmp.mul(q2, 0);
  r.add(-tmp.hi, -tmp.lo); //  r -= (q2 * self);

  q3 = r.hi / hi;

  quick_two_sum(q1, q2);
  add(q3, 0);
#endif
  return *this;
}

dd_real &dd_real::sqrt()
{
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to 1/sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */

  /*if (hi<0) {
    //dbgPoint();
    hi=0;
    lo_=0;
    return *this;
  };*/

  if (hi<=0)
  {
    hi=0;
    lo_=0;
    return *this;
  };


  double x = 1.0 / std::sqrt(hi);
  double ax = hi * x;
  dd_real result;
  result.hi=ax;
  result.lo_=0;
  result.sqr();
  result.chs();
  add(result);
  hi*=(x*0.5);
  lo_=0;
  add_double(ax);
  //return dd_real::add(ax, (a - dd_real::sqr(ax)).x[0] * (x * 0.5));
  return *this;
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

std::strong_ordering dd_real::compare(dd_real const &other) const
{
  if (hi<other.hi)
    return std::strong_ordering::less;//-1;
  else if (hi>other.hi)
    return std::strong_ordering::greater;
  else if (lo_<other.lo_)
    return std::strong_ordering::less;
  else if (lo_>other.lo_)
    return std::strong_ordering::greater;
  else
    return std::strong_ordering::equal;
}

/*int dd_real::compare(double other_h, double other_l) const
{
  if (hi<other_h)
    return -1;
  else if (hi>other_h)
    return +1;
  else if (lo_<other_l)
    return -1;
  else if (lo_>other_l)
    return +1;
  else
    return 0;
}*/

bool dd_real::isequal(dd_real const &other) const
{
  return (hi==other.hi) && (lo_==other.lo_);
}

bool dd_real::is0() const
{
  return hi==0;
}

bool dd_real::isle(dd_real const &other) const
{
  if (hi!=other.hi)
    return hi<other.hi;
  return lo_<=other.lo_;
}

bool dd_real::isle0() const
{
  return hi<=0;
}

bool dd_real::isl0() const
{
  return hi<0;
}

bool dd_real::isl1() const
{
  return hi<1;
}

dd_real &dd_real::min(dd_real const &other)
{
  if (!isle(other))
  {
    assign(other);
  };
  return *this;
}

QString dd_real::toString() const
{
  return QString("dd(%1,%2)").arg(hi, 0, 'f', 16).arg(lo_, 0, 'g', 16);
}

int dd_real::toRound() const
{
  return std::floor(hi+0.5)+std::floor(lo_+0.5);
}

double dd_real::toDouble() const
{
  return hi;
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
