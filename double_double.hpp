#ifndef MANDELMATH_DOUBLE_DOUBLE_HPP
#define MANDELMATH_DOUBLE_DOUBLE_HPP

#include <bit>
#include <cmath>
#include <compare>
#include <cstdint>
#include <QString>

namespace MandelMath {

void fpu_fix_start(unsigned int *old_cw);
void fpu_fix_end(unsigned int *old_cw);

// Decimal <-> __float128 via libquadmath, always with '.' regardless of LC_NUMERIC.
bool float128_fromString(char const *str, __float128 &result); // false unless the whole string is a number
QString float128_toString(__float128 x, int prec);
QString float128_toString(__float128 x); // shortest that parses back to x

// x*2^n, bit-identical to std::ldexp; avoids the libm call while 2^n is a normal double
inline double ldexp2(double x, int n)
{
  if (n>=-1022 && n<=1023)
    return x*std::bit_cast<double>(uint64_t(1023+n)<<52);
  return std::ldexp(x, n);
}

// 2^ilogb(x) by masking the exponent bits, bit-identical to std::ldexp(1.0, std::ilogb(x)); libm only for 0, subnormal, inf, NaN
inline double radixfloor2(double x)
{
  uint64_t e=std::bit_cast<uint64_t>(x) & 0x7FF0000000000000ull;
  if (e!=0 && e!=0x7FF0000000000000ull)
    return std::bit_cast<double>(e);
  return std::ldexp(1.0, std::ilogb(x));
}

class dd_real {
protected:
  static constexpr double QD_SPLITTER=134217729.0;               // = 2^27 + 1
  static constexpr double QD_SPLIT_THRESH=6.69692879491417e+299; // = 2^996
  void quick_two_sum(double a, double b);
  void two_sum(double a, double b);
  static double two_diff(double a, double b, double &err);
  void two_prod(double a, double b);
  void two_sqr(double a);
public:
  void split(double a);
  dd_real(): hi(0), lo_(0) { }
  dd_real(double h, double l): hi(h), lo_(l) { }
  double hi;
  double lo_; // !!can have other sign than hi!!
  double eps2() const;
  double eps234() const;
  dd_real &assign(const dd_real &src) noexcept { hi=src.hi; lo_=src.lo_; return *this; }
  dd_real &zero(double v=0) { hi=v; lo_=0; return *this; }
  dd_real &chs() { hi=-hi; lo_=-lo_; return *this; }
  dd_real &lshift(int exp); //*=2^exp
  dd_real &add_double(double h2);
  dd_real &mul_double(double h2);
  //dd_real &add(double h2, double l2); //TODO: -> add(dd_real)
  dd_real &add(dd_real const &other);
  dd_real &sub(dd_real const &other);
  dd_real &rsub(dd_real const &other) { chs(); return add(other); }
  //dd_real &mul(double h2, double l2);
  dd_real &mul(dd_real const &other);
  dd_real &sqr();
  double radixfloor() const { return radixfloor2(hi); } //nearest smaller power of 2 (1.5->1->1)
  dd_real &recip();
  dd_real &sqrt();
  dd_real &round();
  dd_real &frac();
  dd_real &mod1();
  bool reduce_angle();
  dd_real &add_pi(double x);
  //int compare(const dd_real *other) const; //return -1 if <, 0 if =, +1 if >
  std::strong_ordering compare(dd_real const &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  //int compare(double other_h, double other_l) const;
  bool isequal(dd_real const &other) const;
  bool is0() const;
  bool isle(const dd_real &other) const;
  bool isle0() const;
  bool isl0() const;
  bool isl1() const;
  dd_real &min(dd_real const &other);

  //explicit operator double();
  QString toString() const; // decimal value of hi+lo_, shortest that parses back to the same dd
  bool fromString(char const *str);
  int toRound() const;
  double toDouble() const;
};

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline void dd_real::quick_two_sum(double a, double b)
{
  hi = a + b;
  lo_ = b + (a - hi);
}

/* Computes fl(a+b) and err(a+b).  */
inline void dd_real::two_sum(double a, double b)
{
  hi = a + b;
  double bb = hi - a;
  lo_ = (a - (hi - bb)) + (b - bb);
}

/* Computes fl(a-b) and err(a-b).  */
inline double dd_real::two_diff(double a, double b, double &err)
{
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

/* Computes high word and lo word of a */
inline void dd_real::split(double a)
{
  double temp;
  if (a > QD_SPLIT_THRESH || a < -QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo_ = a - hi;
    hi *= 268435456.0;       // 2^28
    lo_ *= 268435456.0;      // 2^28
  } else {
    temp = QD_SPLITTER * a;  //a=a1*2^27+a2*1  temp=a1*2^54+(a1+a2)*2^27 +(a2)*1
    hi = temp - (temp - a);  //temp-a=a1*2^54+a2*2^27    temp-(temp-a)=a1*2^27
    lo_ = a - hi;
  }
}

/* Computes fl(a*b) and err(a*b). */
inline void dd_real::two_prod(double a, double b)
{
#ifdef FP_FAST_FMA
  hi = a * b;
  lo_ = std::fma(a, b, -hi);
#else
  dd_real aa;
  dd_real bb;
  hi = a * b;
  aa.split(a);
  bb.split(b);
  lo_ = ((aa.hi * bb.hi - hi) + aa.hi * bb.lo_ + aa.lo_ * bb.hi) + aa.lo_ * bb.lo_;
#endif
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline void dd_real::two_sqr(double a)
{
#ifdef FP_FAST_FMA
  hi = a * a;
  lo_ = std::fma(a, a, -hi);
#else
  dd_real aa;
  hi = a * a;
  aa.split(a);
  lo_ = ((aa.hi * aa.hi - hi) + 2.0 * aa.hi * aa.lo_) + aa.lo_ * aa.lo_;
#endif
}

inline double dd_real::eps2() const
{
  return 6.07717e-64; // 2^-(2*(53+52)) rounded up
}

inline double dd_real::eps234() const
{
  return 3.87e-48; // eps2^(3/4)
}

inline dd_real &dd_real::lshift(int exp)
{
  hi=ldexp2(hi, exp);
  lo_=ldexp2(lo_, exp);
  return *this;
}

inline dd_real &dd_real::add_double(double h2)
{
  dd_real se;
  se.two_sum(hi, h2);
  se.lo_ += lo_;
  quick_two_sum(se.hi, se.lo_);
  return *this;
}

inline dd_real &dd_real::mul_double(double h2)
{
  dd_real p;
  p.two_prod(hi, h2);
#ifdef FP_FAST_FMA
  p.lo_ = std::fma(lo_, h2, p.lo_); // DWTimesFP3, Joldes-Muller-Popescu 2017, error <= 2u^2
#else
  p.lo_ += lo_*h2;
#endif
  quick_two_sum(p.hi, p.lo_);
  return *this;
}

inline dd_real &dd_real::add(dd_real const &other)
{
//dqfun.f90#dqadd is the same as dd_inline.h#dd_real::sloppy_add
  dd_real se;
  se.two_sum(hi, other.hi);
  se.lo_ += (lo_ + other.lo_);
  quick_two_sum(se.hi, se.lo_);
  return *this;
}

inline dd_real &dd_real::sub(dd_real const &other)
{
  dd_real se;
  se.two_sum(hi, -other.hi);
  se.lo_ += (lo_ - other.lo_);
  quick_two_sum(se.hi, se.lo_);
  return *this;
}

inline dd_real &dd_real::mul(dd_real const &other)
{
  dd_real p;
  p.two_prod(hi, other.hi);
#ifdef FP_FAST_FMA
  // DWTimesDW3, Joldes-Muller-Popescu 2017, error <= 4u^2
  double t = std::fma(hi, other.lo_, lo_ * other.lo_);
  p.lo_ += std::fma(lo_, other.hi, t);
#else
  p.lo_ += (hi * other.lo_ + lo_ * other.hi);
#endif
  quick_two_sum(p.hi, p.lo_);
  return *this;
}

inline dd_real &dd_real::sqr()
{
  dd_real p;
  p.two_sqr(hi);
#ifdef FP_FAST_FMA
  p.lo_ += std::fma(2.0 * hi, lo_, lo_ * lo_); // DWTimesDW3 with x==y
#else
  p.lo_ += 2.0 * hi * lo_;
  p.lo_ += lo_ * lo_;
#endif
  quick_two_sum(p.hi, p.lo_);
  return *this;
}

inline dd_real &dd_real::recip()
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

inline dd_real &dd_real::sqrt()
{
  /* Strategy:  Use Karp's trick:  if x is an approximation
     to 1/sqrt(a), then

        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)

     The approximation is accurate to twice the accuracy of x.
     Also, the multiplication (a*x) and [-]*x can be done with
     only half the precision.
  */
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
  return *this;
}

inline std::strong_ordering dd_real::compare(dd_real const &other) const
{
  if (hi<other.hi)
    return std::strong_ordering::less;
  else if (hi>other.hi)
    return std::strong_ordering::greater;
  else if (lo_<other.lo_)
    return std::strong_ordering::less;
  else if (lo_>other.lo_)
    return std::strong_ordering::greater;
  else
    return std::strong_ordering::equal;
}

inline bool dd_real::isequal(dd_real const &other) const
{
  return (hi==other.hi) && (lo_==other.lo_);
}

inline bool dd_real::is0() const
{
  return hi==0;
}

inline bool dd_real::isle(dd_real const &other) const
{
  if (hi!=other.hi)
    return hi<other.hi;
  return lo_<=other.lo_;
}

inline bool dd_real::isle0() const
{
  return hi<=0;
}

inline bool dd_real::isl0() const
{
  return hi<0;
}

inline bool dd_real::isl1() const
{
  return hi<1;
}

inline dd_real &dd_real::min(dd_real const &other)
{
  if (!isle(other))
  {
    assign(other);
  };
  return *this;
}

inline int dd_real::toRound() const
{
  return std::floor(hi+0.5)+std::floor(lo_+0.5);
}

inline double dd_real::toDouble() const
{
  return hi;
}

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_DOUBLE_HPP
