#include "double_exp.hpp"
#include <cctype>
#include <charconv>
#include <cstring>
#include <limits>
#include <string>
#include <quadmath.h>
#include <QtGlobal>

namespace MandelMath {

using namespace exp_detail;

namespace {

constexpr double PI_HI=7074237752028440.0/2251799813685248.0; // dd_real pi, as in dd_real::add_pi
constexpr double PI_LO=2483878800010755.0/20282409603651670423947251286016.0;

// m*2^e with 1<=|m|<2: binary128 precision over an int64 exponent range, for decimal conversion
struct wide_real
{
  __float128 m;
  int64_t e;
  void normalize()
  {
    if (m==0 || !finiteq(m))
      return;
    int k=ilogbq(m);
    m=scalbnq(m, -k);
    e+=k;
  }
};

wide_real wide_mul(wide_real a, wide_real const &b)
{
  a.m*=b.m;
  a.e+=b.e;
  a.normalize();
  return a;
}

// 10^k by binary powering; relative error about 2*log2(|k|) binary128 ulps
wide_real wide_pow10(int64_t k)
{
  wide_real result{1, 0};
  wide_real base{10, 0};
  base.normalize();
  for (uint64_t n=k<0?-uint64_t(k):uint64_t(k); n; n>>=1)
  {
    if (n & 1)
      result=wide_mul(result, base);
    base=wide_mul(base, base);
  }
  if (k<0)
  {
    result.m=1/result.m;
    result.e=-result.e;
    result.normalize();
  }
  return result;
}

// Decimal exponents in this range go through binary128 directly (its range is about 1e+-4932).
constexpr int64_t DIRECT_DEC_EXP=4000;
constexpr int64_t DIRECT_BIN_EXP=13000;
// Decimal exponents are clamped here; 10^(+-2e9) is far outside EXP_LIMIT anyway.
constexpr int64_t MAX_DEC_EXP=2000000000;

// Parses str as m*2^e; false unless the whole string is a number.
bool wide_fromString(char const *str, wide_real &result)
{
  // [space][sign]digits[.digits][(e|E)[sign]digits][space]; anything else (inf, nan, hex) goes to libquadmath
  char const *p=str;
  while (std::isspace((unsigned char)*p))
    p++;
  if (*p=='+' || *p=='-')
    p++;
  int digits=0;
  while (std::isdigit((unsigned char)*p))
    p++, digits++;
  if (*p=='.')
  {
    p++;
    while (std::isdigit((unsigned char)*p))
      p++, digits++;
  }
  char const *mant_end=p;
  int64_t dec_exp=0;
  bool decimal=digits>0;
  if (decimal && (*p=='e' || *p=='E'))
  {
    p++;
    bool neg=false;
    if (*p=='+' || *p=='-')
      neg=(*p++=='-');
    if (!std::isdigit((unsigned char)*p))
      decimal=false;
    while (std::isdigit((unsigned char)*p))
    {
      if (dec_exp<MAX_DEC_EXP)
        dec_exp=dec_exp*10+(*p-'0');
      p++;
    }
    if (neg)
      dec_exp=-dec_exp;
  }
  char const *rest=p;
  while (std::isspace((unsigned char)*rest))
    rest++;
  if (!decimal || *rest!=0 || (dec_exp>=-DIRECT_DEC_EXP && dec_exp<=DIRECT_DEC_EXP))
  {
    result.e=0;
    if (!float128_fromString(str, result.m))
      return false;
    result.normalize();
    return true;
  }
  std::string mant(str, mant_end);
  if (!float128_fromString(mant.c_str(), result.m))
    return false;
  result.e=0;
  result.normalize();
  if (result.m!=0)
    result=wide_mul(result, wide_pow10(dec_exp));
  return true;
}

// Like float128_toString(x, prec) but for any exponent.
std::string wide_toString(wide_real x, int prec)
{
  if (x.m==0 || !finiteq(x.m) || (x.e>-DIRECT_BIN_EXP && x.e<DIRECT_BIN_EXP))
    return float128_toString(scalbnq(x.m, int(x.e)), prec).toUtf8().constData();
  // value = y*10^k with y near [1, 10)
  int64_t k=(int64_t)std::floor((double(x.e)+std::log2(std::fabs((double)x.m)))*0.30102999566398119521);
  wide_real y=wide_mul(x, wide_pow10(-k));
  std::string s=float128_toString(scalbnq(y.m, int(y.e)), prec).toUtf8().constData();
  size_t epos=s.find_first_of("eE");
  if (epos!=std::string::npos)
  {
    k+=std::strtoll(s.c_str()+epos+1, nullptr, 10);
    s.resize(epos);
  }
  return s+"e"+std::to_string(k);
}

// std::from_chars for a finite normal double, accepting the same blanks and '+' as number_traits<double>::fromString
bool double_fromString(char const *str, double &result)
{
  while (std::isspace((unsigned char)*str))
    str++;
  if (*str=='+')
    str++;
  char const *end=str+std::strlen(str);
  auto [ptr, ec]=std::from_chars(str, end, result);
  if (ec!=std::errc())
    return false;
  while (std::isspace((unsigned char)*ptr))
    ptr++;
  return ptr==end && std::isnormal(result);
}

} // namespace


void de_real::set_special(double x, int64_t ex)
{
  if (x==0)
  {
    m=x;
    e=EXP_ZERO;
  }
  else if (!std::isfinite(x))
  {
    m=x;
    e=EXP_SPECIAL;
  }
  else
  {
    int q;
    double f=std::frexp(x, &q); // subnormal: 0.5<=|f|<1
    set(2*f, ex+q-1);
  }
}

void de_real::set_range(int64_t ex)
{
  if (ex>0)
  {
    m=std::copysign(std::numeric_limits<double>::infinity(), m);
    e=EXP_SPECIAL;
  }
  else
  {
    m=std::copysign(0.0, m);
    e=EXP_ZERO;
  }
}

de_real &de_real::round()
{
  if (e>=-1022 && e<=1023)
    set(std::round(toDouble()), 0);
  else if (e<0) // |x|<0.5
    set(std::copysign(0.0, m), 0);
  // else an integer, inf or NaN
  return *this;
}

de_real &de_real::frac()
{
  if (e>=-1022 && e<=1023)
  {
    double s=toDouble();
    if (s<0)
      s-=std::ceil(s);
    else
      s-=std::floor(s);
    set(s, 0);
  }
  else if (e>0)
    set(m-m, 0); // integer -> 0, inf -> NaN
  return *this;
}

de_real &de_real::mod1()
{
  if (e>=-1022 && e<=1023)
  {
    double s=toDouble();
    set(s-std::floor(s), 0);
  }
  else if (e>0)
    set(m-m, 0);
  else if (m<0)
    add_double(1);
  return *this;
}

bool de_real::reduce_angle()
{
  double s=toDouble();
  if (s<-3.14159265358979323846)
  {
    add_double(2*3.14159265358979323846);
    return true;
  }
  else if (s>=3.14159265358979323846)
  {
    add_double(-2*3.14159265358979323846);
    return true;
  }
  return false;
}

int de_real::toRound() const
{
  return qRound(toDouble());
}

__float128 de_real::toFloat128() const
{
  return scalbnq(m, e);
}

de_real de_real::fromFloat128(__float128 x)
{
  wide_real w{x, 0};
  w.normalize();
  return de_real((double)w.m, w.e);
}

QString de_real::toString() const
{
  char buf[32];
  if ((e>=-1022 && e<=1023) || e==EXP_ZERO || e==EXP_SPECIAL)
  {
    *std::to_chars(buf, buf+sizeof(buf)-1, toDouble()).ptr=0;
    return QString(buf);
  }
  std::string s;
  for (int prec=1; prec<=17; prec++)
  {
    s=wide_toString({m, e}, prec);
    de_real back;
    if (back.fromString(s.c_str()) && back.isequal(*this))
      break;
  }
  return QString::fromStdString(s);
}

bool de_real::fromString(char const *str)
{
  double d;
  if (double_fromString(str, d)) // correctly rounded, like number_traits<double>::fromString
  {
    set(d, 0);
    return true;
  }
  wide_real w;
  if (!wide_fromString(str, w))
    return false;
  set((double)w.m, w.e);
  return true;
}


void dde_real::set_special(double h, double l, int64_t ex)
{
  if (h==0)
  {
    hi=h;
    lo=0;
    e=EXP_ZERO;
  }
  else if (!std::isfinite(h))
  {
    hi=h;
    lo=0;
    e=EXP_SPECIAL;
  }
  else
  {
    int q;
    double f=std::frexp(h, &q);
    set(2*f, std::ldexp(l, 1-q), ex+q-1);
  }
}

void dde_real::set_range(int64_t ex)
{
  lo=0;
  if (ex>0)
  {
    hi=std::copysign(std::numeric_limits<double>::infinity(), hi);
    e=EXP_SPECIAL;
  }
  else
  {
    hi=std::copysign(0.0, hi);
    e=EXP_ZERO;
  }
}

// e range in which toDD() is exact (lo stays normal) and integers are below 2^1000
static bool dde_fits_dd(int32_t e)
{
  return e>=-960 && e<=1000;
}

dde_real &dde_real::round()
{
  if (dde_fits_dd(e))
  {
    dd_real p=toDD();
    p.round();
    set(p, 0);
  }
  else if (e<0)
    set(std::copysign(0.0, hi), 0, 0);
  return *this;
}

dde_real &dde_real::frac()
{
  if (dde_fits_dd(e))
  {
    dd_real p=toDD();
    p.frac();
    set(p, 0);
  }
  else if (e>0)
    set(hi-hi, 0, 0);
  return *this;
}

dde_real &dde_real::mod1()
{
  if (dde_fits_dd(e))
  {
    dd_real p=toDD();
    p.mod1();
    set(p, 0);
  }
  else if (e>0)
    set(hi-hi, 0, 0);
  else if (hi<0)
    add_double(1);
  return *this;
}

bool dde_real::reduce_angle()
{
  static const dde_real pi(dd_real(PI_HI, PI_LO));
  static const dde_real minus_pi(dd_real(-PI_HI, -PI_LO));
  static const dde_real two_pi(dd_real(2*PI_HI, 2*PI_LO));
  if (compare(minus_pi)==std::strong_ordering::less)
  {
    add(two_pi);
    return true;
  }
  else if (compare(pi)!=std::strong_ordering::less)
  {
    sub(two_pi);
    return true;
  }
  return false;
}

dde_real &dde_real::add_pi(double x)
{
  dde_real pi(dd_real(PI_HI, PI_LO));
  pi.mul_double(x);
  return add(pi);
}

int dde_real::toRound() const
{
  if (dde_fits_dd(e))
    return toDD().toRound();
  return qRound(toDouble());
}

__float128 dde_real::toFloat128() const
{
  return scalbnq((__float128)hi+lo, e);
}

dde_real dde_real::fromFloat128(__float128 x)
{
  wide_real w{x, 0};
  w.normalize();
  double h=(double)w.m;
  dde_real result;
  result.set(h, (double)(w.m-h), w.e);
  return result;
}

QString dde_real::toString() const
{
  if (e==EXP_ZERO || e==EXP_SPECIAL)
    return dd_real(hi, lo).toString();
  wide_real w{(__float128)hi+lo, e}; // exact, as in dd_real::toString
  std::string s;
  for (int prec=31; prec<=36; prec++)
  {
    s=wide_toString(w, prec);
    dde_real back;
    if (back.fromString(s.c_str()) && back.isequal(*this))
      break;
  }
  return QString::fromStdString(s);
}

bool dde_real::fromString(char const *str)
{
  wide_real w;
  if (!wide_fromString(str, w))
    return false;
  // same split as dd_real::fromString
  double h=(double)w.m;
  set(h, (double)(w.m-h), w.e);
  return true;
}

} // namespace MandelMath
