#ifndef MANDELMATH_DOUBLE_EXP_HPP
#define MANDELMATH_DOUBLE_EXP_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <compare>
#include <cstdint>
#include <QString>

#include "double_double.hpp"

namespace MandelMath {

// "floatexp" types: a double / double-double mantissa normalized to 1<=|m|<2 plus an int32_t
// binary exponent.  Inside the double range every operation scales exactly, so results are
// bit-identical to double / dd_real; outside it the range extends to 2^(+-2^30).
namespace exp_detail
{
  constexpr uint64_t EXP_MASK=0x7FF0000000000000ull;
  constexpr uint64_t EXP_ONE=0x3FF0000000000000ull; // exponent field of 1.0
  constexpr int64_t EXP_LIMIT=int64_t(1)<<30;       // larger exponents overflow to inf, smaller underflow to 0
  constexpr int32_t EXP_SPECIAL=INT32_MAX;          // inf and NaN: above every finite exponent
  constexpr int32_t EXP_ZERO=-INT32_MAX;            // zero: below every finite exponent, so add() and compare() need no zero test
  constexpr int32_t LARGE_EXP=int32_t(EXP_LIMIT/2)-16; // isTooLarge() from 2^LARGE_EXP: products of two smaller values stay far below the limit
  inline double pow2(int n) { return std::bit_cast<double>(uint64_t(1023+n)<<52); } // -1022<=n<=1023
}

class de_real
{
  void set_special(double x, int64_t ex); // x is 0, subnormal, inf or NaN
  void set_range(int64_t ex);             // |ex|>EXP_LIMIT
  void set_exp(int64_t ex)
  {
    if (uint64_t(ex+exp_detail::EXP_LIMIT)>uint64_t(2*exp_detail::EXP_LIMIT)) [[unlikely]]
      return set_range(ex);
    e=int32_t(ex);
  }
  de_real &add_signed(double om, int32_t oe)
  {
    using exp_detail::pow2;
    int64_t d=int64_t(e)-oe;
    // beyond 63 bits the smaller operand cannot change the rounded sum
    if (uint64_t(d+63)>126) [[unlikely]]
    {
      if (d<0)
      {
        m=om;
        e=oe;
      }
      return *this;
    }
    // branch-free: which operand is larger is unpredictable; both products are exact
    int sd=int(d);
    set(m*pow2(sd<0?sd:0)+om*pow2(sd>0?-sd:0), sd>=0?e:oe);
    return *this;
  }
public:
  double m;  // 0, inf, NaN or 1<=|m|<2
  int32_t e; // value is m*2^e
  de_real(): m(0), e(exp_detail::EXP_ZERO) { }
  explicit de_real(double x) { set(x, 0); }
  de_real(double x, int64_t ex) { set(x, ex); } // x*2^ex
  void set(double x, int64_t ex)
  {
    uint64_t bits=std::bit_cast<uint64_t>(x);
    uint64_t be=bits & exp_detail::EXP_MASK;
    if ((be>>52)-1>=0x7FE) [[unlikely]] // 0, subnormal, inf or NaN
      return set_special(x, ex);
    m=std::bit_cast<double>((bits & ~exp_detail::EXP_MASK) | exp_detail::EXP_ONE);
    set_exp(ex+int64_t(be>>52)-1023);
  }

  de_real &assign(de_real const &src) noexcept { m=src.m; e=src.e; return *this; }
  de_real &zero(double v=0) { set(v, 0); return *this; }
  de_real &chs() { m=-m; return *this; }
  de_real &lshift(int n)
  {
    if (e>=-exp_detail::EXP_LIMIT && e<=exp_detail::EXP_LIMIT)
      set_exp(int64_t(e)+n);
    return *this;
  }
  de_real &add_double(double x) { return add(de_real(x)); }
  de_real &mul_double(double x) { return mul(de_real(x)); }
  de_real &add(de_real const &o) { return add_signed(o.m, o.e); }
  de_real &sub(de_real const &o) { return add_signed(-o.m, o.e); }
  de_real &rsub(de_real const &o) { m=-m; return add_signed(o.m, o.e); }
  de_real &mul(de_real const &o) { set(m*o.m, int64_t(e)+o.e); return *this; }
  de_real &sqr() { set(m*m, 2*int64_t(e)); return *this; }
  double radixfloor() const { return ldexp2(1.0, e); }
  de_real &recip() { set(1.0/m, -int64_t(e)); return *this; }
  de_real &sqrt()
  {
    if (m<=0 || e==exp_detail::EXP_SPECIAL)
    {
      set(std::sqrt(m), 0);
      return *this;
    }
    int64_t ex=e;
    double mm=m;
    if (ex & 1)
    {
      mm*=2;
      ex--;
    }
    set(std::sqrt(mm), ex/2);
    return *this;
  }
  de_real &round();
  de_real &frac();
  de_real &mod1();
  bool reduce_angle();
  de_real &add_pi(double x) { return add_double(x*3.14159265358979323846); }
  std::strong_ordering compare(de_real const &o) const
  {
    int s=(m>0)-(m<0);
    int os=(o.m>0)-(o.m<0);
    if (s!=os)
      return s<os?std::strong_ordering::less:std::strong_ordering::greater;
    if (s==0)
      return std::strong_ordering::equal;
    if (e!=o.e)
      return ((e<o.e)==(s>0))?std::strong_ordering::less:std::strong_ordering::greater;
    if (m<o.m)
      return std::strong_ordering::less;
    if (m>o.m)
      return std::strong_ordering::greater;
    return std::strong_ordering::equal;
  }
  bool isequal(de_real const &o) const { return m==o.m && e==o.e; }
  bool is0() const { return m==0; }
  bool isle(de_real const &o) const { return compare(o)!=std::strong_ordering::greater; }
  bool isle0() const { return m<=0; }
  bool isl0() const { return m<0; }
  bool isl1() const { return m<1 || e<0; }
  de_real &min(de_real const &o) { if (!isle(o)) assign(o); return *this; }
  bool isTooLarge() const { return e>=exp_detail::LARGE_EXP; }
  static bool isMagTooLarge(de_real const &re, de_real const &im) { return std::max(re.e, im.e)>=exp_detail::LARGE_EXP/2; }

  QString toString() const; // shortest decimal that parses back to the same value
  bool fromString(char const *str);
  int toRound() const;
  double toDouble() const { return ldexp2(m, e); }
  __float128 toFloat128() const;
  static de_real fromFloat128(__float128 x);
};

class dde_real
{
  void set_special(double h, double l, int64_t ex); // h is 0, subnormal, inf or NaN
  void set_range(int64_t ex);                       // |ex|>EXP_LIMIT
  void set_exp(int64_t ex)
  {
    if (uint64_t(ex+exp_detail::EXP_LIMIT)>uint64_t(2*exp_detail::EXP_LIMIT)) [[unlikely]]
      return set_range(ex);
    e=int32_t(ex);
  }
  // past ~2^-960 relative the scaled smaller operand would lose bits to underflow anyway
  static constexpr int ADD_CUTOFF=960;
  dde_real &add_signed(double oh, double ol, int32_t oe)
  {
    int64_t d=int64_t(e)-oe;
    if (uint64_t(d+ADD_CUTOFF)>uint64_t(2*ADD_CUTOFF)) [[unlikely]]
    {
      if (d<0)
      {
        hi=oh;
        lo=ol;
        e=oe;
      }
      return *this;
    }
    // branch-free like de_real::add_signed; two_sum is symmetric so the order does not matter
    int sd=int(d);
    double sa=exp_detail::pow2(sd<0?sd:0);
    double sb=exp_detail::pow2(sd>0?-sd:0);
    dd_real p(hi*sa, lo*sa);
    p.add(dd_real(oh*sb, ol*sb));
    set(p, sd>=0?e:oe);
    return *this;
  }
public:
  double hi; // 0, inf, NaN or 1<=|hi|<2
  double lo; // low word of the dd_real mantissa, scaled like hi
  int32_t e; // value is (hi+lo)*2^e
  dde_real(): hi(0), lo(0), e(exp_detail::EXP_ZERO) { }
  explicit dde_real(double x) { set(x, 0, 0); }
  explicit dde_real(dd_real const &x, int64_t ex=0) { set(x.hi, x.lo_, ex); }
  void set(double h, double l, int64_t ex)
  {
    uint64_t bits=std::bit_cast<uint64_t>(h);
    uint64_t be=bits & exp_detail::EXP_MASK;
    if ((be>>52)-1>=0x7FE) [[unlikely]] // 0, subnormal, inf or NaN
      return set_special(h, l, ex);
    int sh=int(be>>52)-1023;
    hi=std::bit_cast<double>((bits & ~exp_detail::EXP_MASK) | exp_detail::EXP_ONE);
    lo=ldexp2(l, -sh);
    set_exp(ex+sh);
  }
  void set(dd_real const &p, int64_t ex) { set(p.hi, p.lo_, ex); }
  dd_real mantissa() const { return dd_real(hi, lo); }
  dd_real toDD() const { return dd_real(ldexp2(hi, e), ldexp2(lo, e)); }

  dde_real &assign(dde_real const &src) noexcept { hi=src.hi; lo=src.lo; e=src.e; return *this; }
  dde_real &zero(double v=0) { set(v, 0, 0); return *this; }
  dde_real &chs() { hi=-hi; lo=-lo; return *this; }
  dde_real &lshift(int n)
  {
    if (e>=-exp_detail::EXP_LIMIT && e<=exp_detail::EXP_LIMIT)
      set_exp(int64_t(e)+n);
    return *this;
  }
  dde_real &add_double(double x) { return add(dde_real(x)); }
  dde_real &mul_double(double x)
  {
    de_real xx(x);
    dd_real p(hi, lo);
    p.mul_double(xx.m);
    set(p, int64_t(e)+xx.e);
    return *this;
  }
  dde_real &add(dde_real const &o) { return add_signed(o.hi, o.lo, o.e); }
  dde_real &sub(dde_real const &o) { return add_signed(-o.hi, -o.lo, o.e); }
  dde_real &rsub(dde_real const &o) { chs(); return add_signed(o.hi, o.lo, o.e); }
  dde_real &mul(dde_real const &o)
  {
    dd_real p(hi, lo);
    p.mul(dd_real(o.hi, o.lo));
    set(p, int64_t(e)+o.e);
    return *this;
  }
  dde_real &sqr()
  {
    dd_real p(hi, lo);
    p.sqr();
    set(p, 2*int64_t(e));
    return *this;
  }
  double radixfloor() const { return ldexp2(1.0, e); }
  dde_real &recip()
  {
    if (hi==0 || e==exp_detail::EXP_SPECIAL)
    {
      set(1.0/hi, 0, 0);
      return *this;
    }
    dd_real p(hi, lo);
    p.recip();
    set(p, -int64_t(e));
    return *this;
  }
  dde_real &sqrt()
  {
    if (e==exp_detail::EXP_SPECIAL)
    {
      set(std::sqrt(hi), 0, 0);
      return *this;
    }
    if (hi<=0) // same as dd_real::sqrt
    {
      *this=dde_real();
      return *this;
    }
    int64_t ex=e;
    dd_real p(hi, lo);
    if (ex & 1)
    {
      p.hi*=2;
      p.lo_*=2;
      ex--;
    }
    p.sqrt();
    set(p, ex/2);
    return *this;
  }
  dde_real &round();
  dde_real &frac();
  dde_real &mod1();
  bool reduce_angle();
  dde_real &add_pi(double x);
  std::strong_ordering compare(dde_real const &o) const
  {
    int s=(hi>0)-(hi<0);
    int os=(o.hi>0)-(o.hi<0);
    if (s!=os)
      return s<os?std::strong_ordering::less:std::strong_ordering::greater;
    if (s==0)
      return std::strong_ordering::equal;
    if (e!=o.e)
      return ((e<o.e)==(s>0))?std::strong_ordering::less:std::strong_ordering::greater;
    if (hi<o.hi)
      return std::strong_ordering::less;
    if (hi>o.hi)
      return std::strong_ordering::greater;
    if (lo<o.lo)
      return std::strong_ordering::less;
    if (lo>o.lo)
      return std::strong_ordering::greater;
    return std::strong_ordering::equal;
  }
  bool isequal(dde_real const &o) const { return hi==o.hi && lo==o.lo && e==o.e; }
  bool is0() const { return hi==0; }
  bool isle(dde_real const &o) const { return compare(o)!=std::strong_ordering::greater; }
  bool isle0() const { return hi<=0; }
  bool isl0() const { return hi<0; }
  bool isl1() const { return hi<1 || e<0; }
  dde_real &min(dde_real const &o) { if (!isle(o)) assign(o); return *this; }
  bool isTooLarge() const { return e>=exp_detail::LARGE_EXP; }
  static bool isMagTooLarge(dde_real const &re, dde_real const &im) { return std::max(re.e, im.e)>=exp_detail::LARGE_EXP/2; }

  QString toString() const; // shortest decimal that parses back to the same value
  bool fromString(char const *str);
  int toRound() const;
  double toDouble() const { return ldexp2(hi, e); }
  __float128 toFloat128() const;
  static dde_real fromFloat128(__float128 x);
};

} // namespace MandelMath

#endif // MANDELMATH_DOUBLE_EXP_HPP
