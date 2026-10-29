#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <tuple>
#include <type_traits>
#include <variant>
#include <bit>
#include <compare>
#include <QString>
#include <cassert>

#include "double_double.hpp"
#include "dd_heap.hpp"
//#include "multiprec.hpp"

#define NUMBER_DOUBLE_EXISTS 1
#define NUMBER_DOUBLE_ONLY 0

#ifndef M_PI
#define M_PI 3.14159265358979323846
//constexpr double M_PI=4*atan2(1., 0);//acos(-1);
#endif

void dbgPoint();
void dbgPointII(int x1, int x2);
#if QT_NO_DEBUG
#define working_assert(x)
#else
//#define working_assert(x) { if (!(x)) dbgPoint(); }
#define working_assert(x) assert(x)
#endif

template <typename Tout, typename Tin> //typechecked hard cast
Tout specific_cast(Tin in) { return (Tout)in; }

namespace MandelMath {

int gcd(int m, int n);
int ctz16(int x);
template <int total, int block> int ReverseBits(int val);
bool is2kof(int big, int small); //big == small*2^k

template <typename E, E ...vals>
bool enum_is_one_of(E value)
{
  return ((value==vals) || ...);
}

//void atomic_min(std::atomic<int> &a, int val) noexcept;
//void atomic_max(std::atomic<int> &a, int val) noexcept;

template<int N, typename... Ts> using NthTypeOf =
    typename std::tuple_element<N, std::tuple<Ts...>>::type;
template<typename ...Options>
auto visit(auto callable, std::variant<Options...> &vari, auto &...params)
{
  //using cases=std::make_index_sequence<sizeof...(Options)>;
  auto index=vari.index(); //size_t
  if constexpr (sizeof...(Options)>0)
    if (index==0) return callable((NthTypeOf<0, Options...> &)vari, params...); //__detail::__variant::__get<_Np>(__v);
  if constexpr (sizeof...(Options)>1)
    //if (index==1) return callable(std::get<1>(vari), params...);
    if (index==1) return callable((NthTypeOf<1, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>2)
    if (index==2) return callable((NthTypeOf<2, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>3)
    if (index==3) return callable((NthTypeOf<3, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>4)
    if (index==4) return callable((NthTypeOf<4, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>5)
    if (index==5) return callable((NthTypeOf<5, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>6)
    if (index==6) return callable((NthTypeOf<6, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>7)
    if (index==7) return callable((NthTypeOf<7, Options...> &)vari, params...);
  static_assert(sizeof...(Options)<=8, "extend visit");
  using RT=decltype(callable((NthTypeOf<0, Options...> &)vari, params...));
  [[unlikely]]
  if constexpr (std::is_same_v<RT, std::strong_ordering>)
    return std::strong_ordering::equal; //even more failure from C++ commitee
  else
    return RT();
}

template<typename ...Options>
auto visit(auto callable, std::variant<Options...> const &vari, auto &...params)
{
  //using cases=std::make_index_sequence<sizeof...(Options)>;
  auto index=vari.index(); //size_t
  if constexpr (sizeof...(Options)>0)
    if (index==0) return callable((NthTypeOf<0, Options...> &)vari, params...); //__detail::__variant::__get<_Np>(__v);
  if constexpr (sizeof...(Options)>1)
    //if (index==1) return callable(std::get<1>(vari), params...);
    if (index==1) return callable((NthTypeOf<1, Options...> const &)vari, params...);
  if constexpr (sizeof...(Options)>2)
    if (index==2) return callable((NthTypeOf<2, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>3)
    if (index==3) return callable((NthTypeOf<3, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>4)
    if (index==4) return callable((NthTypeOf<4, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>5)
    if (index==5) return callable((NthTypeOf<5, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>6)
    if (index==6) return callable((NthTypeOf<6, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>7)
    if (index==7) return callable((NthTypeOf<7, Options...> &)vari, params...);
  static_assert(sizeof...(Options)<=8, "extend visit");
  using RT=decltype(callable((NthTypeOf<0, Options...> &)vari, params...));
  [[unlikely]]
  if constexpr (std::is_same_v<RT, std::strong_ordering>)
    return std::strong_ordering::equal; //even more failure from C++ commitee
  else
    return RT();
}

//never seen before... yet another failure from the C++ commitee
static_assert(sizeof(std::strong_ordering::less)==1, "change the following cast");
constexpr int strong_ordering_cast(std::strong_ordering val) { return std::bit_cast<int8_t>(val); }
typedef enum
{
  less=strong_ordering_cast(std::strong_ordering::less),
  equal=strong_ordering_cast(std::strong_ordering::equal),
  greater=strong_ordering_cast(std::strong_ordering::greater),
} strong_ordering;

enum NumberType { typeEmpty
#if NUMBER_DOUBLE_EXISTS
            , typeDouble
#endif
#if !NUMBER_DOUBLE_ONLY
            , typeFloat128
            , typeDDouble
            , typeQDouble
            , typeReal642
            , typeDDHeap
#endif
          };

//typedef dd_real dq_real;
class dq_real: public dd_real
{
  public:
  dq_real(): dd_real() {} //hi(0), lo_(0) { }
  dq_real(double h, double l): dd_real(h, l) {} //hi(h), lo_(l) { }
};

struct real642
{
  double val1;
  __float128 val2;
};

template <typename BASE> struct NumberTypeFromBase { static constexpr NumberType ntype=typeEmpty; };
template <> struct NumberTypeFromBase<double> { static constexpr NumberType ntype=typeDouble; };
#if !NUMBER_DOUBLE_ONLY
template <> struct NumberTypeFromBase<__float128> { static constexpr NumberType ntype=typeFloat128; };
template <> struct NumberTypeFromBase<dd_real> { static constexpr NumberType ntype=typeDDouble; };
template <> struct NumberTypeFromBase<dq_real> { static constexpr NumberType ntype=typeQDouble; };
template <> struct NumberTypeFromBase<real642> { static constexpr NumberType ntype=typeReal642; };
template <> struct NumberTypeFromBase<dd_heap> { static constexpr NumberType ntype=typeDDHeap; };
#endif

extern int CreatedInstancesOfNumber;

template<typename BASE>
class number;

template<typename BASE>
class complex;

using number_any=std::variant<std::monostate, number<double>, number<__float128>, number<dd_real>, number<dq_real>, number<real642>, number<dd_heap>>;

template<typename T>
concept IMandelNumber = requires(T a, T b)
{
  //no worky using namespace std;
  //simple way:
  //{ a.add(b) } -> std::same_as<T>;
  //strict way:
  //requires std::same_as<decltype ( &T::add ), T &(T::*)(T const &)>;

  { a.eps2() } -> std::same_as<double>;
  requires std::same_as<decltype ( &T::eps234 ), double(T::*)() const>;
  requires std::same_as<decltype ( &T::sqrtlarge ), double(T::*)() const>;

  //requires std::same_as<decltype ( &T::zero ), T &(T::*)(double val)>;
  requires std::same_as<decltype ( static_cast<T &(T::*)(double x)>(&T::zero) ), T &(T::*)(double)>;
  requires std::same_as<decltype ( &T::assign ), T &(T::*)(T const &) noexcept>;
  //will solve overloads later...maybe  requires std::same_as<decltype ( &T::assign_across ), T &(T::*)(number_a const &)>;
  requires std::same_as<decltype ( &T::chs ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::lshift ), T &(T::*)(int shoft)>;
  requires std::same_as<decltype ( &T::add_double ), T &(T::*)(double val)>;
  requires std::same_as<decltype ( &T::mul_double ), T &(T::*)(double val)>;
  requires std::same_as<decltype ( &T::add ), T &(T::*)(T const &)>;
  requires std::same_as<decltype ( &T::sub ), T &(T::*)(T const &)>;
  requires std::same_as<decltype ( &T::rsub ), T &(T::*)(T const &)>;
  //requires std::same_as<decltype ( &T::mul ), T &(T::*)(T const &)>;
  //select the right obverload
  requires std::same_as<decltype ( static_cast<T &(T::*)(T const &)>(&T::mul) ), T &(T::*)(T const &)>;
  requires std::same_as<decltype ( &T::sqr ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::recip ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::sqrt ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::compare ), std::strong_ordering (T::*)(T const &other) const>; //operator<=>
  requires std::same_as<decltype ( &T::isequal ), bool (T::*)(T const &other) const>; //operator==
  requires std::same_as<decltype ( &T::is0 ), bool (T::*)() const>; // ==0
  requires std::same_as<decltype ( &T::isle ), bool(T::*)(T const &other) const>; //operator<=
  requires std::same_as<decltype ( &T::isle0 ), bool (T::*)() const>; // <=0
  requires std::same_as<decltype ( &T::isl0 ), bool (T::*)() const>; // <0

  requires std::same_as<decltype ( &T::toString ), QString (T::*)() const>;

  /* maybe later
  virtual NumberType ntype() const=0;// { return this->_ntype; }
  virtual NumberType raw_ntype() const=0;// { return this->_ntype; }
  virtual void readFrom(void *storage, int index)=0; //BASE storage[index]
  virtual void writeTo(void *storage, int index) const=0; //BASE storage[index]
*/
  requires std::is_nothrow_assignable_v<T, T>;
  requires std::is_nothrow_copy_assignable_v<T>;
  requires std::is_nothrow_move_constructible_v<T>;
  requires std::is_nothrow_copy_constructible_v<T>;
  requires std::is_move_assignable_v<T>;
  //requires std::is_trivially_copyable_v<T>;
};

template<typename T>
concept IMandelReal = requires(T a, T b)
{
  requires IMandelNumber<T>;
  requires std::same_as<decltype ( &T::round ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::frac ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::mod1 ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::radixfloor ), double (T::*)() const>;
  requires std::same_as<decltype ( &T::reduce_angle ), bool (T::*)()>;
  requires std::same_as<decltype ( &T::add_pi ), T &(T::*)(double val)>;
  requires std::same_as<decltype ( &T::isl1 ), bool (T::*)() const>; // <1
  requires std::same_as<decltype ( &T::min ), T &(T::*)(T const &)>;

  requires std::same_as<decltype ( &T::toRound ), int (T::*)() const>;
  requires std::same_as<decltype ( &T::toDouble ), double (T::*)() const>;
};

template<typename T>
concept IMandelComplex = requires(T a, T b)
{
  requires IMandelNumber<T>;
};

inline double sqr_double(double x) //no one ever needed this function before year 2022, right?
{
  return x*x;
}

// s*2^n by adjusting the exponent bits; memcpy keeps it valid under strict aliasing once inlined
inline void lshift_float128(__float128 &s, int n)
{
  if (s != 0) //don't play with 0's exponent or you get NaN
  {
    constexpr size_t top=sizeof(s)-sizeof(uint16_t); // sign+exponent word, little-endian
    uint16_t exponent;
    std::memcpy(&exponent, (char *)&s+top, sizeof(exponent));
    exponent+=n;
    //TODO: should check for overflow/underflow but unlikely to happen
    std::memcpy((char *)&s+top, &exponent, sizeof(exponent));
  }
}

// Type-specific operations on the raw store of number<BASE>; one specialisation per concrete BASE.
template<typename BASE>
struct number_traits;

template<>
struct number_traits<double>
{
  static double eps2() { return 1.232596e-32; } // 2^-(2*53) rounded up
  static double eps234() { return 1.17e-24; } // eps2^(3/4)
  static double sqrtlarge() { return 1.915e153; } // sqrt(DBL_MAX)/7
  static void zero(double &s, double v) { s=v; }
  static void chs(double &s) { s=-s; }
  static void lshift(double &s, int n) { s=ldexp2(s, n); }
  static void round(double &s) { s=std::round(s); }
  static void frac(double &s)
  {
    if (s<0)
      s -= std::ceil(s);
    else
      s -= std::floor(s);
  }
  static void mod1(double &s) { s -= std::floor(s); }
  static void add_double(double &s, double x) { s += x; }
  static void mul_double(double &s, double x) { s *= x; }
  static void add(double &s, double o) { s += o; }
  static void sub(double &s, double o) { s -= o; }
  static void rsub(double &s, double o) { s = o-s; }
  static void mul(double &s, double o) { s *= o; }
  static void sqr(double &s) { s *= s; }
  static double radixfloor(double s) { return radixfloor2(s); }
  static void recip(double &s) { s = 1.0/s; }
  static void sqrt(double &s) { s = std::sqrt(s); }
  static bool reduce_angle(double &s)
  {
    if (s<-M_PI)
    {
      s+=2*M_PI;
      return true;
    }
    else if (s>=M_PI)
    {
      s-=2*M_PI;
      return true;
    }
    else
      return false;
  }
  static void add_pi(double &s, double x) { s+=x*M_PI; }
  static std::strong_ordering compare(double a, double b)
  {
    if (a == b)
      return std::strong_ordering::equal;
    else if (a < b)
      return std::strong_ordering::less;
    else
      return std::strong_ordering::greater;
  }
  static bool isequal(double a, double b) { return a==b; }
  static bool is0(double s) { return s==0; }
  static bool isle(double a, double b) { return a<=b; }
  static bool isle0(double s) { return s<=0; }
  static bool isl0(double s) { return s<0; }
  static bool isl1(double s) { return s<1; }
  static void min(double &s, double o) { if (s>o) s=o; }
  static QString toString(double store);
  static bool fromString(double &s, char const *str);
  static int toRound(double store);
  static double toDouble(double s) { return s; }
  static double *convert_block(NumberType old_type, const void *old_data, int count);
};

template<>
struct number_traits<__float128>
{
  static constexpr __float128 M_PIq=3.141592653589793238462643383279502884197Q;
  static double eps2() { return 9.27302e-69; } // 2^-(2*113) rounded up
  static double eps234() { return 9.45e-52; } // eps2^(3/4)
  static double sqrtlarge() { return 1e300; } // sqrt(FLT128_MAX)/7 does not fit; limits are compared in double
  static void zero(__float128 &s, double v) { s=v; }
  static void chs(__float128 &s) { s=-s; }
  static void lshift(__float128 &s, int n) { lshift_float128(s, n); }
  static void round(__float128 &s)
  {
    __float128 remainder=s;
    s=0;
    for (;;)
    {
      double some=std::round((double)remainder);
      if (some==0)
        break;
      s+=some;
      remainder-=some;
    }
  }
  static void frac(__float128 &s)
  {
    if (s<0)
    {
      double some=std::ceil((double)s);
      while (some)
      {
        s -= some;
        some=std::ceil((double)s);
      }
    }
    else
    {
      double some=std::floor((double)s);
      while (some)
      {
        s -= some;
        some=std::floor((double)s);
      }
    }
  }
  static void mod1(__float128 &s)
  {
    double some=std::floor((double)s);
    while (some)
    {
      s -= some;
      some=std::floor((double)s);
    }
  }
  static void add_double(__float128 &s, double x) { s += x; }
  static void mul_double(__float128 &s, double x) { s *= x; }
  static void add(__float128 &s, __float128 o) { s += o; }
  static void sub(__float128 &s, __float128 o) { s -= o; }
  static void rsub(__float128 &s, __float128 o) { s = o-s; }
  static void mul(__float128 &s, __float128 o) { s *= o; }
  static void sqr(__float128 &s) { s *= s; }
  static double radixfloor(__float128 s) { return radixfloor2((double)s); }
  static void recip(__float128 &s) { s = __float128(1.0)/s; }
  static void sqrt(__float128 &s)
  {
    double x=1.0/std::sqrt((double)s);
    __float128 ax=x*(double)s;
    double tmp=s - ax*ax;
    double tmp2=tmp*x/2;
    //one more iteration; error >30 eps with just one because float128 is 112 bits but 2*double is 104 bits
    ax+=tmp2;
    tmp=s - ax*ax;
    tmp2=tmp*x/2;
    s = ax + tmp2;
  }
  static bool reduce_angle(__float128 &s)
  {
    if (s<-M_PIq)
    {
      s+=2*M_PIq;
      return true;
    }
    else if (s>=M_PIq)
    {
      s-=2*M_PIq;
      return true;
    }
    else
      return false;
  }
  static void add_pi(__float128 &s, double x) { s+=x*M_PIq; }
  static std::strong_ordering compare(__float128 a, __float128 b)
  {
    if (a == b)
      return std::strong_ordering::equal;
    else if (a < b)
      return std::strong_ordering::less;
    else
      return std::strong_ordering::greater;
  }
  static bool isequal(__float128 a, __float128 b) { return a==b; }
  static bool is0(__float128 s) { return s==0; }
  static bool isle(__float128 a, __float128 b) { return a<=b; }
  static bool isle0(__float128 s) { return s<=0; }
  static bool isl0(__float128 s) { return s<0; }
  static bool isl1(__float128 s) { return s<1; }
  static void min(__float128 &s, __float128 o) { if (s>o) s=o; }
  static QString toString(__float128 store) { return float128_toString(store); }
  static bool fromString(__float128 &s, char const *str) { return float128_fromString(str, s); }
  static int toRound(__float128 store);
  static double toDouble(__float128 s) { return (double)s; }
  static __float128 *convert_block(NumberType old_type, const void *old_data, int count);
};

template<>
struct number_traits<dd_real>
{
  static double eps2() { return 6.07717e-64; } // 2^-(2*(53+52)) rounded up
  static double eps234() { return 3.87e-48; } // eps2^(3/4)
  static double sqrtlarge() { return 1.915e153; } // sqrt(DBL_MAX)/7
  static void zero(dd_real &s, double v) { s.zero(v); }
  static void chs(dd_real &s) { s.chs(); }
  static void lshift(dd_real &s, int n) { s.lshift(n); }
  static void round(dd_real &s) { s.round(); }
  static void frac(dd_real &s) { s.frac(); }
  static void mod1(dd_real &s) { s.mod1(); }
  static void add_double(dd_real &s, double x) { s.add_double(x); }
  static void mul_double(dd_real &s, double x) { s.mul_double(x); }
  static void add(dd_real &s, dd_real const &o) { s.add(o); }
  static void sub(dd_real &s, dd_real const &o) { s.sub(o); }
  static void rsub(dd_real &s, dd_real const &o) { s.rsub(o); }
  static void mul(dd_real &s, dd_real const &o) { s.mul(o); }
  static void sqr(dd_real &s) { s.sqr(); }
  static double radixfloor(dd_real const &s) { return s.radixfloor(); }
  static void recip(dd_real &s) { s.recip(); }
  static void sqrt(dd_real &s) { s.sqrt(); }
  static bool reduce_angle(dd_real &s) { return s.reduce_angle(); }
  static void add_pi(dd_real &s, double x) { s.add_pi(x); }
  static std::strong_ordering compare(dd_real const &a, dd_real const &b) { return a.compare(b); }
  static bool isequal(dd_real const &a, dd_real const &b) { return a.isequal(b); }
  static bool is0(dd_real const &s) { return s.is0(); }
  static bool isle(dd_real const &a, dd_real const &b) { return a.isle(b); }
  static bool isle0(dd_real const &s) { return s.isle0(); }
  static bool isl0(dd_real const &s) { return s.isl0(); }
  static bool isl1(dd_real const &s) { return s.isl1(); }
  static void min(dd_real &s, dd_real const &o) { if (!s.isle(o)) s.assign(o); }
  static QString toString(dd_real const &s) { return s.toString(); }
  static bool fromString(dd_real &s, char const *str) { return s.fromString(str); }
  static int toRound(dd_real const &s) { return s.toRound(); }
  static double toDouble(dd_real const &s) { return s.hi; }
  static dd_real *convert_block(NumberType old_type, const void *old_data, int count);
};

template<>
struct number_traits<dq_real>: number_traits<dd_real>
{
  static dq_real *convert_block(NumberType old_type, const void *old_data, int count);
};

template<>
struct number_traits<real642>
{
  static double eps2() { return 1.232596e-32; } // 2^-(2*53) rounded up
  static double eps234() { return 1.17e-24; } // eps2^(3/4)
  static double sqrtlarge() { return 1.915e153; } // sqrt(DBL_MAX)/7, limited by val1
  static void zero(real642 &s, double v) { s.val1=v; s.val2=v; }
  static void chs(real642 &s) { s.val1=-s.val1; s.val2=-s.val2; }
  static void lshift(real642 &s, int n)
  {
    s.val1=ldexp2(s.val1, n);
    lshift_float128(s.val2, n);
  }
  static void round(real642 &s)
  {
    s.val1=std::round(s.val1);
    __float128 remainder=s.val2;
    s.val2=0;
    for (;;)
    {
      double some=std::round((double)remainder);
      if (some==0)
        break;
      s.val2+=some;
      remainder-=some;
    }
  }
  static void frac(real642 &s)
  {
    if (s.val1<0)
      s.val1 -= std::ceil(s.val1);
    else
      s.val1 -= std::floor(s.val1);
    if (s.val2<0)
    {
      double some=std::ceil((double)s.val2);
      while (some)
      {
        s.val2 -= some;
        some=std::ceil((double)s.val2);
      }
    }
    else
    {
      double some=std::floor((double)s.val2);
      while (some)
      {
        s.val2 -= some;
        some=std::floor((double)s.val2);
      }
    }
  }
  static void mod1(real642 &s)
  {
    s.val1 -= std::floor(s.val1);
    double some=std::floor((double)s.val2);
    while (some)
    {
      s.val2 -= some;
      some=std::floor((double)s.val2);
    }
  }
  static void add_double(real642 &s, double x) { s.val1 += x; s.val2 += x; }
  static void mul_double(real642 &s, double x) { s.val1 *= x; s.val2 *= x; }
  static void add(real642 &s, real642 const &o) { s.val1 += o.val1; s.val2 += o.val2; }
  static void sub(real642 &s, real642 const &o) { s.val1 -= o.val1; s.val2 -= o.val2; }
  static void rsub(real642 &s, real642 const &o) { s.val1 = o.val1-s.val1; s.val2 = o.val2-s.val2; }
  static void mul(real642 &s, real642 const &o) { s.val1 *= o.val1; s.val2 *= o.val2; }
  static void sqr(real642 &s) { s.val1 *= s.val1; s.val2 *= s.val2; }
  static double radixfloor(real642 const &s) { return radixfloor2(s.val1); }
  static void recip(real642 &s) { s.val1 = 1.0/s.val1; s.val2 = __float128(1)/s.val2; }
  static void sqrt(real642 &s)
  {
    s.val1 = std::sqrt(s.val1);
    double x=1/std::sqrt((double)s.val2);
    __float128 ax=x*(double)s.val2;
    ax += double(s.val2 - ax*ax)*x/2;
    s.val2 = ax + double(s.val2 - ax*ax)*x/2;
  }
  static bool reduce_angle(real642 &s)
  {
    if (s.val1<-M_PI)
    {
      s.val1+=2*M_PI;
      s.val2+=2.0Q*3.141592653589793238462643383279502884197Q;
      return true;
    }
    else if (s.val1>=M_PI)
    {
      s.val1-=2*M_PI;
      s.val2-=2.0Q*3.141592653589793238462643383279502884197Q;
      return true;
    }
    else
      return false;
  }
  static void add_pi(real642 &s, double x)
  {
    s.val1+=x*M_PI;
    s.val2+=x*3.141592653589793238462643383279502884197Q;
  }
  static std::strong_ordering compare(real642 const &a, real642 const &b)
  {
    if (a.val1 == b.val1)
      return std::strong_ordering::equal;
    else if (a.val1 < b.val1)
      return std::strong_ordering::less;
    else
      return std::strong_ordering::greater;
  }
  static bool isequal(real642 const &a, real642 const &b) { return a.val1==b.val1; }
  static bool is0(real642 const &s) { return s.val1==0; }
  static bool isle(real642 const &a, real642 const &b) { return a.val1<=b.val1; }
  static bool isle0(real642 const &s) { return s.val1<=0; }
  static bool isl0(real642 const &s) { return s.val1<0; }
  static bool isl1(real642 const &s) { return s.val1<1; }
  static void min(real642 &s, real642 const &o) { if (s.val1>o.val1) s=o; }
  static QString toString(real642 const &store) { return float128_toString(store.val2); }
  static bool fromString(real642 &s, char const *str) { return number_traits<double>::fromString(s.val1, str) && float128_fromString(str, s.val2); }
  static int toRound(real642 const &store);
  static double toDouble(real642 const &s) { return s.val1; }
  static real642 *convert_block(NumberType old_type, const void *old_data, int count);
};

template<>
struct number_traits<dd_heap>
{
  using D=number_traits<dd_real>;
  static double eps2() { return D::eps2(); }
  static double eps234() { return D::eps234(); }
  static double sqrtlarge() { return D::sqrtlarge(); }
  static void zero(dd_heap &s, double v) { D::zero(s.value(), v); }
  static void chs(dd_heap &s) { D::chs(s.value()); }
  static void lshift(dd_heap &s, int n) { D::lshift(s.value(), n); }
  static void round(dd_heap &s) { D::round(s.value()); }
  static void frac(dd_heap &s) { D::frac(s.value()); }
  static void mod1(dd_heap &s) { D::mod1(s.value()); }
  static void add_double(dd_heap &s, double x) { D::add_double(s.value(), x); }
  static void mul_double(dd_heap &s, double x) { D::mul_double(s.value(), x); }
  static void add(dd_heap &s, dd_heap const &o) { D::add(s.value(), o.value()); }
  static void sub(dd_heap &s, dd_heap const &o) { D::sub(s.value(), o.value()); }
  static void rsub(dd_heap &s, dd_heap const &o) { D::rsub(s.value(), o.value()); }
  static void mul(dd_heap &s, dd_heap const &o) { D::mul(s.value(), o.value()); }
  static void sqr(dd_heap &s) { D::sqr(s.value()); }
  static double radixfloor(dd_heap const &s) { return D::radixfloor(s.value()); }
  static void recip(dd_heap &s) { D::recip(s.value()); }
  static void sqrt(dd_heap &s) { D::sqrt(s.value()); }
  static bool reduce_angle(dd_heap &s) { return D::reduce_angle(s.value()); }
  static void add_pi(dd_heap &s, double x) { D::add_pi(s.value(), x); }
  static std::strong_ordering compare(dd_heap const &a, dd_heap const &b) { return D::compare(a.value(), b.value()); }
  static bool isequal(dd_heap const &a, dd_heap const &b) { return D::isequal(a.value(), b.value()); }
  static bool is0(dd_heap const &s) { return D::is0(s.value()); }
  static bool isle(dd_heap const &a, dd_heap const &b) { return D::isle(a.value(), b.value()); }
  static bool isle0(dd_heap const &s) { return D::isle0(s.value()); }
  static bool isl0(dd_heap const &s) { return D::isl0(s.value()); }
  static bool isl1(dd_heap const &s) { return D::isl1(s.value()); }
  static void min(dd_heap &s, dd_heap const &o) { D::min(s.value(), o.value()); }
  static QString toString(dd_heap const &s) { return D::toString(s.value()); }
  static bool fromString(dd_heap &s, char const *str) { return D::fromString(s.value(), str); }
  static int toRound(dd_heap const &s) { return D::toRound(s.value()); }
  static double toDouble(dd_heap const &s) { return D::toDouble(s.value()); }
  static dd_real *convert_block(NumberType old_type, const void *old_data, int count) { return D::convert_block(old_type, old_data, count); }
};

// Element type of point arrays and convert_block; differs from BASE only for heap-backed types.
template<typename BASE> struct number_storage { using type=BASE; };
template<> struct number_storage<dd_heap> { using type=dd_real; };

// Value-returning operators and free functions, written in terms of the in-place methods of N.
template<typename N>
class number_operators
{
public:
  friend N &operator+=(N &a, N const &b) { return a.add(b); }
  friend N &operator-=(N &a, N const &b) { return a.sub(b); }
  friend N &operator*=(N &a, N const &b) { return a.mul(b); }
  friend N &operator+=(N &a, double b) { return a.add_double(b); }
  friend N &operator-=(N &a, double b) { return a.add_double(-b); }
  friend N &operator*=(N &a, double b) { return a.mul_double(b); }
  friend N operator+(N a, N const &b) { a.add(b); return a; }
  friend N operator-(N a, N const &b) { a.sub(b); return a; }
  friend N operator*(N a, N const &b) { a.mul(b); return a; }
  friend N operator+(N a, double b) { a.add_double(b); return a; }
  friend N operator-(N a, double b) { a.add_double(-b); return a; }
  friend N operator*(N a, double b) { a.mul_double(b); return a; }
  friend N operator-(N a) { a.chs(); return a; }
  friend bool operator==(N const &a, N const &b) { return a.isequal(b); }
  friend N sqr(N a) { a.sqr(); return a; }
  friend N ldexp(N a, int n) { a.lshift(n); return a; } // TODO: rename to lshift
  friend N recip(N a) { a.recip(); return a; }
  friend N sqrt(N a) { a.sqrt(); return a; }
};

template<typename BASE>
class number: public number_operators<number<BASE>>
{
  public:
  using traits=number_traits<BASE>;
  using storage_t=typename number_storage<BASE>::type;
  protected:
  BASE store;
  template<typename> friend class number;
  template<typename> friend class complex;
  public:
  double eps2() const { return traits::eps2(); }
  double eps234() const { return traits::eps234(); }
  double sqrtlarge() const { return traits::sqrtlarge(); } // limit for squared magnitudes so that one more squaring cannot overflow
  number(): store() { }
  number(NumberType): store() { } // for the same constructor calls as number<number_any>
  void swap(number &other) noexcept
  {
    using std::swap; // use std::swap unless more specific `friend void swap` exists
    swap(store, other.store);
  }
  friend void swap(number &a, number &b) noexcept { a.swap(b); } // used by std::variant::swap in number_any
  NumberType ntype() const { return NumberTypeFromBase<BASE>::ntype; }
  NumberType raw_ntype() const { return NumberTypeFromBase<BASE>::ntype; }
  void readFrom(void *storage, int index) { store=((storage_t *)storage)[index]; } //storage_t storage[index]
  void writeTo(void *storage, int index) const { ((storage_t *)storage)[index]=storage_t(store); } //storage_t storage[index]
  static storage_t *convert_block(NumberType old_type, const void *old_data, int count) { return traits::convert_block(old_type, old_data, count); }

  number &zero(double val=0) { traits::zero(store, val); return *this; }
  number &assign(const number &src) noexcept { store=src.store; return *this; }
  template <typename OTHER_BASE>
  number &assign_across(const number<OTHER_BASE> &src); //1. wrap number<BASE> into number<number_any>, 2. unwrap, 3. convert unrelated number<OTHER>
  number &chs() { traits::chs(store); return *this; }
  number &lshift(int shoft) { traits::lshift(store, shoft); return *this; } // self <<= shoft; 1 lshift -10000 = 0 not error
  number &round() { traits::round(store); return *this; }
  number &frac() { traits::frac(store); return *this; } //-1<result<1
  number &mod1() { traits::mod1(store); return *this; } //0<=result<1
  number &add_double(double x) { traits::add_double(store, x); return *this; }
  number &mul_double(double x) { traits::mul_double(store, x); return *this; }
  number &add(const number &other) { traits::add(store, other.store); return *this; }
  number &sub(const number &other) { traits::sub(store, other.store); return *this; }
  number &rsub(const number &other) { traits::rsub(store, other.store); return *this; }
  number &mul(const number &other) { traits::mul(store, other.store); return *this; }
  number &sqr() { traits::sqr(store); return *this; }
  double radixfloor() const { return traits::radixfloor(store); } //nearest smaller power of 2 (1.5->1->1)
  number &recip() { traits::recip(store); return *this; }
  number &sqrt() { traits::sqrt(store); return *this; }
  bool reduce_angle() { return traits::reduce_angle(store); }
  number &add_pi(double x) { traits::add_pi(store, x); return *this; }
  std::strong_ordering compare(const number &other) const { return traits::compare(store, other.store); }
  std::strong_ordering operator<=>(number const &other) const { return compare(other); }
  bool isequal(const number &other) const { return traits::isequal(store, other.store); }
  bool is0() const { return traits::is0(store); }
  bool isle(const number &other) const { return traits::isle(store, other.store); } //return store<=other
  bool isle0() const { return traits::isle0(store); } //return store<=0
  bool isl0() const { return traits::isl0(store); } //return store<0
  bool isl1() const { return traits::isl1(store); } //return store<1
  number &min(const number &other) { traits::min(store, other.store); return *this; }

  QString toString() const { return traits::toString(store); }
  bool fromString(QString const &str) { return traits::fromString(store, str.toUtf8().constData()); } // decimal; false and unspecified value on error
  int toRound() const { return traits::toRound(store); }
  double toDouble() const { return traits::toDouble(store); }
};

static_assert(std::is_trivially_copyable_v<number<double>>);
static_assert(std::is_trivially_copyable_v<number<__float128>>);
static_assert(std::is_trivially_copyable_v<number<dd_real>>);
static_assert(std::is_trivially_copyable_v<number<dq_real>>);
static_assert(std::is_trivially_copyable_v<number<real642>>);

// Runtime-selected precision; all members are out-of-line in MandelMath_var.cpp.
template<>
class number<number_any>: public number_operators<number<number_any>>
{
  protected:
  number_any store;
  template<typename> friend class number;
  template<typename> friend class complex;
  public:
  double eps2() const;
  double eps234() const;
  double sqrtlarge() const;
  number(): store() { } //for ShareableViewInfo
  number(const number &x) noexcept: store(x.store) { }
  number(number &&x) noexcept: store() { x.swap(*this); };
  number(NumberType ntype); // ntype selects the variant alternative
  ~number();
  void constructLateBecauseQtIsAwesome(const number &src); //for Qt: default-construct, then assign
  void swap(number &other) noexcept { std::swap(store, other.store); }
  NumberType ntype() const;
  NumberType raw_ntype() const;
  void readFrom(void *storage, int index); //BASE storage[index]
  void writeTo(void *storage, int index) const; //BASE storage[index]

  number &zero(double val=0);
  number &assign(const number &src) noexcept;
  number &operator=(number src) noexcept { swap(src); return *this; }
  template <typename OTHER_BASE>
  number &assign_across(const number<OTHER_BASE> &src);
  number &chs();
  number &lshift(int shoft);
  number &round();
  number &frac();
  number &mod1();
  number &add_double(double x);
  number &mul_double(double x);
  number &add(const number &other);
  number &sub(const number &other);
  number &rsub(const number &other);
  number &mul(const number &other);
  number &sqr();
  double radixfloor() const;
  number &recip();
  number &sqrt();
  bool reduce_angle();
  number &add_pi(double x);
  std::strong_ordering compare(const number &other) const;
  std::strong_ordering operator<=>(number const &other) const { return compare(other); }
  bool isequal(const number &other) const;
  bool is0() const;
  bool isle(const number &other) const;
  bool isle0() const;
  bool isl0() const;
  bool isl1() const;
  number &min(const number &other);

  QString toString() const;
  bool fromString(QString const &str);
  int toRound() const;
  double toDouble() const;
};

template <typename BASE>
class complex
{
public:
  number<BASE> re;
  number<BASE> im;
  double eps2() const { return re.eps2(); }
  double eps234() const { return re.eps234(); }
  double sqrtlarge() const { return re.sqrtlarge(); }
  complex(): re(), im() { }
  complex(NumberType ntype): re(ntype), im(ntype) { } // ntype selects variant for number_any; no-op for concrete types
  complex(complex const &src) = default; // no implicit move: a moved-from complex<number_any> must keep its value
  complex &operator=(complex const &src) = default;
  void constructLateBecauseQtIsAwesome(const complex &src);
  //complex(BASE re, BASE im): re(re), im(im) { }
  void readFrom(void *storage, int index); //BASE storage[index], [index+1]
  void writeTo(void *storage, int index) const; //BASE storage[index], [index+1]
  complex &zero(double r=0);
  complex &zero(double r, double i);
  complex &assign(const complex &other) noexcept;
  template <typename OTHER_BASE>
  complex &assign_across(const complex<OTHER_BASE> &src);
  //complex &assign_across_from(const complex<number_any> &src);
  //complex &assign_across_to(complex<number_any> &src) const;
  //void assign_across(const BASE re, const BASE im);
  complex &lshift(int shoft);
  //add_double(r), add_double(r, i)
  complex &add_double(double a) { re.add_double(a); return *this; }
  complex &mul_double(double m);
  complex &add(const complex &other);
  complex &chs();
  complex &sub(const complex &other); //this=this-other
  complex &rsub(const complex &other); //this=other-this
  complex &sqr();
  complex &mul(const number<BASE> &other); //maybe called "scale"
  complex &mul(const complex &other);
  complex &recip();
  complex &recip_prepared(number<BASE> mag); //mag=getMag()
  complex &sqrt();
  complex &pow_uint(unsigned int n); //this^n
  complex &root_approx(int n); //this^(1/n), most real root; only in double precision
  complex &ln_approx(); //this:=ln(this), only in double precision
  complex &exp_approx(); //this:=exp(this), only in double precision
  complex &sign(); //this/=sqrt(this.mag())
  complex &cossin(number<BASE> &angle);
  void arctan2(number<BASE> *result) const;
  double getMag_double() const;
  number<BASE> getMag() const; //re*re+im*im
  number<BASE> getMag1() const; //getMag()-1
  number<BASE> getDist1() const; //|this-1|^2
  number<BASE> mulreT(const complex &other) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  number<BASE> ccw(const complex &other) const; //"counterclockwise" Im(other/this)*|this|^2 = re*o.im-im*o.re
  number<BASE> dist2(const complex &other) const;
  double dist2_double(const complex &other) const;
  complex &from_pmdist(const complex &one, const complex &second); //re:=|o-s|^2, im:=|o+s|^2
  std::strong_ordering compare(complex const &other) const; //compare as vector, im more important
  bool isequal(const complex &other) const;
  bool isle(const complex &other) const { switch(strong_ordering_cast(im<=>other.im)) { case less: return true;
                                                                  case greater: return false;
                                                                  default: return re.isle(other.re); } }
  bool is0() const;
  bool isle0() const; // { return isl0() || is0(); }
  bool isl0() const; //was isNegative() im<0 || im==0 && re<0
  std::strong_ordering mag_cmp_1() const; //(mag() <=> 1) -> -1,0,+1, better return (double)(mag()-1)
  QString toString() const;
  bool fromString(QString const &re_str, QString const &im_str) { return re.fromString(re_str) & im.fromString(im_str); }

  complex &operator+=(complex const &b) { return add(b); }
  complex &operator-=(complex const &b) { return sub(b); }
  complex &operator*=(complex const &b) { return mul(b); }
  complex &operator*=(number<BASE> const &b) { return mul(b); }
  complex &operator+=(double b) { return add_double(b); } //real part only
  complex &operator-=(double b) { return add_double(-b); }
  complex &operator*=(double b) { return mul_double(b); }
  friend complex operator+(complex a, complex const &b) { a.add(b); return a; }
  friend complex operator-(complex a, complex const &b) { a.sub(b); return a; }
  friend complex operator*(complex a, complex const &b) { a.mul(b); return a; }
  friend complex operator*(complex a, number<BASE> const &b) { a.mul(b); return a; }
  friend complex operator+(complex a, double b) { a.add_double(b); return a; }
  friend complex operator-(complex a, double b) { a.add_double(-b); return a; }
  friend complex operator*(complex a, double b) { a.mul_double(b); return a; }
  friend complex operator-(complex a) { a.chs(); return a; }
  friend complex sqr(complex a) { a.sqr(); return a; }
  friend complex ldexp(complex a, int n) { a.lshift(n); return a; }
  friend complex recip(complex a) { a.recip(); return a; }
  friend complex sqrt(complex a) { a.sqrt(); return a; }
};

static_assert(sizeof(number<double>)==8);
static_assert(sizeof(number<__float128>)==16);
static_assert(sizeof(number<dd_real>)==16);
static_assert(sizeof(number<dq_real>)==8*2);
static_assert(sizeof(number<real642>)==8+8+16);
static_assert(sizeof(number_any)==16+32);
static_assert(std::is_trivially_copyable_v<complex<double>>);
static_assert(std::is_trivially_copyable_v<complex<dd_real>>);

template<typename BASE>
void complex<BASE>::readFrom(void *storage, int index)
{
  re.readFrom(storage, index);
  im.readFrom(storage, index+1);
}

template<typename BASE>
void complex<BASE>::writeTo(void *storage, int index) const
{
  re.writeTo(storage, index);
  im.writeTo(storage, index+1);
}

template<typename BASE>
complex<BASE> &complex<BASE>::zero(double r)
{
  re.zero(r);
  im.zero(0);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::zero(double r, double i)
{
  re.zero(r);
  im.zero(i);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::assign(complex const &other) noexcept
{
  re.assign(other.re);
  im.assign(other.im);
  return *this;
}

template<typename BASE>
void complex<BASE>::constructLateBecauseQtIsAwesome(const complex &src)
{
  if constexpr (std::is_same_v<BASE, number_any>)
  {
    re.constructLateBecauseQtIsAwesome(src.re);
    im.constructLateBecauseQtIsAwesome(src.im);
  }
  else
  {
    working_assert(false);
  }
}

template<typename BASE>
complex<BASE> &complex<BASE>::lshift(int shoft)
{
  re.lshift(shoft);
  im.lshift(shoft);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::mul_double(double m)
{
  re.mul_double(m);
  im.mul_double(m);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::add(complex const &other)
{
  re.add(other.re);
  im.add(other.im);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::chs()
{
  re.chs();
  im.chs();
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::sub(complex const &other)
{
  re.sub(other.re);
  im.sub(other.im);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::rsub(complex const &other)
{
  re.rsub(other.re);
  im.rsub(other.im);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::sqr()
{
  //r:=(r+i)*(r-i)
  //i=2*r*i
  //should be more precise than r*r+i*i (+ fewer mul but sqr is cheaper than mul so not really)
  number<BASE> t1(re);
  t1.sub(im);    //r-i
  number<BASE> t2(im);
  im.mul(re);
  im.lshift(1); //2*r*i
  re.add(t2);    //r+i
  re.mul(t1);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::mul(const number<BASE> &other)
{
  re.mul(other);
  im.mul(other);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::mul(complex const &other)
{
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  number<BASE> t1(re);
  number<BASE> t2(im);
  re.mul(other.im);
  im.mul(other.re);
  t1.mul(other.re); //r1*r2
  t2.mul(other.im); //i1*i2
  t1.sub(t2); //r1*r2-i1*i2
  im.add(re); //i1*r2+r1*i2
  re.assign(t1);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::recip()
{
  return recip_prepared(getMag());
}

template<typename BASE>
complex<BASE> &complex<BASE>::recip_prepared(number<BASE> mag)
{ // 1/(re+i*im) = (re-i*im)/((re+i*im)*(re-i*im)) = (re-i*im)/(re*re+im*im)
  mag.recip();
  im.chs();
  re.mul(mag);
  im.mul(mag);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::sqrt()
{
  number<BASE> t1=getMag();
  t1.sqrt();
  if (!re.isle0())
  {
    t1.add(re);
    t1.lshift(-1);
    t1.sqrt(); //t1=sqrt((sqrt(re*re+im*im)+re)/2);
    re.assign(t1); //re=t1
    t1.lshift(1);
    t1.recip();
    im.mul(t1); //im=im/(2*t1);
  }
  else
  {
    t1.sub(re);
    t1.lshift(-1);
    t1.sqrt(); //t1=sqrt((sqrt(re*re+im*im)-re)/2);
    if (t1.isle0()) //t1==0
    {
      re.zero();
      im.zero();
    }
    else
    {
      re.assign(im);
      im.assign(t1); //new im=t1
      t1.lshift(1);
      t1.recip();
      re.mul(t1); //re=old im/(2*t1);
      if (re.isl0())
      {
        re.chs();
        im.chs();
      };
    }
  };
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::pow_uint(unsigned int n)
{
  if (n==1 || is0())
    return *this;
  if (n<=0)
  {
    zero(1, 0);
    return *this;
  }
  if (n==2)
  {
    return sqr();
  };
  number<BASE> t3(re);
  number<BASE> t4(im);

  //works for n>=2
  unsigned int n2=n>>2;
  unsigned int mask=1;
  while (n2)
  {
    mask<<=1;
    n2>>=1;
  }
  //same as mask:=radixfloor(n>>1) ?
  working_assert(mask==1u<<(std::numeric_limits<unsigned int>::digits-std::countl_zero(n>>2)));

  while (mask>0)
  {
    sqr();
    if (n&mask)
    {
      number<BASE> t1(re);
      number<BASE> t2(im);
      re.mul(t4);        //      r1*i2
      im.mul(t3);        //i1*r2
      t1.mul(t3); //r1*r2
      t2.mul(t4); //      i1*i2
      im.add(re);               //i1*r2+r1*i2
      t1.sub(t2); //r1*r2-i1*i2
      re.assign(t1);
    };
    mask>>=1;
  }
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::root_approx(int n)
{
  if (!is0())
  {
    double this_re=re.toDouble();
    double this_im=im.toDouble();
    double r_abs=std::exp(std::log(this_re*this_re+this_im*this_im)/(2*n));
    double r_phi=std::atan2(this_im, this_re)/n;
    double r_re=r_abs*std::cos(r_phi);
    double r_im=r_abs*std::sin(r_phi);
    zero(r_re, r_im);
  };
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::ln_approx()
{
  double r=re.toDouble();
  double i=im.toDouble();
  return zero(std::log(r*r+i*i)/2, std::atan2(i, r));
}

template<typename BASE>
complex<BASE> &complex<BASE>::exp_approx()
{
  double d=std::exp(re.toDouble());
  double p=im.toDouble();
  return zero(d*std::cos(p), d*std::sin(p));
}

template<typename BASE>
complex<BASE> &complex<BASE>::sign()
{
  number<BASE> t1=getMag();
  t1.sqrt();
  t1.recip();
  re.mul(t1);
  im.mul(t1);
  return *this;
}

template<typename BASE>
complex<BASE> &complex<BASE>::cossin(number<BASE> &angle)
{
  return zero(std::cos(angle.toDouble()), std::sin(angle.toDouble()));
}

template<typename BASE>
void complex<BASE>::arctan2(number<BASE> *result) const
{
  if (im.is0() && re.isl0())
    result->zero(std::atan2(-0, 1));
  else
    result->zero(std::atan2(im.toDouble(), re.toDouble()));
}

template<typename BASE>
double complex<BASE>::getMag_double() const
{
  return sqr_double(re.toDouble())+sqr_double(im.toDouble());
}

template<typename BASE>
number<BASE> complex<BASE>::getMag() const
{
  number<BASE> t1(re);
  t1.sqr();
  number<BASE> t2(im);
  t2.sqr();
  t1.add(t2);
  return t1;
}

template<typename BASE>
number<BASE> complex<BASE>::getMag1() const
{
  number<BASE> t1=getMag();
  t1.add_double(-1);
  return t1;
}

template<typename BASE>
number<BASE> complex<BASE>::getDist1() const
{
  number<BASE> t1(re);
  t1.add_double(-1);
  t1.sqr();
  number<BASE> t2(im);
  t2.sqr();
  t1.add(t2);
  return t1;
}

template<typename BASE>
number<BASE> complex<BASE>::mulreT(complex const &other) const
{
  number<BASE> t1(re);
  number<BASE> t2(im);
  t1.mul(other.re);
  t2.mul(other.im);
  t1.add(t2);
  return t1;
}

template<typename BASE>
number<BASE> complex<BASE>::ccw(complex const &other) const
{
  number<BASE> t1(re);
  number<BASE> t2(im);
  t1.mul(other.im);
  t2.mul(other.re);
  t1.sub(t2); //re*o.im-im*o.re
  return t1;
}

template<typename BASE>
double complex<BASE>::dist2_double(complex const &other) const
{
  number<BASE> t1(re);
  number<BASE> t2(im);
  t1.sub(other.re);
  t2.sub(other.im);
  return sqr_double(t1.toDouble())+sqr_double(t2.toDouble());
}

template<typename BASE>
number<BASE> complex<BASE>::dist2(complex const &other) const
{
  number<BASE> t1(re);
  number<BASE> t2(im);
  t1.sub(other.re);
  t2.sub(other.im);
  t1.sqr();
  t2.sqr();
  t1.add(t2);
  return t1;
}

template<typename BASE>
complex<BASE> &complex<BASE>::from_pmdist(const complex &one, const complex &second)
{
  //re:=(a-c)^2+(b-d)^2=a^2-2*a*c+c^2+b^2-2*b*d+d^2=a^2+b^2+c^2+d^2 -2*a*c-2*b*d
  //im:=(a+c)^2+(b+d)^2=a^2+2*a*c+c^2+b^2+2*b*d+d^2=a^2+b^2+c^2+d^2 +2*a*c+2*b*d
  number<BASE> t1(one.re);
  im.assign(one.im);
  t1.mul(second.re);
  im.mul(second.im);
  t1.add(im); //a*c+b*d
  t1.lshift(1);

  re.assign(one.re);
  re.sqr();
  im.assign(one.im);
  im.sqr();
  re.add(im);
  im.assign(second.re);
  im.sqr();
  re.add(im);
  im.assign(second.im);
  im.sqr();
  re.add(im); //a^2+b^2+c^2+d^2

  im.assign(re);
  re.sub(t1);
  im.add(t1);
  return *this;
}

template<typename BASE>
bool complex<BASE>::isequal(complex const &other) const
{
  return re.isequal(other.re) &&
         im.isequal(other.im);
}

template<typename BASE>
bool complex<BASE>::is0() const
{
  return re.is0() && im.is0();
}

template<typename BASE>
bool complex<BASE>::isle0() const
{
  //im<0 || im==0 && re<=0
  return im.isl0() || (im.is0() && re.isle0());
}

template<typename BASE>
bool complex<BASE>::isl0() const
{
  //im<0 || im==0 && re<0
  return im.isl0() || (im.is0() && re.isl0());
}

template<typename BASE>
std::strong_ordering complex<BASE>::mag_cmp_1() const
{
  //first reduce fz_re, fz_im to 2nd octant:
  //  fz.re>=0, fz.im>=0, fz.im>=fz.re
  number<BASE> t1(re);
  number<BASE> t2(im);
  if (t1.isl0()) t1.chs();
  if (t2.isl0()) t2.chs();
  if (!t1.isle(t2))
  {
    number<BASE> t3(t1);
    t1.assign(t2);
    t2.assign(t3);
  };
  //im>=re, now check re^2+im^2>1
  if ((t1.toDouble()>=0.71) || //im>=re>=0.71 -> mag>=1.0082
      (t2.toDouble()>=1.01))   //im>=1.01 -> mag>=1.01
  {
    return std::strong_ordering::greater;
  }
  else if (t2.toDouble()<=0.70) //0.70>=im>=re -> mag<=0.98
    return std::strong_ordering::less;
  else
  { //   re*re+im*im>1
    //   gets inaccurate for re~0, im~1
    //   im*im>1-re*re ; im>sqrt(1-re*re) ; im-1>sqrt(1-re*re)-1
    //>> x=Sqrt(1-r*r)-1
    //   one cycle of Newton should work with any precision
    //   (x+1)*(x+1)=(1-r*r)
    //   x*x+2*x+r*r=0  f'=2*x+2
    //   x2=x-(x*x+2*x+r*r)/(2*x+2)
    //>> x2=(x*x-r*r)/(2*x+2)
    t1.sqr(); //r^2
    t1.chs();
    number<BASE> t3(t1);
    t3.add_double(1);
    t3.sqrt();              //sqrt(1-r*r) = x+1
    number<BASE> t4(t3);
    t4.add_double(-1);        //x
    t4.sqr();
    t1.add(t4);    //x*x-r*r
    t3.recip(); //-1..-0.51 -> -1..-1.96
    t1.mul(t3);
    t1.lshift(-1);           //(x*x-r*r)/(x+1)/2 = better x
    t2.add_double(-1);
    return t2.compare(t1); //im-1 vs limit
  };
}

void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im); //res_re>=0
double radixfloor_double(double x1, double x2);
void real_double_quadratic(double *res, double a, double b2, double c);
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);

} // namespace MandelMath
#endif // MANDELMATH_NUMBER_HPP
