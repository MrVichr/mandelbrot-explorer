#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
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
#define working_assert(x) { if (!(x)) dbgPoint(); }
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

template<int N, typename... Ts> using NthTypeOf =
    typename std::tuple_element<N, std::tuple<Ts...>>::type;
template<typename ...Options>
auto mvisit(auto callable, std::variant<Options...> &vari, auto &...params)
{
  //using cases=std::make_index_sequence<sizeof...(Options)>;
  auto index=vari.index(); //size_t
  if constexpr (sizeof...(Options)>0)
      if (index==0) return callable(std::get<0>(vari), params...); //__detail::__variant::__get<_Np>(__v);
  if constexpr (sizeof...(Options)>1)
      //if (index==1) return callable(std::get<1>(vari), params...);
      if (index==1) return callable((NthTypeOf<1, Options...> &)vari, params...);
  if constexpr (sizeof...(Options)>2)
      if (index==2) return callable(std::get<2>(vari), params...);
  if constexpr (sizeof...(Options)>3)
      if (index==3) return callable(std::get<3>(vari), params...);
  if constexpr (sizeof...(Options)>4)
      if (index==4) return callable(std::get<4>(vari), params...);
  if constexpr (sizeof...(Options)>5)
      if (index==5) return callable(std::get<5>(vari), params...);
  if constexpr (sizeof...(Options)>6)
      if (index==6) return callable(std::get<6>(vari), params...);
  using RT=decltype(callable((NthTypeOf<0, Options...> &)vari, params...));
  [[unlikely]]
  if constexpr (std::is_same_v<RT, std::strong_ordering>)
    return std::strong_ordering::equal; //even more failure from C++ commitee
  else
    return RT();
}

template<typename ...Options>
auto mvisit(auto callable, std::variant<Options...> const &vari, auto &...params)
{
  //using cases=std::make_index_sequence<sizeof...(Options)>;
  auto index=vari.index(); //size_t
  if constexpr (sizeof...(Options)>0)
      if (index==0) return callable(std::get<0>(vari), params...); //__detail::__variant::__get<_Np>(__v);
  if constexpr (sizeof...(Options)>1)
      //if (index==1) return callable(std::get<1>(vari), params...);
      if (index==1) return callable((NthTypeOf<1, Options...> const &)vari, params...);
  if constexpr (sizeof...(Options)>2)
      if (index==2) return callable(std::get<2>(vari), params...);
  if constexpr (sizeof...(Options)>3)
      if (index==3) return callable(std::get<3>(vari), params...);
  if constexpr (sizeof...(Options)>4)
      if (index==4) return callable(std::get<4>(vari), params...);
  if constexpr (sizeof...(Options)>5)
      if (index==5) return callable(std::get<5>(vari), params...);
  if constexpr (sizeof...(Options)>6)
      if (index==6) return callable(std::get<6>(vari), params...);
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
//template <> struct NumberTypeFromBase<dq_real> { static constexpr NumberType ntype=typeQDouble; };
template <> struct NumberTypeFromBase<real642> { static constexpr NumberType ntype=typeReal642; };
#endif

extern int CreatedInstancesOfNumber;

template<typename BASE>
class number;

using number_any=std::variant<std::monostate, number<double>, number<__float128>, number<dd_real>, number<dq_real>, number<real642>>;

template<typename T>
concept IMandelNumber_ = requires(T a, T b)
{
  //no worky using namespace std;
  //simple way:
  //{ a.add(b) } -> std::same_as<T>;
  //strict way:
  //requires std::same_as<decltype ( &T::add ), T &(T::*)(T const &)>;

  { a.eps2() } -> std::same_as<double>;
  requires std::same_as<decltype ( &T::eps234 ), double(T::*)() const>;

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
  requires IMandelNumber_<T>;
  requires std::same_as<decltype ( &T::round ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::frac ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::mod1 ), T &(T::*)()>;
  requires std::same_as<decltype ( &T::radixfloor ), double (T::*)(T const &other) const>;
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
  requires IMandelNumber_<T>;
};

#if 0
class number_a
{
protected:
  /*friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<real642>;
  friend class number<number<bool> *>;*/
public:
  virtual double eps2() const=0;
  virtual double eps234() const=0;
  number_a();//: store() {}
  //number(NumberType ntype);
  //number(BASE store): store(store) {}
  virtual ~number_a();
  virtual NumberType ntype() const=0;// { return this->_ntype; }
  virtual NumberType raw_ntype() const=0;// { return this->_ntype; }
  virtual void readFrom(void *storage, int index)=0; //BASE storage[index]
  virtual void writeTo(void *storage, int index) const=0; //BASE storage[index]
  //static void *convert_block(NumberType old_type, const void *old_data, int count);

  virtual void zero(double val=0)=0;
  virtual void assign(const number_a &src)=0;
  virtual void assign_across(const number_a *src)=0;
  virtual void chs()=0;
  virtual void lshift(int shoft)=0; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual void round()=0;
  virtual void frac()=0; //-1<result<1
  virtual void mod1()=0; //0<=result<1
  virtual void add_double(double x)=0;
  virtual void mul_double(double x)=0;
  virtual void add(const number_a &other)=0;
  virtual void sub(const number_a &other)=0;
  virtual void rsub(const number_a &other)=0;
  virtual void mul(const number_a &other)=0;
  virtual void sqr()=0;
  virtual double radixfloor(const number_a &other) const=0; //nearest smaller power of 2 (1.5->1->1)
  virtual void recip()=0;
  virtual void sqrt()=0;
  virtual bool reduce_angle()=0; //one step towards -pi<=x<pi, true if changed
  virtual void add_pi(double x)=0;
  virtual int compare(const number_a &other) const=0; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  virtual bool isequal(const number_a &other) const=0; //return store==other
  virtual bool is0() const=0;
  virtual bool isle(const number_a &other) const=0; //return store<=other
  virtual bool isle0() const=0; //return store<=0
  virtual bool isl0() const=0; //return store<0
  virtual bool isl1() const=0; //return store<1
  virtual void min(const number_a &other)=0;

  virtual QString toString() const=0;
  virtual int toRound() const=0;
  virtual double toDouble() const=0;
};
#endif

template<typename BASE>
class number
{
  protected:
  BASE store; //either "double" for special optimized template
      //or "double *" for universal implementation
      //going meta by using BASE="number<actual base> *"
  friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<dq_real>;
  friend class number<real642>;
  //friend class number<number<bool> *>;
  public:
  double eps2() const;
  double eps234() const;
  number(): store() { } //for ShareableViewInfo
  number(const number &x) noexcept: store(x.store) { };
  number(number &&x) noexcept: store() { x.swap(*this); }
  number(NumberType ntype);
  void constructLateBecauseQtIsAwesome(NumberType ntype);
  void swap(number &other) noexcept { std::swap(store, other.store); }
  //number(BASE store): store(store) {}
  ~number();
  NumberType ntype() const;// { return this->_ntype; }
  NumberType raw_ntype() const;
  void readFrom(void *storage, int index); //BASE storage[index]
  void writeTo(void *storage, int index) const; //BASE storage[index]
  static BASE *convert_block(NumberType old_type, const void *old_data, int count);

  number &zero(double val=0);
  number &assign(const number &src) noexcept;// { store=src.store; return *this; } //or { BASE tmp(src.store); std::swap(tmp, store); return *this; }
  //number &operator=(number &src) noexcept { swap(src); return *this; } //required by std::variant<>
  number &operator=(number src) noexcept { swap(src); return *this; }
  //void assign(const MandelNumber &src);
  //number &assign_across(number_a const &src);
  template<typename OTHER>
  number<BASE> &assign_across(number<OTHER> const &src);
  number &chs();
  number &lshift(int shoft); // self <<= shoft; 1 lshift -10000 = 0 not error
  number &round();
  number &frac(); //-1<result<1
  number &mod1(); //0<=result<1
  number &add_double(double x);
  number &mul_double(double x);
  //void add(const number_a &other);
  number &add(number const &other);
  //void sub(const number_a &other);
  number &sub(number const &other);
  //virtual void rsub(const number_a &other) override;
  number &rsub(number const &other);
  //virtual void mul(const number_a &other) override;
  number &mul(number const &other);
  number &sqr();
  //virtual double radixfloor(const number_a &other) const override;
  double radixfloor(number const &other) const; //nearest smaller power of 2 (1.5->1->1)
  number &recip();
  number &sqrt();
  bool reduce_angle();
  number &add_pi(double x);
  //virtual int compare(const number_a &other) const override; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering compare(number const &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering operator<=>(number const &other) const { return compare(other); }
  //virtual bool isequal(const number_a &other) const override; //return store==other
  bool isequal(number const &other) const; //return store==other
  bool is0() const;
  //virtual bool isle(const number_a &other) const override; //return store<=other
  bool isle(number const &other) const; //return store<=other
  bool isle0() const; //return store<=0
  bool isl0() const; //return store<0
  bool isl1() const; //return store<1
  //virtual void min(const number_a &other) override;
  number &min(number const &other);

  QString toString() const;
  int toRound() const;
  double toDouble() const;

  //static_assert(std::is_nothrow_assignable_v<number<BASE>, number<BASE>>, "why how");
  //static_assert(std::is_nothrow_copy_assignable_v<number<BASE>>, "why how");
  //static_assert(std::is_nothrow_move_constructible_v<number<BASE>>, "why how");
  //static_assert(std::is_nothrow_copy_constructible_v<number<BASE>>, "why how");
};

static_assert(IMandelReal<number<double>>, "fix that");
static_assert(sizeof(number<double>)==sizeof(double), "fix this");
/*static_assert(std::is_nothrow_assignable_v<number<double>, number<double>>, "why how");
static_assert(std::is_nothrow_copy_assignable_v<number<double>>, "why how");
static_assert(std::is_nothrow_move_constructible_v<number<double>>, "why how");
static_assert(std::is_nothrow_copy_constructible_v<number<double>>, "why how");*/

//static_assert(Number_double requires IMandelNumber, "fix that");
static_assert(IMandelReal<number<__float128>>, "fix that");
static_assert(sizeof(number<__float128>)==sizeof(__float128), "fix this");

static_assert(IMandelReal<number<dd_real>>, "fix that");
static_assert(IMandelReal<number<real642>>, "fix that");
static_assert(IMandelReal<number<number_any>>, "fix that");

#if 0
template <IMandelNumber BASE>
class number: public number_a
{
protected:
  BASE store; //either "double" for special optimized template
              //or "double *" for universal implementation
              //going meta by using BASE="number<actual base> *"
  friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<real642>;
  friend class number<number_a *>;
public:
  virtual double eps2() const override;
  virtual double eps234() const override;
  number(): number_a(), store() { } //for ShareableViewInfo
  number(number &x)=delete;
  number(number &&x)=delete;
  number(NumberType ntype);
  void constructLateBecauseQtIsAwesome(NumberType ntype);
  //number(BASE store): store(store) {}
  virtual ~number() override;
  virtual NumberType ntype() const override;// { return this->_ntype; }
  virtual NumberType raw_ntype() const override;
  virtual void readFrom(void *storage, int index) override; //BASE storage[index]
  virtual void writeTo(void *storage, int index) const override; //BASE storage[index]
  static BASE *convert_block(NumberType old_type, const void *old_data, int count);
  //static worker_multi *create(Type ntype, int capacity);
  //static worker_multi *create(Type ntype, Allocator<worker_multi> *allocator);
  //Allocator<worker_multi> *getAllocator() { return &allocator; }
  //virtual void assign_block(int dst, worker_multi *src_worker, int src, int len)=0; //for now, assert(this.ntype==src.ntype)
    //or assign_block(int dst, allocator *src) + assign_block(allocator *dst, int src) + assign_block(int dst, int src, int len)
  //virtual void assign_block(int dst, const Allocator<worker_multi> *src)
   //{ int first, capac; src->_getFirstCapac(first, capac); assign_block(dst, src->worker, first, capac); }

  virtual void zero(double val=0) override;
  virtual void assign(const number_a &src) override;
  void assign(const number<BASE> &src);
  virtual void assign_across(const number_a *src) override;
  //template <typename OTHER_BASE>
  //void assign_across(const number<OTHER_BASE> *src);
  void assign_across(const number<BASE> *src);
  //void assign_across(BASE src);
  virtual void chs() override;
  virtual void lshift(int shoft) override; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual void round() override;
  virtual void frac() override; //-1<result<1
  virtual void mod1() override; //0<=result<1
  virtual void add_double(double x) override;
  virtual void mul_double(double x) override;
  virtual void add(const number_a &other) override;
  void add(const number<BASE> &other);
  virtual void sub(const number_a &other) override;
  void sub(const number<BASE> &other);
  virtual void rsub(const number_a &other) override;
  void rsub(const number<BASE> &other);
  virtual void mul(const number_a &other) override;
  void mul(const number<BASE> &other);
  virtual void sqr() override;
  virtual double radixfloor(const number_a &other) const override;
  double radixfloor(const number<BASE> &other) const; //nearest smaller power of 2 (1.5->1->1)
  virtual void recip() override;
  virtual void sqrt() override;
  virtual bool reduce_angle() override;
  virtual void add_pi(double x) override;
  virtual int compare(const number_a &other) const override; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  int compare(const number<BASE> &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  virtual bool isequal(const number_a &other) const override; //return store==other
  bool isequal(const number<BASE> &other) const; //return store==other
  virtual bool is0() const override;
  virtual bool isle(const number_a &other) const override; //return store<=other
  bool isle(const number<BASE> &other) const; //return store<=other
  virtual bool isle0() const override; //return store<=0
  virtual bool isl0() const override; //return store<0
  virtual bool isl1() const override; //return store<1
  virtual void min(const number_a &other) override;
  void min(const number<BASE> &other);

  virtual QString toString() const override;
  virtual int toRound() const override;
  virtual double toDouble() const override;
};
#endif

template <typename BASE>
class complex
{
  public:
  struct Scratchpad
  {
      NumberType ntype;
      number<BASE> tmp1;
      number<BASE> tmp2;
      number<BASE> tmp3;
      number<BASE> tmp4;
      Scratchpad(): ntype(NumberType::typeEmpty), tmp1(), tmp2(), tmp3(), tmp4() { }
      Scratchpad(NumberType ntype): ntype(ntype), tmp1(ntype), tmp2(ntype), tmp3(ntype), tmp4(ntype) { }
  };
protected:
  Scratchpad *tmp;
public:
  number<BASE> re;
  number<BASE> im;
  double eps2() const { return re.eps2(); }
  double eps234() const { return re.eps234(); }
  complex(): re(), im() { }
  //complex(NumberType ntype): re(ntype), im(ntype) { }
  complex(Scratchpad *spad): tmp(spad), re(spad->ntype), im(spad->ntype) { }
  complex(complex const &src) noexcept: tmp(src.tmp), re(src.re), im(src.im) { }
  //complex(BASE re, BASE im): re(re), im(im) { }
  void readFrom(void *storage, int index); //BASE storage[index], [index+1]
  void writeTo(void *storage, int index) const; //BASE storage[index], [index+1]
  complex &zero(double r=0);
  complex &zero(double r, double i);
  complex &assign(const complex &other) noexcept;
  template <typename OTHER_BASE>
  void assign_across(const complex<OTHER_BASE> &src);
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
  complex &recip_prepared();
  complex &sqrt();
  complex &pow_uint(unsigned int n); //this^n
  complex &root_approx(int n); //this^(1/n), most real root; only in double precision
  complex &ln_approx(); //this:=ln(this), only in double precision
  complex &exp_approx(); //this:=exp(this), only in double precision
  complex &sign(); //this/=sqrt(this.mag())
  complex &cossin(number<BASE> &angle);
  void arctan2(number<BASE> *result) const;
  double getMag_double() const;
  const number<BASE> *getMag_tmp() const;
  const number<BASE> *getMag1_tmp() const;
  const number<BASE> *getDist1_tmp() const;
  const number<BASE> *mulreT_tmp(const complex &other) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  const number<BASE> *ccw_tmp(const complex &other) const; //"counterclockwise" Im(other/this)*|this|^2 = re*o.im-im*o.re
  double dist2_double(const complex &other) const;
  const number<BASE> *dist2_tmp(const complex &other) const;
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
};

static_assert(IMandelComplex<complex<double>>, "oh well");

double sqr_double(double x); //no one ever needed this function before year 2022, right?
void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im); //res_re>=0
double radixfloor_double(double x1, double x2);
void real_double_quadratic(double *res, double a, double b2, double c);
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im);

} // namespace MandelMath

namespace std { namespace __detail { namespace __variant {
/*template<> struct _Never_valueless_alt<MandelMath::number<double>>: std::bool_constant<true> { };
template<> struct _Never_valueless_alt<MandelMath::number<__float128>>: std::bool_constant<true> { };
template<> struct _Never_valueless_alt<MandelMath::number<MandelMath::dd_real>>: std::bool_constant<true> { };
template<> struct _Never_valueless_alt<MandelMath::number<MandelMath::dq_real>>: std::bool_constant<true> { };
template<> struct _Never_valueless_alt<MandelMath::number<MandelMath::real642>>: std::bool_constant<true> { };*/
template<> constexpr bool __never_valueless<std::monostate, MandelMath::number<double>, MandelMath::number<__float128>, MandelMath::number<MandelMath::dd_real>, MandelMath::number<MandelMath::dq_real>, MandelMath::number<MandelMath::real642>>() { return true; }
} } }
#endif // MANDELMATH_NUMBER_HPP
