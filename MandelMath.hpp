#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
//#include "multiprec.hpp"

#define NUMBER_DOUBLE_EXISTS 1
#define NUMBER_DOUBLE_ONLY 0
#define NUMBER_IS_VARIANT 1

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

void atomic_min(std::atomic<int> &a, int val) noexcept;
void atomic_max(std::atomic<int> &a, int val) noexcept;

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
template <> struct NumberTypeFromBase<dq_real> { static constexpr NumberType ntype=typeQDouble; };
template <> struct NumberTypeFromBase<real642> { static constexpr NumberType ntype=typeReal642; };
#endif

extern int CreatedInstancesOfNumber;

template<typename BASE>
class number;

#if NUMBER_IS_VARIANT
using number_any=std::variant<std::monostate, number<double>, number<__float128>, number<dd_real>, number<dq_real>, number<real642>>;
#endif

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

#if NUMBER_IS_VARIANT
template<typename BASE>
class number
{
  public:
  struct Scratchpad
  {
    NumberType ntype;
    void *inner;
    number<BASE> tmp1;
    number<BASE> tmp2;
    number<BASE> tmp3;
    number<BASE> tmp4;
    //Scratchpad(): ntype(NumberType::typeEmpty), inner(nullptr), tmp1(), tmp2(), tmp3(), tmp4() { }
    Scratchpad(NumberType ntype);//: ntype(ntype), tmp1(this), tmp2(this), tmp3(this), tmp4(this), inner(nullptr) { }
    //Scratchpad(Scratchpad *inner): ntype(inner->ntype), tmp1(this), tmp2(this), tmp3(this), tmp4(this), inner(inner) { }
    ~Scratchpad();
  };
  protected:
  BASE store; //either "double" for special optimized template
      //or "double *" for universal implementation
      //going meta by using BASE="number<actual base> *"
  Scratchpad *tmp;
  //friend class ::ShareableViewInfo;
  friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<dq_real>;
  friend class number<real642>;
  friend class number<number_any>;
  public:
  double eps2() const;
  double eps234() const;
  number(): store(), tmp(nullptr) { } //for ShareableViewInfo
  number(const number &x) noexcept: store(x.store), tmp(x.tmp) { }
  number(number &&x) noexcept: store(), tmp(x.tmp) { x.swap(*this); };
  number(Scratchpad *spad);
  ~number();
  void constructLateBecauseQtIsAwesome(const number &arc); //starts to look like copy_construct
  void swap(number &other) noexcept { std::swap(store, other.store); }
  //number(BASE store): store(store) {}
  NumberType ntype() const;// { return this->_ntype; }
  NumberType raw_ntype() const;
  void readFrom(void *storage, int index); //BASE storage[index]
  void writeTo(void *storage, int index) const; //BASE storage[index]
  static BASE *convert_block(NumberType old_type, const void *old_data, int count); //TODO: probably needs to introduce number<OTHER>::serialized, and accept variant<serialized1 *, serialized2 *,...>
  //static worker_multi *create(Type ntype, int capacity);
  //static worker_multi *create(Type ntype, Allocator<worker_multi> *allocator);
  //Allocator<worker_multi> *getAllocator() { return &allocator; }
  //virtual void assign_block(int dst, worker_multi *src_worker, int src, int len)=0; //for now, assert(this.ntype==src.ntype)
  //or assign_block(int dst, allocator *src) + assign_block(allocator *dst, int src) + assign_block(int dst, int src, int len)
  //virtual void assign_block(int dst, const Allocator<worker_multi> *src)
  //{ int first, capac; src->_getFirstCapac(first, capac); assign_block(dst, src->worker, first, capac); }

  number &zero(double val=0);
  //number &assign(const number_any &src);
  number &assign(const number<BASE> &src) noexcept;
  number &operator=(number<BASE> src) noexcept { swap(src); return *this; }
  template <typename OTHER_BASE>
  number &assign_across(const number<OTHER_BASE> &src); //1. wrap number<BASE> into number<number_any>, 2. unwrap, 3. convert unrelated number<OTHER>
  //template <typename OTHER_BASE>
  //void assign_across(const number<OTHER_BASE> *src);
  //number &assign_across(const number<BASE> &src);
  //void assign_across(BASE src);
  number &chs();
  number &lshift(int shoft); // self <<= shoft; 1 lshift -10000 = 0 not error
  number &round();
  number &frac(); //-1<result<1
  number &mod1(); //0<=result<1
  number &add_double(double x);
  number &mul_double(double x);
  number &add(const number_any &other);
  number &add(const number<BASE> &other);
  number &sub(const number_any &other);
  number &sub(const number<BASE> &other);
  number &rsub(const number_any &other);
  number &rsub(const number<BASE> &other);
  number &mul(const number_any &other);
  number &mul(const number<BASE> &other);
  number &sqr();
  double radixfloor() const;
  //double radixfloor(const number<BASE> &other) const; //nearest smaller power of 2 (1.5->1->1)
  number &recip();
  number &sqrt();
  bool reduce_angle();
  number &add_pi(double x);
  //std::strong_ordering compare(const number_any &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering compare(const number<BASE> &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  //std::strong_ordering operator<=>(number_any const &other) const { return compare(other); }
  std::strong_ordering operator<=>(number<BASE> const &other) const { return compare(other); }
  bool isequal(const number_any &other) const; //return store==other
  bool isequal(const number<BASE> &other) const; //return store==other
  bool is0() const;
  bool isle(const number_any &other) const; //return store<=other
  bool isle(const number<BASE> &other) const; //return store<=other
  bool isle0() const; //return store<=0
  bool isl0() const; //return store<0
  bool isl1() const; //return store<1
  number &min(const number_any &other);
  number &min(const number<BASE> &other);

  QString toString() const;
  int toRound() const;
  double toDouble() const;
};
#else
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

  virtual number_a &zero(double val=0)=0;
  virtual number_a &assign(const number_a &src)=0;
  virtual number_a &assign_across(const number_a &src)=0;
  virtual number_a &chs()=0;
  virtual number_a &lshift(int shoft)=0; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual number_a &round()=0;
  virtual number_a &frac()=0; //-1<result<1
  virtual number_a &mod1()=0; //0<=result<1
  virtual number_a &add_double(double x)=0;
  virtual number_a &mul_double(double x)=0;
  virtual number_a &add(const number_a &other)=0;
  virtual number_a &sub(const number_a &other)=0;
  virtual number_a &rsub(const number_a &other)=0;
  virtual number_a &mul(const number_a &other)=0;
  virtual number_a &sqr()=0;
  virtual double radixfloor() const=0; //nearest smaller power of 2 (1.5->1->1)
  virtual number_a &recip()=0;
  virtual number_a &sqrt()=0;
  virtual bool reduce_angle()=0; //one step towards -pi<=x<pi, true if changed
  virtual number_a &add_pi(double x)=0;
  virtual std::strong_ordering compare(const number_a &other) const=0; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering operator<=>(number_a const &other) const { return compare(other); }
  virtual bool isequal(const number_a &other) const=0; //return store==other
  virtual bool is0() const=0;
  virtual bool isle(const number_a &other) const=0; //return store<=other
  virtual bool isle0() const=0; //return store<=0
  virtual bool isl0() const=0; //return store<0
  virtual bool isl1() const=0; //return store<1
  virtual number_a &min(const number_a &other)=0;

  virtual QString toString() const=0;
  virtual int toRound() const=0;
  virtual double toDouble() const=0;
};

using number_any=number_a *;

template<typename BASE>
class number: public number_a
{
  public:
  struct Scratchpad
  {
      NumberType ntype;
      void *inner;
      number<BASE> tmp1;
      number<BASE> tmp2;
      number<BASE> tmp3;
      number<BASE> tmp4;
      //Scratchpad(): ntype(NumberType::typeEmpty), inner(nullptr), tmp1(), tmp2(), tmp3(), tmp4() { }
      Scratchpad(NumberType ntype);//: ntype(ntype), tmp1(this), tmp2(this), tmp3(this), tmp4(this), inner(nullptr) { }
      //Scratchpad(Scratchpad *inner): ntype(inner->ntype), tmp1(this), tmp2(this), tmp3(this), tmp4(this), inner(inner) { }
      ~Scratchpad();
  };
protected:
  BASE store; //either "double" for special optimized template
              //or "double *" for universal implementation
              //going meta by using BASE="number<actual base> *"
  Scratchpad *tmp;
  //friend class ::ShareableViewInfo;
  friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<dq_real>;
  friend class number<real642>;
  friend class number<number_any>;
public:
  virtual double eps2() const override;
  virtual double eps234() const override;
  number(): number_a(), store(), tmp(nullptr) { } //for ShareableViewInfo
  number(const number &x) noexcept;
  number(number &&x) noexcept: store(), tmp(x.tmp) { x.swap(*this); };
  number(Scratchpad *spad);
  void constructLateBecauseQtIsAwesome(const number &arc); //starts to look like copy_construct
  void swap(number &other) noexcept { std::swap(store, other.store); }
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

  virtual number &zero(double val=0) override;
  virtual number &assign(const number_a &src) override;
  number &assign(const number<BASE> &src) noexcept;
  virtual number &assign_across(const number_a &src) override;
  //template <typename OTHER_BASE>
  //void assign_across(const number<OTHER_BASE> *src);
  number &assign_across(const number<BASE> &src);
  //void assign_across(BASE src);
  virtual number &chs() override;
  virtual number &lshift(int shoft) override; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual number &round() override;
  virtual number &frac() override; //-1<result<1
  virtual number &mod1() override; //0<=result<1
  virtual number &add_double(double x) override;
  virtual number &mul_double(double x) override;
  virtual number &add(const number_a &other) override;
  number &add(const number<BASE> &other);
  virtual number &sub(const number_a &other) override;
  number &sub(const number<BASE> &other);
  virtual number &rsub(const number_a &other) override;
  number &rsub(const number<BASE> &other);
  virtual number &mul(const number_a &other) override;
  number &mul(const number<BASE> &other);
  virtual number &sqr() override;
  virtual double radixfloor() const override;
  //double radixfloor(const number<BASE> &other) const; //nearest smaller power of 2 (1.5->1->1)
  virtual number &recip() override;
  virtual number &sqrt() override;
  virtual bool reduce_angle() override;
  virtual number &add_pi(double x) override;
  virtual std::strong_ordering compare(const number_a &other) const override; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering compare(const number<BASE> &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  std::strong_ordering operator<=>(number_a const &other) const { return compare(other); }
  std::strong_ordering operator<=>(number<BASE> const &other) const { return compare(other); }
  virtual bool isequal(const number_a &other) const override; //return store==other
  bool isequal(const number<BASE> &other) const; //return store==other
  virtual bool is0() const override;
  virtual bool isle(const number_a &other) const override; //return store<=other
  bool isle(const number<BASE> &other) const; //return store<=other
  virtual bool isle0() const override; //return store<=0
  virtual bool isl0() const override; //return store<0
  virtual bool isl1() const override; //return store<1
  virtual number &min(const number_a &other) override;
  number &min(const number<BASE> &other);

  virtual QString toString() const override;
  virtual int toRound() const override;
  virtual double toDouble() const override;
};
#endif

template <typename BASE>
class complex
{
protected:
  using Scratchpad=typename number<BASE>::Scratchpad;
  Scratchpad *tmp;
public:
  number<BASE> re;
  number<BASE> im;
  double eps2() const { return re.eps2(); }
  double eps234() const { return re.eps234(); }
  complex(): tmp(nullptr), re(), im() { }
  //complex(NumberType ntype): re(ntype), im(ntype) { }
  complex(Scratchpad *spad): tmp(spad), re(spad), im(spad) { }
  complex(complex const &src) noexcept: tmp(src.tmp), re(src.re), im(src.im) { }
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

//static_assert(IMandelReal<dd_real>, "maybe?");

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
#endif // MANDELMATH_NUMBER_HPP
