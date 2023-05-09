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
  virtual double radixfloor(const number_a &other) const=0; //nearest smaller power of 2 (1.5->1->1)
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
protected:
  BASE store; //either "double" for special optimized template
              //or "double *" for universal implementation
              //going meta by using BASE="number<actual base> *"
  friend class number<double>;
  friend class number<__float128>;
  friend class number<dd_real>;
  friend class number<dq_real>;
  friend class number<real642>;
  friend class number<number_a *>;
public:
  virtual double eps2() const override;
  virtual double eps234() const override;
  number(): number_a(), store() { } //for ShareableViewInfo
  number(const number &x) noexcept;
  number(number &&x) noexcept: store() { x.swap(*this); };
  number(NumberType ntype);
  void constructLateBecauseQtIsAwesome(NumberType ntype);
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
  virtual double radixfloor(const number_a &other) const override;
  double radixfloor(const number<BASE> &other) const; //nearest smaller power of 2 (1.5->1->1)
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
  complex &assign_across(const complex<OTHER_BASE> &src);
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
