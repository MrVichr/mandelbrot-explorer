#ifndef MANDELMATH_NUMBER_HPP
#define MANDELMATH_NUMBER_HPP

#include <cstdint>
#include <QString>

#include "double_double.hpp"
#include "multiprec.hpp"

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

typedef dd_real dq_real;

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
  virtual bool reduceAngle()=0; //one step towards -pi<=x<pi, true if changed
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

template <typename BASE>
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
  virtual bool reduceAngle() override;
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

template <typename BASE>
class complex
{
public:
  struct Scratchpad
  {
    number<BASE> tmp1;
    number<BASE> tmp2;
    number<BASE> tmp3;
    number<BASE> tmp4;
    Scratchpad(): tmp1(), tmp2(), tmp3(), tmp4() { }
    Scratchpad(NumberType ntype): tmp1(ntype), tmp2(ntype), tmp3(ntype), tmp4(ntype) { }
  };
  number<BASE> re;
  number<BASE> im;
  complex(): re(), im() { }
  complex(NumberType ntype): re(ntype), im(ntype) { }
  //complex(BASE re, BASE im): re(re), im(im) { }
  void readFrom(void *storage, int index); //BASE storage[index], [index+1]
  void writeTo(void *storage, int index) const; //BASE storage[index], [index+1]
  void zero(double r=0, double i=0);
  void assign(const complex *other);
  template <typename OTHER_BASE>
  void assign_across(const complex<OTHER_BASE> *src);
  //void assign_across(const BASE re, const BASE im);
  void lshift(int shoft);
  //add_double(r), add_double(r, i)
  void mul_double(double m);
  void add(const complex *const other);
  void chs();
  void sub(const complex *other); //this=this-other
  void rsub(const complex *other); //this=other-this
  void sqr(Scratchpad *tmp);
  void mul(const number<BASE> &other); //maybe called "scale"
  void mul(const complex *other, Scratchpad *tmp);
  void recip(Scratchpad *tmp);
  void recip_prepared(Scratchpad *tmp);
  void sqrt(Scratchpad *tmp);
  void pow_int(int n, Scratchpad *tmp); //this^n
  void root_approx(int n); //this^(1/n), most real root; only in double precision
  void ln_approx(); //this:=ln(this), only in double precision
  void exp_approx(); //this:=exp(this), only in double precision
  void sign(Scratchpad *tmp); //this/=sqrt(this.mag())
  void cossin(number<BASE> &angle, Scratchpad *tmp);
  void arctan2(number<BASE> *result, Scratchpad *tmp) const;
  double getMag_double() const;
  const number<BASE> *getMag_tmp(Scratchpad *tmp) const;
  const number<BASE> *getMag1_tmp(Scratchpad *tmp) const;
  const number<BASE> *getDist1_tmp(Scratchpad *tmp) const;
  const number<BASE> *mulreT_tmp(const complex *other, Scratchpad *tmp) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  const number<BASE> *ccw_tmp(const complex *other, Scratchpad *tmp) const; //"counterclockwise" Im(other/this)*|this|^2 = re*o.im-im*o.re
  double dist2_double(const complex *other, Scratchpad *tmp) const;
  const number<BASE> *dist2_tmp(const complex *other, Scratchpad *tmp) const;
  bool isequal(const complex *other) const;
  bool is0() const;
  bool isNegative() const; //im<0 || im==0 && re<0
  int mag_cmp_1(Scratchpad *tmp) const; //(mag() <=> 1) -> -1,0,+1, better return (double)(mag()-1)
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
