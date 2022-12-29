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
  virtual bool isle(const number_a &other)=0; //return store<=other
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
  virtual bool isle(const number_a &other) override; //return store<=other
  bool isle(const number<BASE> &other); //return store<=other
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

#if 0
class worker_multi;

struct number_index
{
  const int index;
  worker_multi *dbg_worker;
  operator int() const { return index; }
  //number_index is const operator =(int index) { this->index=index; }
  number_index(int index, worker_multi *dbg_worker): index(index), dbg_worker(dbg_worker) { }
};

struct complex_index
{
  const int index;
  worker_multi *dbg_worker;
  operator int() const { return index; }
  //number_index is const operator =(int index) { this->index=index; }
  complex_index(int index, worker_multi *dbg_worker): index(index), dbg_worker(dbg_worker) { }
};

//template <class WORKER_MULTI>
class number_place; //always <worker_multi>
//template <class WORKER_MULTI>
class complex_place; //always <worker_multi>

/*struct number_index
{
  int index_into_worker_;
  virtual operator int() { return index_into_worker_; }
#if QT_NO_DEBUG
  number_index(worker_multi *): index_into_worker(-1) {}
  number_index(worker_multi *, int iiw): index_into_worker(iiw) {}
#else
  worker_multi *const dbg_worker;
  number_index(worker_multi *worker): index_into_worker_(-1), dbg_worker(worker) {}
  number_index(worker_multi *worker, int iiw): index_into_worker_(iiw), dbg_worker(worker) {}
#endif
};

template <int INDEX_INTO_WORKER>
struct number_index_fix_: public number_index
{
  //int index_into_worker;
  operator int() override { return INDEX_INTO_WORKER; }
#if QT_NO_DEBUG
  number_index(worker_multi *): index_into_worker(-1) {}
  number_index(worker_multi *, int iiw): index_into_worker(iiw) {}
#else
  worker_multi *const dbg_worker;
  number_index_fix_(worker_multi *worker): dbg_worker(worker) {}
  number_index_fix_(worker_multi *worker, int iiw): dbg_worker(worker) { working_assert(iiw==INDEX_INTO_WORKER); }
#endif
};

template <class WORKER_MULTI>
class number_place;

template <class WORKER_MULTI>
class complex_place;

struct complex_index
{
  int index_into_worker;
#if QT_NO_DEBUG
  complex_index(worker_multi *): index_into_worker(-1) {}
  complex_index(worker_multi *, int iiw): index_into_worker(iiw) {}
  number_index getre() const { return number_index(nullptr, index_into_worker); }
  number_index getim() const { return number_index(nullptr, index_into_worker+1); }
#else
  worker_multi *const dbg_worker;
  complex_index(worker_multi *worker): index_into_worker(-1), dbg_worker(worker) {}
  complex_index(worker_multi *worker, int iiw): //DO NOT turn it to (complex_place ...), we need compiler to know that worker is fixed
    index_into_worker(iiw), dbg_worker(worker) {}
  number_index getre() const { return number_index(dbg_worker, index_into_worker); }
  number_index getim() const { return number_index(dbg_worker, index_into_worker+1); }
#endif
};

template <int INDEX_INTO_WORKER>
struct complex_index_fix
{
#if QT_NO_DEBUG
  complex_index(worker_multi *): index_into_worker(-1) {}
  complex_index(worker_multi *, int iiw): index_into_worker(iiw) {}
  number_index_fix<INDEX_INTO_WORKER> getre() { return number_index(nullptr, INDEX_INTO_WORKER); }
  number_index_fix<INDEX_INTO_WORKER+1> getim() { return number_index(nullptr, INDEX_INTO_WORKER+1); }
#else
  worker_multi *const dbg_worker;
  complex_index_fix(worker_multi *worker): dbg_worker(worker) {}
  complex_index_fix(worker_multi *worker, int iiw): //DO NOT turn it to (complex_place ...), we need compiler to know that worker is fixed
    dbg_worker(worker) { working_assert(iiw==INDEX_INTO_WORKER); }
  number_index_fix_<INDEX_INTO_WORKER> getre() const { return number_index(dbg_worker, INDEX_INTO_WORKER); }
  number_index_fix_<INDEX_INTO_WORKER+1> getim() const { return number_index(dbg_worker, INDEX_INTO_WORKER+1); }
#endif
};



/ * const is no more... maybe in a few years
#if 1 //faster
typedef number_pointer_ number_pointer_c_;
#else //better
union number_pointer_c_
{
  const double *asf64;
#if !ONLY_DOUBLE_WORKER
  const __float128 *asf128;
  const dd_real *asdd;
  const dq_real *asdq;
#endif
  number_pointer_c_(): asf64(nullptr) {}
  number_pointer_c_(double *asf64): asf64(asf64) {}
  number_pointer_c_(const real642 *as642): as642(as642) {}
  number_pointer_c_(const number_pointer_ &src): asf64(src.asf64) {}
};
#endif
*/

/*template <class WORKER_MULTI>
class complex_instance; //forward declaration for friend
template <class WORKER_MULTI, int INDEX_INTO_WORKER>
class complex_instance_fix; //forward declaration for friend*/

class worker_multi
{
public:
  enum Type { typeEmpty
#if NUMBER_DOUBLE_EXISTS
              , typeDouble
#endif
#if !ONLY_DOUBLE_WORKER
              , typeFloat128
              , typeDDouble
              , typeQDouble
              , typeReal642
#endif
            };

  template <class WORKER_MULTI>
  class Allocator
  {
  protected:
    int lowest;
    int allocUpTo;
    int capacity;
    Allocator *above;
  public:
    WORKER_MULTI *worker;
    Allocator(int capacity); //fake constructor
    Allocator(worker_multi *worker, int capacity);
    Allocator(Allocator *allo, int capacity); //alloc a block
    Allocator(Allocator *allo, int lowest, int count); //reuse allocated block
    ~Allocator();
    bool checkFill() { return (allocUpTo==lowest+capacity); }
    bool checkIndex(int index) { return (index>=lowest && index<lowest+allocUpTo); }
    number_index alloc();
    void dealloc(number_index ptr);
    void _getFirstCapac(int &first, int &capacity) const { first=lowest; capacity=this->capacity; }
    friend WORKER_MULTI;
    //friend class complex_instance<WORKER_MULTI>;
  };

protected:
  Type _ntype;
  virtual number_index getNumber(int index)=0; //support for Allocator
  Allocator<worker_multi> allocator;
  //int capacity;

  enum IndexIntoWorker
  {
    iiw_tmp1=0,
    iiw_tmp2=1,
    iiw_tmp3=2,
    iiw_tmp4=3,
    iiw__COUNT=4
  };

  //static constexpr int COUNT_OF_TMPS=4; //tmp1..tmp4
  //virtual void getTmp12(number_index &t1, number_index &t2)=0;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)=0;
  //friend class complex_instance<worker_multi>;
  //template <int IIW>
  //friend class complex_instance_fix<worker_multi, IIW>;

public:
  virtual double eps2() { return 1.232596e-32; }
  virtual double eps234() { return 1.17e-24; }
  worker_multi(Type ntype, int capacity): _ntype(ntype), allocator(this, capacity)/*, capacity(capacity)*/ { }
  virtual ~worker_multi() { /*capacity=0;*/ }
  virtual Type ntype() const { return this->_ntype; }
  static worker_multi *create(Type ntype, int capacity);
  static worker_multi *create(Type ntype, Allocator<worker_multi> *allocator);
  Allocator<worker_multi> *getAllocator() { return &allocator; }
  /*virtual number_pointer_ alloc()=0;
  virtual void alloc_array(int count)=0;
  virtual void dealloc(number_pointer_ store)=0;
  virtual void dealloc_array(int count)=0;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len)=0; //for now, assert(this.ntype==src.ntype)
    //or assign_block(int dst, allocator *src) + assign_block(allocator *dst, int src) + assign_block(int dst, int src, int len)
  virtual void assign_block(int dst, const Allocator<worker_multi> *src)
   { int first, capac; src->_getFirstCapac(first, capac); assign_block(dst, src->worker, first, capac); }

  virtual void zero(const number_index store, double val=0)=0;
  virtual void assign(number_index store, const number_index src)=0;
  virtual void assign_across(number_index store, const number_place *src)=0;
  virtual void chs(const number_index store)=0;
  virtual void lshift(const number_index store, int shoft)=0; // self <<= shoft; 1 lshift -10000 = 0 not error
  virtual void round(const number_index store)=0;
  virtual void frac(const number_index store)=0; //-1<result<1
  virtual void mod1(const number_index store)=0; //0<=result<1
  virtual void add_double(const number_index store, double x)=0;
  virtual void mul_double(const number_index store, double x)=0;
  virtual void add(const number_index store, const number_index other)=0;
  virtual void sub(const number_index store, const number_index other)=0;
  virtual void rsub(const number_index store, const number_index other)=0;
  virtual void mul(const number_index store, const number_index other)=0;
  virtual void sqr(const number_index store)=0;
  virtual double radixfloor(const number_index store1, number_index store2)=0; //nearest smaller power of 2 (1.5->1->1)
  virtual void recip(const number_index store)=0;
  virtual void sqrt(const number_index store)=0;
  virtual int compare(const number_index store, const number_index other)=0; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  virtual bool isequal(const number_index store, const number_index other)=0; //return store==other
  virtual bool is0(const number_index store)=0;
  virtual bool isle(const number_index store, const number_index other)=0; //return store<=other
  virtual bool isle0(const number_index store)=0; //return store<=0
  virtual bool isl0(const number_index store)=0; //return store<0
  virtual bool isl1(const number_index store)=0; //return store<1
  virtual void min(const number_index store, const number_index other)=0;

  virtual QString toString(const number_index store)=0;
  virtual int toRound(const number_index store)=0;
  virtual double toDouble(const number_index store)=0;
};

//TODO: combine worker_multi_xxx into template<dd_real, Type::typeDDouble, 1e-64>
//will need to write a wrapper around double and float128

#if NUMBER_DOUBLE_EXISTS
class worker_multi_double: public worker_multi
{
  double *storage_;
  //double tmp1, tmp2;
  //double tmp3, tmp4;
protected:
  virtual number_index getNumber(int index) override; //support for Allocator
  //virtual void getTmp12(number_index &t1, number_index &t2) override;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4) override;
  friend class Allocator<worker_multi_double>;
  //friend class complex<worker_multi>;
  //friend class complex_instance<worker_multi_double>;
public:
  double eps2() override { return 1.232596e-32;  /* 2^-(2*53) rounded up */ }
  double eps234() override { return 1.17e-24;  /* eps2^(3/4) */ }

  worker_multi_double(int capacity): worker_multi(Type::typeDouble, capacity),
      storage_(new double[iiw__COUNT+capacity]) { }//for (int i=0; i<capacity; i++) storage[i]=std::numeric_limits<double>::signaling_NaN(); }
  //worker_multi_double(worker_multi *source);
  worker_multi_double(Allocator<worker_multi> *source);
  virtual ~worker_multi_double() override;
  virtual Type ntype() const override { return typeDouble; }
  /*virtual number_pointer_ alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer_ store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)
  const double *_getStorage_() { return storage_; }

  void zero(const number_index store, double val=0) override;
  void assign(const number_index store, const number_index src) override;// { store->assign<number_double>(*src); };
  void assign_across(const number_index store, const number_place *src) override;
  void chs(const number_index store) override;
  void lshift(const number_index store, int shoft) override;
  void round(const number_index store) override;
  void frac(const number_index store) override; //-1<result<1
  void mod1(const number_index store) override; //0<=result<1
  void add_double(const number_index store, double x) override;
  void mul_double(const number_index store, double x) override;
  void add(const number_index store, const number_index other) override;
  void sub(const number_index store, const number_index other) override;
  void rsub(const number_index store, const number_index other) override;
  void mul(const number_index store, const number_index other) override;
  void sqr(const number_index store) override;
  double radixfloor(const number_index store1, number_index store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_index store) override;
  void sqrt(const number_index store) override;
  int compare(const number_index store, const number_index other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_index store, const number_index other) override;
  bool is0(const number_index store) override;
  bool isle(const number_index store, const number_index other) override;
  bool isle0(const number_index store) override;
  bool isl0(const number_index store) override;
  bool isl1(const number_index store) override;
  void min(const number_index store, const number_index other) override;

  QString toString(const number_index store) override;
  int toRound(const number_index store) override;
  double toDouble(const number_index store) override;
};
#endif //NUMBER_DOUBLE_EXISTS

#if !ONLY_DOUBLE_WORKER
class worker_multi_float128: public worker_multi
{
  __float128 *storage_;
  //__float128 tmp1, tmp2;
  //__float128 tmp3, tmp4;
protected:
  virtual number_index getNumber(int index) override; //support for Allocator
  //virtual void getTmp12(number_index &t1, number_index &t2) override;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4) override;
  friend class Allocator<worker_multi_float128>;
  //friend class complex_instance<worker_multi_float128>;
public:
  double eps2() override { return 9.27302e-69;  /* 2^-(2*113) rounded up */ }
  double eps234() override { return 9.45e-52;  /* eps2^(3/4) */ }

  worker_multi_float128(int capacity): worker_multi(Type::typeFloat128, capacity),
      storage_(new __float128[iiw__COUNT+capacity]) {}
  worker_multi_float128(Allocator<worker_multi> *source);
  virtual ~worker_multi_float128();
  virtual Type ntype() const override { return typeFloat128; }
  /*virtual number_pointer_ alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer_ store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)
  const __float128 *_getStorage_() { return storage_; }

  void zero(const number_index store, double val=0) override;
  void assign(const number_index store, const number_index src) override;// { store->assign<number_double>(*src); };
  void assign_across(const number_index store, const number_place *src) override;
  void chs(const number_index store) override;
  void lshift(const number_index store, int shoft) override;
  void round(const number_index store) override;
  void frac(const number_index store) override; //-1<result<1
  void mod1(const number_index store) override; //0<=result<1
  void add_double(const number_index store, double x) override;
  void mul_double(const number_index store, double x) override;
  void add(const number_index store, const number_index other) override;
  void sub(const number_index store, const number_index other) override;
  void rsub(const number_index store, const number_index other) override;
  void mul(const number_index store, const number_index other) override;
  void sqr(const number_index store) override;
  double radixfloor(const number_index store1, number_index store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_index store) override;
  void sqrt(const number_index store) override;
  int compare(const number_index store, const number_index other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_index store, const number_index other) override;
  bool is0(const number_index store) override;
  bool isle(const number_index store, const number_index other) override;
  bool isle0(const number_index store) override;
  bool isl0(const number_index store) override;
  bool isl1(const number_index store) override;
  void min(const number_index store, const number_index other) override;

  QString toString(const number_index store) override;
  int toRound(const number_index store) override;
  double toDouble(const number_index store) override;
};

class worker_multi_ddouble: public worker_multi
{
protected:
  dd_real *storage_;
  //dd_real tmp1, tmp2;
  //dd_real tmp3, tmp4;
protected:
  virtual number_index getNumber(int index) override; //support for Allocator
  //virtual void getTmp12(number_index &t1, number_index &t2) override;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4) override;
  friend class Allocator<worker_multi_ddouble>;
  //friend class complex_instance<worker_multi_ddouble>;
public:
  double eps2() override { return 6.07717e-64; /* 2^-(2*(53+52)) rounded up */ }
  double eps234() override { return 3.87e-48;  /* eps2^(3/4) */ }

  worker_multi_ddouble(int capacity): worker_multi(Type::typeDDouble, capacity),
    storage_(new dd_real[iiw__COUNT+capacity]) {}
  worker_multi_ddouble(Allocator<worker_multi> *source);
  virtual ~worker_multi_ddouble();
  //{ assert((store->dbgType==Type::typeDDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() const override { return typeDDouble; }
  /*virtual number_pointer_ alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer_ store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)
  const dd_real *_getStorage_() { return storage_; }

  void zero(const number_index store, double val=0) override;
  void assign(const number_index store, const number_index src) override;// { store->assign<number_double>(*src); };
  void assign_across(const number_index store, const number_place *src) override;
  void chs(const number_index store) override;
  void lshift(const number_index store, int shoft) override;
  void round(const number_index store) override;
  void frac(const number_index store) override; //-1<result<1
  void mod1(const number_index store) override; //0<=result<1
  void add_double(const number_index store, double x) override;
  void mul_double(const number_index store, double x) override;
  void add(const number_index store, const number_index other) override;
  void sub(const number_index store, const number_index other) override;
  void rsub(const number_index store, const number_index other) override;
  void mul(const number_index store, const number_index other) override;
  void sqr(const number_index store) override;
  double radixfloor(const number_index store1, number_index store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_index store) override;
  void sqrt(const number_index store) override;
  int compare(const number_index store, const number_index other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_index store, const number_index other) override;
  bool is0(const number_index store) override;
  bool isle(const number_index store, const number_index other) override;
  bool isle0(const number_index store) override;
  bool isl0(const number_index store) override;
  bool isl1(const number_index store) override;
  void min(const number_index store, const number_index other) override;

  QString toString(const number_index store) override;
  int toRound(const number_index store) override;
  double toDouble(const number_index store) override;
};

class worker_multi_qdouble: public worker_multi
{
protected:
  dq_real *storage_;
  //dq_real tmp1, tmp2;
  //dq_real tmp3, tmp4;
protected:
  virtual number_index getNumber(int index) override; //support for Allocator
  //virtual void getTmp12(number_index &t1, number_index &t2) override;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4) override;
  friend class Allocator<worker_multi_qdouble>;
  //friend class complex_instance<worker_multi_qdouble>;
public:
  double eps2() override { return 6.07717e-64; /* 2^-(2*(53+52)) */ }
  double eps234() override { return 3.87e-48;  /* eps2^(3/4) */ }
  worker_multi_qdouble(int capacity): worker_multi(Type::typeQDouble, capacity),
    storage_(new dq_real[iiw__COUNT+capacity]) {}
  worker_multi_qdouble(Allocator<worker_multi> *source);
  virtual ~worker_multi_qdouble();
  //{ assert((store->dbgType==Type::typeDDouble) ||
  //         (store->dbgType==Type::typeEmpty)); }
  virtual Type ntype() const override { return typeQDouble; }
  /*virtual number_pointer_ alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer_ store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)
  const dq_real *_getStorage_() { return storage_; } //TODO: return storage_+COUNT_OF_TMPS

  void zero(const number_index store, double val=0) override;
  void assign(const number_index store, const number_index src) override;// { store->assign<number_double>(*src); };
  void assign_across(const number_index store, const number_place *src) override;
  void chs(const number_index store) override;
  void lshift(const number_index store, int shoft) override;
  void round(const number_index store) override;
  void frac(const number_index store) override; //-1<result<1
  void mod1(const number_index store) override; //0<=result<1
  void add_double(const number_index store, double x) override;
  void mul_double(const number_index store, double x) override;
  void add(const number_index store, const number_index other) override;
  void sub(const number_index store, const number_index other) override;
  void rsub(const number_index store, const number_index other) override;
  void mul(const number_index store, const number_index other) override;
  void sqr(const number_index store) override;
  double radixfloor(const number_index store1, number_index store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_index store) override;
  void sqrt(const number_index store) override;
  int compare(const number_index store, const number_index other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_index store, const number_index other) override;
  bool is0(const number_index store) override;
  bool isle(const number_index store, const number_index other) override;
  bool isle0(const number_index store) override;
  bool isl0(const number_index store) override;
  bool isl1(const number_index store) override;
  void min(const number_index store, const number_index other) override;

  QString toString(const number_index store) override;
  int toRound(const number_index store) override;
  double toDouble(const number_index store) override;
};

class worker_multi_real642: public worker_multi
{
  real642 *storage_;
  //real642 tmp1, tmp2;
  //real642 tmp3, tmp4;
protected:
  virtual number_index getNumber(int index) override; //support for Allocator
  //virtual void getTmp12(number_index &t1, number_index &t2) override;
  //virtual void getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4) override;
  friend class Allocator<worker_multi_real642>;
  //friend class complex_instance<worker_multi_real642>;
public:
  double eps2() override { return 9.27302e-69;  /* 2^-(2*113) rounded up */ }
  double eps234() override { return 9.45e-52;  /* eps2^(3/4) */ }

  worker_multi_real642(int capacity): worker_multi(Type::typeReal642, capacity),
      storage_(new real642[iiw__COUNT+capacity]) {}
  worker_multi_real642(Allocator<worker_multi> *source);
  virtual ~worker_multi_real642();
  virtual Type ntype() const override { return typeReal642; }
  /*virtual number_pointer_ alloc() override;
  virtual void alloc_array(int count) override;
  virtual void dealloc(number_pointer_ store) override;
  virtual void dealloc_array(int count) override;*/
  virtual void assign_block(int dst, worker_multi *src_worker, int src, int len) override; //for now, assert(this.ntype==src.ntype)
  const real642 *_getStorage_() { return storage_; }

  void zero(const number_index store, double val=0) override;
  void assign(const number_index store, const number_index src) override;// { store->assign<number_double>(*src); };
  void assign_across(const number_index store, const number_place *src) override;
  void chs(const number_index store) override;
  void lshift(const number_index store, int shoft) override;
  void round(const number_index store) override;
  void frac(const number_index store) override; //-1<result<1
  void mod1(const number_index store) override; //0<=result<1
  void add_double(const number_index store, double x) override;
  void mul_double(const number_index store, double x) override;
  void add(const number_index store, const number_index other) override;
  void sub(const number_index store, const number_index other) override;
  void rsub(const number_index store, const number_index other) override;
  void mul(const number_index store, const number_index other) override;
  void sqr(const number_index store) override;
  double radixfloor(const number_index store1, number_index store2) override; //nearest smaller power of 2 (1.5->1->1)
  void recip(const number_index store) override;
  void sqrt(const number_index store) override;
  int compare(const number_index store, const number_index other) override; //return -1 if <, 0 if =, +1 if >
  bool isequal(const number_index store, const number_index other) override;
  bool is0(const number_index store) override;
  bool isle(const number_index store, const number_index other) override;
  bool isle0(const number_index store) override;
  bool isl0(const number_index store) override;
  bool isl1(const number_index store) override;
  void min(const number_index store, const number_index other) override;

  QString toString(const number_index store) override;
  int toRound(const number_index store) override;
  double toDouble(const number_index store) override;
};

#endif //ONLY_DOUBLE_WORKER

/*
New rules:
  - .h files only contain number_pointer and ..._place, never _instance
  - ..._place = worker+number_pointer
  -     in function calls, _place is only used across workers
  -     used in class fields
*/
//template <class WORKER_MULTI>
//class number_instance;

//template <class WORKER_MULTI>
class number_place
{
protected:
  //worker_multi *worker;
  worker_multi::Allocator<worker_multi> *allocator;
  //friend class number_instance<WORKER_MULTI>;
  friend worker_multi;
  friend worker_multi_double;
  friend worker_multi_float128;
  friend worker_multi_ddouble;
  friend worker_multi_qdouble;
  friend worker_multi_real642;
public:
  number_index ptr;
  //or maybe implement cast to number_pointer_c_ but I don't like that
  number_place(worker_multi *worker): allocator(/*(worker_multi::Allocator<WORKER_MULTI> *)*/worker->getAllocator()),
    ptr(allocator->alloc()) {}
  number_place(worker_multi::Allocator<worker_multi> *allocator): allocator(/*(worker_multi::Allocator<WORKER_MULTI> *)*/allocator),
    ptr(allocator->alloc()) {}
  //number(number_pointer_ *src); //swap into me, and clear src, src->freeond=false because already freed
  ~number_place()
  {
    //if (free_storage_on_destroy)
    {
      allocator->dealloc(ptr);
    };
  }
  bool checkWorker(worker_multi *worker) { return worker==allocator->worker; }
};

/*template <class WORKER_MULTI>
class number_instance
{
protected:
  WORKER_MULTI *const worker;
  const number_index _ptr;
  friend class complex_instance<WORKER_MULTI>;
  number_instance(WORKER_MULTI *worker, number_index ptr): worker(worker), _ptr(ptr)
   { }
public:
  //static constexpr int iiw=INDEX_INTO_WORKER;
  number_instance(number_place<WORKER_MULTI> *place): worker(place->allocator->worker), _ptr(place->ptr)
   { }
  number_index getptr() const { return _ptr; }
  void zero(double v=0);
  void assign(const number_index src);
  //void assign_across(const worker_multi *src_worker, const number_pointer_c_ src);
  void assign_across(const number_place<WORKER_MULTI> *src);
  void assign_across(const number_place<WORKER_MULTI> &&src);
  //would be nicer to have assign_across_from_complex_re(const complex *src)... but not by much
  void lshift(int shoft);
  void add(const number_index other);
  void chs();
  void sub(const number_index other);
  void sqr();
  void mul(const number_index other);
  void recip();
  //void recip_prepared();
  void sqrt();
  void add_double(double x);
  bool is0() const;
  void min(const number_index other);
  QString toString() const { return worker->toString(_ptr); }
  int toRound() const;
  double toDouble() const;
};

template <class WORKER_MULTI, int INDEX_INTO_WORKER>
class number_instance_fix
{
protected:
  WORKER_MULTI *const worker;
public:
  //static constexpr int iiw=INDEX_INTO_WORKER;
  //setting worker directly, rather than from place, in hopes of better optimization
  number_instance_fix(WORKER_MULTI *const worker, number_place<WORKER_MULTI> *place): worker(worker)
   {
    working_assert(place->checkWorker(worker));
    if (place->ptr.index_into_worker!=INDEX_INTO_WORKER)
      dbgPointII(place->ptr.index_into_worker, INDEX_INTO_WORKER);
    working_assert(place->ptr.index_into_worker==INDEX_INTO_WORKER);
   }
  number_index_fix_<INDEX_INTO_WORKER> getptr() const { return number_index_fix_<INDEX_INTO_WORKER>(worker); }
  void zero(double v=0)
   { worker->zero(number_pointer_(worker, INDEX_INTO_WORKER), v); }
  void assign(const number_index src)
   { worker->assign(number_pointer_(worker, INDEX_INTO_WORKER), src); }
  //void assign_across(const worker_multi *src_worker, const number_pointer_c_ src)
   //{ worker->assign_across(number_pointer_(worker, INDEX_INTO_WORKER), src_worker, src); }
  void assign_across(const number_place<WORKER_MULTI> *src)
   { worker->assign_across(number_pointer_(worker, INDEX_INTO_WORKER), src->allocator->worker, src->ptr); }
  void lshift(int shoft)
   { worker->lshift(number_pointer_(worker, INDEX_INTO_WORKER), shoft); }
  void add(const number_index other)
   { worker->add(number_pointer_(worker, INDEX_INTO_WORKER), other); }
  void chs()
   { worker->chs(number_pointer_(worker, INDEX_INTO_WORKER)); }
  void sub(const number_index other)
   { worker->sub(number_pointer_(worker, INDEX_INTO_WORKER), other); }
  void sqr()
   { worker->sqr(number_pointer_(worker, INDEX_INTO_WORKER)); }
  void mul(const number_index other)
   { worker->mul(number_pointer_(worker, INDEX_INTO_WORKER), other); }
  void recip()
   { worker->recip(number_pointer_(worker, INDEX_INTO_WORKER)); }
  void sqrt()
   { worker->sqrt(number_pointer_(worker, INDEX_INTO_WORKER)); }
  void add_double(double x)
   { worker->add_double(number_pointer_(worker, INDEX_INTO_WORKER), x); }
  bool is0()
   { return worker->is0(number_pointer_(worker, INDEX_INTO_WORKER)); }
  void min(const number_index other)
   { worker->min(number_pointer_(worker, INDEX_INTO_WORKER), other); }
  QString toString() const
   { return worker->toString(getptr()); }
  int toRound() const
   { return worker->toRound(number_pointer_(worker, INDEX_INTO_WORKER)); }
  double toDouble() const
   { return worker->toDouble(number_pointer_(worker, INDEX_INTO_WORKER)); }
};
*/

//template <class WORKER_MULTI>
class complex_place
{
protected:
  //worker_multi *worker;
  worker_multi::Allocator<worker_multi> *allocator;

public:
  const number_place re;
  const number_place im;
  //complex(number *re, number *im, number *tmp1, number *tmp2): tmp1(tmp1), tmp2(tmp2), re(re), im(im) { }
  complex_place(worker_multi *worker): allocator(/*(worker_multi::Allocator<WORKER_MULTI> *)*/worker->getAllocator()),
    re(allocator), im(allocator) {}
  complex_place(worker_multi::Allocator<worker_multi> *allocator): allocator(/*(worker_multi::Allocator<WORKER_MULTI> *)*/allocator),
    re(allocator), im(allocator) {}
#if QT_NO_DEBUG
  complex_index ptr() const { return complex_index(nullptr, re); }
#else
  complex_index ptr() const { return complex_index(re.ptr, allocator->worker); }
#endif
  //number_place<WORKER_MULTI> re_place() const { return number_place<WORKER_MULTI>(allocator->worker, re); }
  //number_place<WORKER_MULTI> im_place() const { return number_place<WORKER_MULTI>(allocator->worker, im); }
  ~complex_place()
  {
    //allocator->dealloc(im);
    //allocator->dealloc(re);
  }
};

/*template <class WORKER_MULTI>
class complex_instance
{
//protected:
  //in .re and in .im WORKER_MULTI *const worker;
  //const number_pointer_ re;
  //const number_pointer_ im;
public:
  number_instance<WORKER_MULTI> re;
  number_instance<WORKER_MULTI> im;
  complex_index getptr() const { return complex_index(re.worker, re._ptr.index_into_worker); }
  complex_instance(complex_place<WORKER_MULTI> *place): //worker(place->allocator->worker),
    re(place->allocator->worker, place->re), im(place->allocator->worker, place->im) { }
  void zero(double r=0, double i=0);
  void assign(const complex_index other);
  void assign_across(const complex_place<WORKER_MULTI> *src);
  void lshift(int shoft);
  void add(const complex_index other);
  void chs();
  void sub(const complex_index other); //this=this-other
  void rsub(const complex_index other); //this=other-this
  //add_double(r), add_double(r, i)
  void sqr();
  void mul(const complex_index other);
  //mul(number or number_pointer_c_) maybe called "scale"
  //or mul_double, div_double using tmp (usually mul_int, div_int but also *(1/m-1/n))
  void recip();
  void recip_prepared();
  void sqrt();
  void root_approx(int n); //this^(1/n), most real root; only in double precision
  void sign(); //this/=sqrt(this.mag())
  double getMag_double() const;
  const number_index getMag_tmp() const;
  const number_index getMag1_tmp() const;
  const number_index getDist1_tmp() const;
  const number_index mulreT_tmp(const complex_index other) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  double dist2_double(const complex_index other) const;
  const number_index dist2_tmp(const complex_index other) const;
  bool isequal(const complex_index other) const;
  //assign_re(number_pointer_) or assign_re(number) & assign_re_re(complex) & assign_re_im(complex)
  bool is0() const;
  int mag_cmp_1(); //(mag() <=> 1) -> -1,0,+1
  QString toString();
};

template <class WORKER_MULTI, int INDEX_INTO_WORKER>
class complex_instance_fix
{
protected:
  WORKER_MULTI *const worker;
  //const number_pointer_ re;
  //const number_pointer_ im;
public:
  number_instance_fix<WORKER_MULTI, INDEX_INTO_WORKER> re;
  number_instance_fix<WORKER_MULTI, INDEX_INTO_WORKER+1> im;
  complex_instance_fix(WORKER_MULTI *worker, complex_place<WORKER_MULTI> *place): worker(worker),
    re(place->allocator->worker, place->re_place()), im(place->allocator->worker, place->im_place()) { }
  complex_index_fix<INDEX_INTO_WORKER> getptr() const { return complex_index_fix<INDEX_INTO_WORKER>(re.worker, INDEX_INTO_WORKER); }
  void zero(double r=0, double i=0);
  //template <int INDEX_INTO_OTHER>
  void assign(const complex_index other);
  template <int INDEX_INTO_OTHER>
  void assign(const complex_index_fix<INDEX_INTO_OTHER> other)
   { re.assign(other.getre()); im.assign(other.getim()); }
  //template <int INDEX_INTO_OTHER>
  void assign_across(const complex_place<WORKER_MULTI> *other);
  void lshift(int shoft);
  //template <int INDEX_INTO_OTHER>
  void add(const complex_index other);
  void chs();
  //template <int INDEX_INTO_OTHER>
  void sub(const complex_index other); //this=this-other
  //template <int INDEX_INTO_OTHER>
  void rsub(const complex_index other); //this=other-this
  //add_double(r), add_double(r, i)
  void sqr();
  void mul(const number_index other)
   { re.mul(other); im.mul(other); }
  //template <int INDEX_INTO_OTHER>
  void mul(const complex_index other);
  //mul(number or number_pointer_c_) maybe called "scale"
  //or mul_double, div_double using tmp (usually mul_int, div_int but also *(1/m-1/n))
  void recip();
  void recip_prepared();
  void sqrt();
  void root_approx(int n); //this^(1/n), most real root; only in double precision
  void sign(); //this/=sqrt(this.mag())
  double getMag_double() const;
  const number_index getMag_tmp() const;
  const number_index getMag1_tmp() const;
  const number_index getDist1_tmp() const;
  //template <int INDEX_INTO_OTHER>
  const number_index mulreT_tmp(const complex_index other) const; //Re(this*conjugate(other)) = re*o->re+im*o->im
  //template <int INDEX_INTO_OTHER>
  double dist2_double(const complex_index other) const;
  //template <int INDEX_INTO_OTHER>
  const number_index dist2_tmp(const complex_index other) const;
  //template <int INDEX_INTO_OTHER>
  bool isequal(const complex_index other) const;
  //assign_re(number_pointer_) or assign_re(number) & assign_re_re(complex) & assign_re_im(complex)
  bool is0() const;
  int mag_cmp_1(); //(mag() <=> 1) -> -1,0,+1
  QString toString();
};
*/

#endif

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
