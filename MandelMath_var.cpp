//#define working_assert(x) { if (!(x)) dbgPoint(); }
#include "MandelMath.hpp"
#include <math.h>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstring>

//#include <cassert>
#include <cmath>

void doNothing(int &x)
{
  //x++;
  (void)x;
}

void doNothing(void *x)
{
  (void)x;
}

void nop() { }

void dbgPoint()
{
  int x=3;
  doNothing(x);
}

void dbgPointII(int x1, int x2)
{
  (void)x1;
  (void)x2;
  dbgPoint();
}

namespace MandelMath {

/* moved to .h
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
  [[unlikely]]
  return decltype(callable((NthTypeOf<0, Options...> &)vari, params...)) ();
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
  [[unlikely]]
  return decltype(callable((NthTypeOf<0, Options...> const &)vari, params...)) ();
}
*/
int gcd(int m, int n)
{
  if (m==n)
    return m;
  if (m==0)
    return n;
  if (n==0)
    return m;
  int c=0;
  while (((m|n)&1)==0)
  {
    m>>=1; n>>=1; c++;
  }

  /* Dividing n by 2 until n becomes odd */
  while ((n & 1) == 0)
    n >>= 1;

  /* From here on, 'n' is always odd. */
  do
  {
    /* If m is even, remove all factor of 2 in m */
    while ((m & 1) == 0)
      m >>= 1;

    /* Now m and n are both odd.
       Swap if necessary so n <= m,
       then set m = m - n (which is even).*/
    if (n > m)
    { int x=n; n=m; m=x; }

    m = (m - n)>>1;
  } while (m != 0);

  /* restore common factors of 2 */
  return n << c;
}

int ctz16(int x)
{
  int ctzidx=(((0xF65*(x&-x))&0x7800)>>10);
  int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
  return ctz1;
  /*
   * first find deBruijn sequence of length ctz_capacity
  shift, discard top bit and append i, rotate to be max; is i followed by 0 only (or last)? -> negate else not negate the discarded bit -> append
      i is considered 1-epsilon (1 except when a tie with another 1)
  ctz 8 bit
  00011101
  000-00i-i00-y   1       old-with_i-rotated-negate
   001-01i-1i0-y   1
    011-11i-11i-y   1
     111-11i-11i-y   0
      110-10i-i10-n   1
       101-01i-1i0-y   0
        010-10i-i10-n   0
         100-00i-i00-y   0
          000
   00011101=0x1D
   0x1D * 0=0x0000 >>4&7=0                      .......-
   0x1D * 1=0x001D >>4&7=1                      ......0.
   0x1D * 2=0x003A >>4&7=3   MAGIC[3]:=ctz(2)   ....1...
   0x1D * 4=0x0074 >>4&7=7   MAGIC[7]:=ctz(4)   2.......
   0x1D * 8=0x00E8 >>4&7=6                      .3......
   0x1D *16=0x01D0 >>4&7=5                      ..4.....
   0x1D *32=0x03A0 >>4&7=2                      .....5..
   0x1D *64=0x0740 >>4&7=4                      ...6....
   0x1D*128=0x0E80 >>4&7=0   MAGIC[0]=ctz(128)  .......7
   0x1D*256=0x1D00 >>4&7=0
   (0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07 ?= ctz(0..255)

   ctz 16 bit
   0000-000i-i000-y 1
    0001-001i-1i00-y 1
     0011-011i-11i0-y 1
      0111-111i-111i-y 1
       1111-111i-111i-y 0
        1110-110i-i110-n 1
         1101-101i-1i10-n 1
          1011-011i-11i0-y 0
           0110-110i-i110-n 0
            1100-100i-i100-n 1
             1001-001i-1i00-y 0
              0010-010i-(10i0 or i010)-(y or n)-(1 or 0)  (10i0 correct)
               0101-101i-1i10-n 0
                1010-010i-10i0 or i010-y or n-0 or 1 (10i0 correct)
                 0100-100i-i100-n 0
                  1000-000i-i000-y 0
                   0000
   0000111101100101 = 0xF65
   (0x... >> (((0xF65*(i&-i))&0x7800)>>10))&0x0F ?= ctz(0..255)
         FEDCBA9876543210
   MAGIC=34586C9E27BD1A0F
   MAGIC23=0x59EC6C8C; //00 0101 1001 1110 1100 0110 1100 1000 1100
   MAGIC01=0xC486BD63; // 1100 0100 1000 0110 1011 1101 0110 0011
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);

   ctz 32 bit
   00000-0000i-i0000-y 1
    00001-0001i-1i000-y 1                                        0
     00011-0011i-11i00-y 1                                       1
      00111-0111i-111i0-y 1                                      2
       01111-1111i-1111i-y 1                                     3
        11111-1111i-1111i-y 0                                    4
         11110-1110i-i1110-n 1                                   5
          11101-1101i-1i110-n 1                                  6
           11011-1011i-11i10-n 1                                 7
            10111-0111i-111i0-y 0                                8
             01110-1110i-i1110-n 0                               9
              11100-1100i-i1100-n 1                              a
               11001-1001i-1i100-n 1                             b
                10011-0011i-11i00-y 0                            c
                 00110-0110i-110i0-y 1
                  01101-1101i-1i110-n 0
                   11010-1010i-i1010-n 1
                    10101-0101i-1i010-n 1
                     01011-1011i-11i10-n 0
                      10110-0110i-110i0-y 0
                       01100-1100i-i1100-n 0
                        11000-1000i-i1000-n 1
                         10001-0001i-1i000-y 0
                          00010-0010i-10i00-y 1
                           00101-0101i-1i010-n 0
                            01010-1010i-i1010-n 0
                             10100-0100i-i0100-n 1
                              01001-1001i-1i100-n 0
                               10010-0010i-10i00-y 0
                                00100-0100i-i0100-n 0
                                 01000-1000i-i1000-n 0
                                  10000-0000i-i0000-y 0
                                   00000
   00000111110111001101011000101001 = 0000 0111 1101 1100 1101 0110 0010 1001 = 0x7DCD629
         1f 1e 1d 1c 1b 1a 19 18 17 16 15 14 13 12 11 10 0f 0e 0d 0c 0b 0a 09 08 07 06 05 04 03 02 01 00
   MAGIC= 4  5  6  a  7  f  b 14  8 12 10 19  c 1b 15 1e  3  9  e 13 11 18 1a 1d  2  d 17 1c  1 16  0 1f
   MAGIC4=00000001011101110001111100110101 = 0000 0001 0111 0111 0001 1111 0011 0101 0000 = 0x1771f350
   MAGIC3=00010110100111010110011101010001 = 000 1011 0100 1110 1011 0011 1010 1000 1000 =  0xb4eb3a88
   MAGIC2=11101101000010110010000101110101 = .. 1110 1101 0000 1011 0010 0001 0111 0101 = 0xed0b2175; 0x3b42c85d4 won't fit
   MAGIC1=00111110010001011011001010100101 = 0 0111 1100 1000 1011 0110 0101 0100 1010 = 0x7c8b654a
   MAGIC0=01001110000101101101100101101001 = 0100 1010 0001 0110 1101 1001 0110 1001 = 0x4e16d969

    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
   */

  /*for (int i=0; i<256; i++) //0 returns 7
  {
    int ctz1=(0x23461507 >> (((0x1D*(i&-i))&0x70)>>2))&0x07;
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else ctz2=7;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int i=0; i<65536; i++) //0 returns 15
  {
    int ctzidx=(((0xF65*(i&-i))&0x7800)>>10);
    int ctz1=((0x59EC6C8C >> ctzidx)&0x0C) | ((0xC486BD63 >> ctzidx)&0x03);
    int ctz2=0;
    if (i&1) ctz2=0;
    else if (i&0x03) ctz2=1;
    else if (i&0x07) ctz2=2;
    else if (i&0x0F) ctz2=3;
    else if (i&0x1F) ctz2=4;
    else if (i&0x3F) ctz2=5;
    else if (i&0x7F) ctz2=6;
    else if (i&0xFF) ctz2=7;
    else if (i&0x1FF) ctz2=8;
    else if (i&0x3FF) ctz2=9;
    else if (i&0x7FF) ctz2=10;
    else if (i&0xFFF) ctz2=11;
    else if (i&0x1FFF) ctz2=12;
    else if (i&0x3FFF) ctz2=13;
    else if (i&0x7FFF) ctz2=14;
    else if (i&0xFFFF) ctz2=15;
    else ctz2=15;
    if (ctz2!=ctz1)
      dbgPoint();
  }

  for (int ii=0; ii<32; ii++)
  {
    unsigned int i=(1u<<ii);
    int ctzidx=(((0x7DCD629*(i&-i))&0x7C000000)>>26);
    int ctz1=((0x1771f350 >> ctzidx)&0x10) |
             ((0xb4eb3a88 >> ctzidx)&0x08) |
             ((0xed0b2175 >> ctzidx)&1)<<2 |
             ((0x7c8b654a >> ctzidx)&0x02) |
             ((0x4e16d969 >> ctzidx)&0x01);
    if (ii==32)
    {
      if (ctz1!=31)
        dbgPoint();
    }
    else if (ctz1!=ii)
      dbgPoint();
  }*/
}

bool is2kof(int big, int small) //big == small*2^k, assuming big>=small
{
  //wrong if big=1 && small=2
  big>>=ctz16(big);
  big>>=ctz16(big);
  small>>=ctz16(small);
  return big==small;
}

template <>
int ReverseBits<7, 1>(int val) //reverse bottom 7 bits in blocks of 1
{
  int rh=0x73516240>>((val&7)<<2);
  int rl=0x73516240>>((val&0x70)>>2);
  return ((rh&0x07)<<4) | (val&0x08) | (rl&0x07);
}

/*template <int total, int block>
int ReverseBits(int val)
{
  static_assert(false, "no can do");
}*/


/*double two_pow_n(unsigned int n)
{
  working_assert(n<4096);
  double result=1;
  while (n>=31)
  {
    result*=(1u<<31);
    n-=31;
  }
  if (n>0)
    result*=(1u<<n);
  return result;
}*/

//template<> class number<double>;
//template<> class number<number_a *>;

int CreatedInstancesOfNumber=0;

/*number_a::number_a()
{
  CreatedInstancesOfNumber++;
}

number_a::~number_a()
{
  CreatedInstancesOfNumber--;
  if (CreatedInstancesOfNumber<=0)
    nop(); //should be hit exactly once on program end
}*/

//these 2 have to be in front... maybe move all <number_a *> ?
//...no more?

double *number_traits<double>::convert_block(NumberType old_type, const void *old_data, int count)
{
  double *result=new double[count];
  switch (old_type)
  {
    case NumberType::typeEmpty:
    {
      for (int i=0; i<count; i++)
        result[i]=0;
    } break;
    case NumberType::typeDouble:
    {
      const double *src_storage=(const double *)old_data;
      //TODO: memmove but we'll see later
      for (int i=0; i<count; i++)
      {
        result[i]=*src_storage;
        src_storage++;
      }
    } break;
#if !NUMBER_DOUBLE_ONLY
    case NumberType::typeFloat128:
    {
      const __float128 *src_storage=(const __float128 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=*src_storage;
        src_storage++;
      }
    } break;
    case NumberType::typeDDouble:
    case NumberType::typeDDHeap:
    {
      const dd_real *src_storage=(const dd_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->hi;
        src_storage++;
      }
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real *src_storage=(const dq_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->hi;
        src_storage++;
      }
    } break;
    case NumberType::typeReal642:
    {
      const real642 *src_storage=(const real642 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->val1;
        src_storage++;
      }
    } break;
#endif // !NUMBER_DOUBLE_ONLY
  }
  return result;
}

template<> template<typename OTHER_BASE>
number<double> &number<double>::assign_across(number<OTHER_BASE> const &src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;

  //store=src.store;
  //return *this;
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<double> &self) { self.store=0; }
    void operator()(number<double> const &other, number<double> &self) { self.store=other.store; }
    void operator()(number<__float128> const &other, number<double> &self) { self.store=other.store; }
    void operator()(number<dd_real> const &other, number<double> &self) { self.store=other.store.hi; }
    void operator()(number<dq_real> const &other, number<double> &self) { self.store=other.store.hi; }
    void operator()(number<real642> const &other, number<double> &self) { self.store=other.store.val1; }
    void operator()(number<dd_heap> const &other, number<double> &self) { self.store=other.store.value().hi; }
    void operator()(number<number_any> const &other, number<double> &self) { visit(*this, other.store, self); }
  } assigner;
  assigner(src, *this);
  return *this;
}

/*template<> template <>
void number<double>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeDouble);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  store=access->store;
}*/

/*template<>
number<double> &number<double>::assign_across(number_any const &src)
{
  //working_assert(src->ntype()==typeDouble);
  //working_assert(src->raw_ntype()==typeDouble);
  / *visit([](const auto &other, number<double> &self){
    using T=std::decay_t<decltype(other)>;
    if constexpr(std::is_same_v<T, std::monostate>)
      self.store=0;
    else
      self.assign_across(other);
  }, src, *this);* /
  class Assigner
  {
    public:
    void operator()(std::monostate const &, number<double> &self) { self.store=0; }
    void operator()(number<double> const &other, number<double> &self) { self.store=other.store; }
    void operator()(number<__float128> const &other, number<double> &self) { self.store=other.store; }
    void operator()(number<dd_real> const &other, number<double> &self) { self.store=other.store.hi; }
    void operator()(number<dq_real> const &other, number<double> &self) { self.store=other.store.hi; }
    void operator()(number<real642> const &other, number<double> &self) { self.store=other.store.val1; }
  } assigner;
  visit(assigner, src, *this);
  return *this;*/

  /*switch (src.store.index())
  {
    case NumberType::typeEmpty: assign_across(((const number<number_a *> *)src)->store); break;
    case NumberType::typeDouble: store=((const number<double> *)src)->store; break;
    case NumberType::typeFloat128: store=((const number<__float128> *)src)->store; break;
    case NumberType::typeDDouble:
    {
      const dd_real &access=((const number<dd_real> *)src)->store;
      store=access.hi;
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real &access=((const number<dq_real> *)src)->store;
      store=access.hi;
    } break;
    case NumberType::typeReal642:
    {
      const real642 &access=((const number<real642> *)src)->store;
      store=access.val1;
    } break;
  }*/
//}

QString number_traits<double>::toString(double store)
{
  char buf[32];
  *std::to_chars(buf, buf+sizeof(buf)-1, store).ptr=0; // shortest round-trip, locale-independent
  return QString(buf);
}

bool number_traits<double>::fromString(double &s, char const *str)
{
  while (std::isspace((unsigned char)*str))
    str++;
  if (*str=='+')
    str++;
  char const *end=str+std::strlen(str);
  auto [ptr, ec]=std::from_chars(str, end, s);
  if (ec!=std::errc())
    return false;
  while (std::isspace((unsigned char)*ptr))
    ptr++;
  return ptr==end;
}

int number_traits<double>::toRound(double store)
{
  return qRound(store);
}

#if !NUMBER_DOUBLE_ONLY
__float128 *number_traits<__float128>::convert_block(NumberType old_type, const void *old_data, int count)
{
  __float128 *result=new __float128[count];
  switch (old_type)
  {
    case NumberType::typeEmpty:
    {
      for (int i=0; i<count; i++)
        result[i]=0;
    } break;
#if NUMBER_DOUBLE_EXISTS
    case NumberType::typeDouble:
    {
      const double *src_storage=(const double *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=*src_storage;
        src_storage++;
      }
    } break;
#endif
    case NumberType::typeFloat128:
    {
      const __float128 *src_storage=(const __float128 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=*src_storage;
        src_storage++;
      }
    } break;
    case NumberType::typeDDouble:
    case NumberType::typeDDHeap:
    {
      const dd_real *src_storage=(const dd_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real *src_storage=(const dq_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case NumberType::typeReal642:
    {
      const real642 *src_storage=(const real642 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i]=src_storage->val1;
        src_storage++;
      }
    } break;
  }
  return result;
}

template<> template<typename OTHER_BASE>
number<__float128> &number<__float128>::assign_across(number<OTHER_BASE> const &src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;

  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<__float128> &self) { self.store=0; }
    void operator()(number<double> const &other, number<__float128> &self) { self.store=other.store; }
    void operator()(number<__float128> const &other, number<__float128> &self) { self.store=other.store; }
    void operator()(number<dd_real> const &other, number<__float128> &self) { self.store=other.store.hi+__float128(other.store.lo_); }
    void operator()(number<dq_real> const &other, number<__float128> &self) { self.store=other.store.hi+__float128(other.store.lo_); }
    void operator()(number<real642> const &other, number<__float128> &self) { self.store=other.store.val1; }
    void operator()(number<dd_heap> const &other, number<__float128> &self) { self.store=other.store.value().hi+__float128(other.store.value().lo_); }
    void operator()(number<number_any> const &other, number<__float128> &self) { visit(*this, other.store, self); }
  } assigner;
  assigner(src, *this);
  return *this;
/*
  using T=std::remove_cv_t<std::remove_reference_t<OTHER>>;
  if constexpr (std::is_same_v<T, double>)
    store=src.store;
  else if constexpr (std::is_same_v<T, __float128>)
    store=src.store;
  else if constexpr (std::is_same_v<T, dd_real>)
    //const dd_real &access=((const number<dd_real> *)src)->store;
    //store=access.hi;
    store=src.store.hi+__float128(src.store.lo_);
  else if constexpr (std::is_same_v<T, dq_real>)
    store=src.store.hi+__float128(src.store.lo_);
  else if constexpr (std::is_same_v<T, real642>)
    store=src.store.val1;
  else
    //static_assert(false, "unknown number type");
    static_assert(!sizeof(T), "unknown number type");*/
  return *this;
}

/*template<> template <>
void number<__float128>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeFloat128);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<__float128> *access=specific_cast<number<__float128> *, number_a *>(src->store);
  store=access->store;
}*/

/*template<>
number<__float128> &number<__float128>::assign_across(number_any const &src)
{
  //working_assert(src->ntype()==typeFloat128);
  //working_assert(src->raw_ntype()==typeFloat128);

  / *visit([](auto const &other, number<__float128> &self){
      using T=std::decay_t<decltype(other)>;
      if constexpr(std::is_same_v<T, std::monostate>)
          self.store=0;
      else
          self.assign_across(other);
  }, src, *this);* /
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<__float128> &self) { self.store=0; }
    void operator()(number<double> const &other, number<__float128> &self) { self.store=other.store; }
    void operator()(number<__float128> const &other, number<__float128> &self) { self.store=other.store; }
    void operator()(number<dd_real> const &other, number<__float128> &self) { self.store=other.store.hi+__float128(other.store.lo_); }
    void operator()(number<dq_real> const &other, number<__float128> &self) { self.store=other.store.hi+__float128(other.store.lo_); }
    void operator()(number<real642> const &other, number<__float128> &self) { self.store=other.store.val1; }
  } assigner;
  visit(assigner, src, *this);
  return *this;*/

  /*switch (src->raw_ntype())
  {
    case NumberType::typeEmpty: assign_across(((const number<number_a *> *)src)->store); break;
    case NumberType::typeDouble: store=((const number<double> *)src)->store; break;
    case NumberType::typeFloat128: store=((const number<__float128> *)src)->store; break;
    case NumberType::typeDDouble:
    {
      const dd_real &access=((const number<dd_real> *)src)->store;
      store=access.hi+__float128(access.lo_);
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real &access=((const number<dq_real> *)src)->store;
      store=access.hi+__float128(access.lo_);
    } break;
    case NumberType::typeReal642:
    {
      const real642 &access=((const number<real642> *)src)->store;
      store=access.val1;
    } break;
  }*/
//}

int number_traits<__float128>::toRound(__float128 store)
{
  return qRound((double)store);
}

dd_real *number_traits<dd_real>::convert_block(NumberType old_type, const void *old_data, int count)
{
  dd_real *result=new dd_real[count];
  switch (old_type)
  {
    case NumberType::typeEmpty:
    {
      for (int i=0; i<count; i++)
        result[i].zero(0);
    } break;
    case NumberType::typeDouble:
    {
      const double *src_storage=(const double *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].zero(*src_storage);
        src_storage++;
      }
    } break;
    case NumberType::typeFloat128:
    {
      const __float128 *src_storage=(const __float128 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].hi=(double)*src_storage;
        result[i].lo_=*src_storage-result[i].hi;
        src_storage++;
      }
    } break;
    case NumberType::typeDDouble:
    case NumberType::typeDDHeap:
    {
      const dd_real *src_storage=(const dd_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real *src_storage=(const dq_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case NumberType::typeReal642:
    {
      const real642 *src_storage=(const real642 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].hi=src_storage->val1;
        result[i].lo_=0;
        src_storage++;
      }
    } break;
  }
  return result;
}

/*template<> //template <>
number<dd_real> &number<dd_real>::assign_across(const number<dd_real> &src)
{
  store=src.store;
  return *this;
}*/

/*template<> template <>
void number<dd_real>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeDDouble);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<dd_real> *access=specific_cast<number<dd_real> *, number_a *>(src->store);
  store=access->store;
}*/

template<> template<typename OTHER_BASE>
number<dd_real> &number<dd_real>::assign_across(number<OTHER_BASE> const &src)
{
  //working_assert(src->ntype()==typeDDouble);
  //working_assert(src->raw_ntype()==typeDDouble);
  //store=specific_cast<const number<dd_real> *, const number_a *>(src)->store;
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<dd_real> &self) { self.store.zero(); }
    void operator()(number<double> const &other, number<dd_real> &self) { self.store.zero(other.store); }
    void operator()(number<__float128> const &other, number<dd_real> &self) { self.store.hi=other.store; self.store.lo_=other.store-self.store.hi; }
    void operator()(number<dd_real> const &other, number<dd_real> &self) { self.store=other.store; }
    void operator()(number<dq_real> const &other, number<dd_real> &self) { self.store=other.store; }
    void operator()(number<real642> const &other, number<dd_real> &self) { self.store.zero(other.store.val1); }
    void operator()(number<dd_heap> const &other, number<dd_real> &self) { self.store=other.store.value(); }
    void operator()(number<number_any> const &other, number<dd_real> &self) { visit(*this, other.store, self); }
  } assigner;
  assigner(src, *this);
  //return *this;

/*
  using T=std::remove_cv_t<std::remove_reference_t<OTHER>>;
  //could do if (monostate) store.hi=store.lo_=0;
  if constexpr (std::is_same_v<T, double>)
    { store.hi=src.store; store.lo_=0; }
  else if constexpr (std::is_same_v<T, __float128>)
  {
    store.hi=src.store;
    store.lo_=src.store-store.hi;
  }
  else if constexpr (std::is_same_v<T, dd_real>)
      //const dd_real &access=((const number<dd_real> *)src)->store;
      //store=access.hi;
      store=src.store;
  else if constexpr (std::is_same_v<T, dq_real>)
      store=src.store;
  else if constexpr (std::is_same_v<T, real642>)
  {
    store.hi=src.store.val1;
    store.lo_=0;
  }
  else
      //static_assert(false, "unknown number type");
    static_assert(!sizeof(T), "unknown number type");*/
  return *this;
}

/*template<>
number<dd_real> &number<dd_real>::assign_across(number_any const &src)
{
  //working_assert(src->ntype()==typeFloat128);
  //working_assert(src->raw_ntype()==typeFloat128);

  / *mvisit([](auto const &other, number<dd_real> &self){
      using T=std::decay_t<decltype(other)>;
      if constexpr(std::is_same_v<T, std::monostate>)
        { self.store.hi=self.store.lo_=0; }
      else
        self.assign_across(other);
  }, src.store, *this);* /
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<dd_real> &self) { self.store.zero(); }
    void operator()(number<double> const &other, number<dd_real> &self) { self.store.zero(other.store); }
    void operator()(number<__float128> const &other, number<dd_real> &self) { self.store.hi=other.store; self.store.lo_=other.store-self.store.hi; }
    void operator()(number<dd_real> const &other, number<dd_real> &self) { self.store=other.store; }
    void operator()(number<dq_real> const &other, number<dd_real> &self) { self.store=other.store; }
    void operator()(number<real642> const &other, number<dd_real> &self) { self.store.zero(other.store.val1); }
  } assigner;
  visit(assigner, src, *this);
  return *this;
}*/

dq_real *number_traits<dq_real>::convert_block(NumberType old_type, const void *old_data, int count)
{
  dq_real *result=new dq_real[count];
  switch (old_type)
  {
  case NumberType::typeEmpty:
  {
    for (int i=0; i<count; i++)
        result[i].zero(0);
  } break;
  case NumberType::typeDouble:
  {
    const double *src_storage=(const double *)old_data;
    for (int i=0; i<count; i++)
    {
        result[i].zero(*src_storage);
        src_storage++;
    }
  } break;
  case NumberType::typeFloat128:
  {
    const __float128 *src_storage=(const __float128 *)old_data;
    for (int i=0; i<count; i++)
    {
        result[i].hi=(double)*src_storage;
        result[i].lo_=*src_storage-result[i].hi;
        src_storage++;
    }
  } break;
  case NumberType::typeDDouble:
  case NumberType::typeDDHeap:
  {
    const dd_real *src_storage=(const dd_real *)old_data;
    for (int i=0; i<count; i++)
    {
        result[i].assign(*src_storage);
        src_storage++;
    }
  } break;
  case NumberType::typeQDouble:
  {
    const dq_real *src_storage=(const dq_real *)old_data;
    for (int i=0; i<count; i++)
    {
        result[i].assign(*src_storage);
        src_storage++;
    }
  } break;
  case NumberType::typeReal642:
  {
    const real642 *src_storage=(const real642 *)old_data;
    for (int i=0; i<count; i++)
    {
        result[i].hi=src_storage->val1;
        result[i].lo_=0;
        src_storage++;
    }
  } break;
  }
  return result;
}

/*template<> //template <>
number<dq_real> &number<dq_real>::assign_across(const number<dq_real> &src)
{
  store=src.store;
  return *this;
}*/

/*template<> template <>
void number<dq_real>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeDDouble);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<dq_real> *access=specific_cast<number<dq_real> *, number_a *>(src->store);
  store=access->store;
}*/

template<> template<typename OTHER_BASE>
number<dq_real> &number<dq_real>::assign_across(number<OTHER_BASE> const &src)
{
  //working_assert(src->ntype()==typeDDouble);
  //working_assert(src->raw_ntype()==typeDDouble);
  //store=specific_cast<const number<dq_real> *, const number_a *>(src)->store;
  class Assigner
  {
public:
    void operator()(std::monostate const &, number<dq_real> &self) { self.store.zero(); }
    void operator()(number<double> const &other, number<dq_real> &self) { self.store.zero(other.store); }
    void operator()(number<__float128> const &other, number<dq_real> &self) { self.store.hi=other.store; self.store.lo_=other.store-self.store.hi; }
    void operator()(number<dd_real> const &other, number<dq_real> &self) { self.store.hi=other.store.hi; self.store.lo_=other.store.lo_; }
    void operator()(number<dq_real> const &other, number<dq_real> &self) { self.store=other.store; }
    void operator()(number<real642> const &other, number<dq_real> &self) { self.store.zero(other.store.val1); }
    void operator()(number<dd_heap> const &other, number<dq_real> &self) { self.store.hi=other.store.value().hi; self.store.lo_=other.store.value().lo_; }
    void operator()(number<number_any> const &other, number<dq_real> &self) { visit(*this, other.store, self); }
  } assigner;
  assigner(src, *this);

  /*using T=std::remove_cv_t<std::remove_reference_t<OTHER>>;
  //could do if (monostate) store.hi=store.lo_=0;
  if constexpr (std::is_same_v<T, double>)
  { store.hi=src.store; store.lo_=0; }
  else if constexpr (std::is_same_v<T, __float128>)
  {
    store.hi=src.store;
    store.lo_=src.store-store.hi;
  }
  else if constexpr (std::is_same_v<T, dd_real>)
  {
    //const dd_real &access=((const number<dd_real> *)src)->store;
    //store=access.hi;
    store.hi=src.store.hi;
    store.lo_=src.store.lo_;
  }
  else if constexpr (std::is_same_v<T, dq_real>)
    store=src.store;
  else if constexpr (std::is_same_v<T, real642>)
  {
    store.hi=src.store.val1;
    store.lo_=0;
  }
  else
      //static_assert(false, "unknown number type");
    static_assert(!sizeof(T), "unknown number type");*/
  return *this;
}

/*template<>
number<dq_real> &number<dq_real>::assign_across(number_any const &src)
{
  //working_assert(src->ntype()==typeFloat128);
  //working_assert(src->raw_ntype()==typeFloat128);

  / *mvisit([](auto const &other, number<dq_real> &self){
      using T=std::decay_t<decltype(other)>;
      if constexpr(std::is_same_v<T, std::monostate>)
      { self.store.hi=self.store.lo_=0; }
      else
        self.assign_across(other);
  }, src.store, *this);* /
  class Assigner
  {
public:
    void operator()(std::monostate const &, number<dq_real> &self) { self.store.zero(); }
    void operator()(number<double> const &other, number<dq_real> &self) { self.store.zero(other.store); }
    void operator()(number<__float128> const &other, number<dq_real> &self) { self.store.hi=other.store; self.store.lo_=other.store-self.store.hi; }
    void operator()(number<dd_real> const &other, number<dq_real> &self) { self.store.hi=other.store.hi; self.store.lo_=other.store.lo_; }
    void operator()(number<dq_real> const &other, number<dq_real> &self) { self.store=other.store; }
    void operator()(number<real642> const &other, number<dq_real> &self) { self.store.zero(other.store.val1); }
  } assigner;
  visit(assigner, src, *this);
  return *this;
}*/









template<> template<typename OTHER_BASE>
number<dd_heap> &number<dd_heap>::assign_across(number<OTHER_BASE> const &src)
{
  number<dd_real> value;
  value.assign_across(src);
  store=value.store;
  return *this;
}

real642 *number_traits<real642>::convert_block(NumberType old_type, const void *old_data, int count)
{
  real642 *result=new real642[count];
  switch (old_type)
  {
    case NumberType::typeEmpty:
    {
      for (int i=0; i<count; i++)
      {
        result[i].val1=0;
        result[i].val2=0;
      }
    } break;
    case NumberType::typeDouble:
    {
      const double *src_storage=(const double *)old_data;
      //TODO: memmove but we'll see later
      for (int i=0; i<count; i++)
      {
        result[i].val1=*src_storage;
        result[i].val2=*src_storage;
        src_storage++;
      }
    } break;
    case NumberType::typeFloat128:
    {
      const __float128 *src_storage=(const __float128 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].val1=*src_storage;
        result[i].val2=*src_storage;
        src_storage++;
      }
    } break;
    case NumberType::typeDDouble:
    case NumberType::typeDDHeap:
    {
      const dd_real *src_storage=(const dd_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].val1=src_storage->hi;
        result[i].val2=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real *src_storage=(const dq_real *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].val1=src_storage->hi;
        result[i].val2=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case NumberType::typeReal642:
    {
      const real642 *src_storage=(const real642 *)old_data;
      for (int i=0; i<count; i++)
      {
        result[i].val1=src_storage->val1;
        result[i].val2=src_storage->val2;
        src_storage++;
      }
    } break;
  }
  return result;
}

/*template<> template <typename OTHER>
number<real642> &number<real642>::assign_across(number<OTHER> &src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;
  store=src->store;
}*/

/*template<> template <>
void number<real642>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeReal642);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<real642> *access=specific_cast<number<real642> *, number_a *>(src->store);
  store=access->store;
}*/

template<> template<typename OTHER_BASE>
number<real642> &number<real642>::assign_across(number<OTHER_BASE> const &src)
{
  //working_assert(src->ntype()==typeReal642);
  //working_assert(src->raw_ntype()==typeReal642);
  //store=specific_cast<const number<real642> *, const number_a *>(src)->store;
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<real642> &self) { self.store.val2=self.store.val1=0; }
    void operator()(number<double> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store; }
    void operator()(number<__float128> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store; }
    void operator()(number<dd_real> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store.hi; }
    void operator()(number<dq_real> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store.hi; }
    void operator()(number<real642> const &other, number<real642> &self) { self.store=other.store; }
    void operator()(number<dd_heap> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store.value().hi; }
    void operator()(number<number_any> const &other, number<real642> &self) { visit(*this, other.store, self); }
  } assigner;
  assigner(src, *this);

  /*using T=std::remove_cv_t<std::remove_reference_t<OTHER>>;
  if constexpr (std::is_same_v<T, double>)
    store.val2=store.val1=src.store;
  else if constexpr (std::is_same_v<T, __float128>)
    store.val2=store.val1=src.store;
  else if constexpr (std::is_same_v<T, dd_real>)
      //const dd_real &access=((const number<dd_real> *)src)->store;
      //store=access.hi;
    store.val2=store.val1=src.store.hi;
  else if constexpr (std::is_same_v<T, dq_real>)
    store.val2=store.val1=src.store.hi;
  else if constexpr (std::is_same_v<T, real642>)
    store=src.store;
  else
    //static_assert(false, "unknown number type");
    static_assert(!sizeof(T), "unknown number type");*/
  return *this;
}

/*template<>
number<real642> &number<real642>::assign_across(number_any const &src)
{
  //working_assert(src->ntype()==typeFloat128);
  //working_assert(src->raw_ntype()==typeFloat128);

  / *mvisit([](auto const &other, number<real642> &self){
      using T=std::decay_t<decltype(other)>;
      if constexpr(std::is_same_v<T, std::monostate>)
      { self.store.val2=self.store.val1=0; }
      else
        self.assign_across(other);
  }, src.store, *this);* /
  class Assigner
  {
  public:
    void operator()(std::monostate const &, number<real642> &self) { self.store.val2=self.store.val1=0; }
    void operator()(number<double> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store; }
    void operator()(number<__float128> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store; }
    void operator()(number<dd_real> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store.hi; }
    void operator()(number<dq_real> const &other, number<real642> &self) { self.store.val2=self.store.val1=other.store.hi; }
    void operator()(number<real642> const &other, number<real642> &self) { self.store=other.store; }
  } assigner;
  visit(assigner, src, *this);
  return *this;
}*/

int number_traits<real642>::toRound(real642 const &store)
{
  return qRound(store.val1);
}

#endif // !NUMBER_DOUBLE_ONLY

double number<number_any>::eps2() const
{
  return visit([](auto const &val){
      using T=std::remove_cv_t<std::remove_reference_t<decltype(val)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return 1.0;
      else
          return val.eps2();
  }, store);
}

double number<number_any>::eps234() const
{
  return visit([](auto &val){
      using T=std::remove_cv_t<std::remove_reference_t<decltype(val)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return 2.0;
      else
          return val.eps234();
  }, store);
}

double number<number_any>::sqrtlarge() const
{
  return visit([](auto &val){
      using T=std::remove_cv_t<std::remove_reference_t<decltype(val)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return 1.915e153;
      else
          return val.sqrtlarge();
  }, store);
}

NumberType number<number_any>::ntype() const
{
  return visit([](auto &val){
      //using T = std::decay_t<decltype(val)>;
      using T = std::remove_cv_t<std::remove_reference_t<decltype(val)>>;
      //using T = decltype(val);
      if constexpr (std::is_same_v<T, std::monostate>)
          return NumberType::typeEmpty;
      else
          return val.ntype();
  }, store);
  //return store->ntype();
}

NumberType number<number_any>::raw_ntype() const
{
  return NumberType::typeEmpty;
}

/*template<>
number<number_any>::number(number<number_any> const &src) noexcept: store(src.store)
{

}*/

number<number_any>::number(NumberType ntype): store()
{
  static_assert(std::is_copy_assignable_v<std::monostate>, "woot");
  static_assert(std::is_copy_assignable_v<number<double>>, "woot");
  static_assert(std::is_copy_assignable_v<number<__float128>>, "woot");
  static_assert(std::is_copy_assignable_v<number<dd_real>>, "woot");
  static_assert(std::is_copy_assignable_v<number<dq_real>>, "woot");
  static_assert(std::is_copy_assignable_v<number<real642>>, "woot");
  static_assert(std::is_copy_assignable_v<number_any>, "woot");
  switch (ntype)
  {
    //case NumberType::typeEmpty: dbgPoint(); store=new number<double>(); break;
    case NumberType::typeDouble: store.emplace<number<double>>(); break;
  //case NumberType::typeDouble: store.emplace<1>(wtf); break;
#if !NUMBER_DOUBLE_ONLY
    case NumberType::typeFloat128: store.emplace<number<__float128>>(); break;
    case NumberType::typeDDouble: store.emplace<number<dd_real>>(); break;
    case NumberType::typeQDouble: store.emplace<number<dq_real>>(); break;
    case NumberType::typeReal642: store.emplace<number<real642>>(); break;
    case NumberType::typeDDHeap: store.emplace<number<dd_heap>>(); break;
#endif
    default:
      dbgPoint();
      working_assert(false); // ntype is typeEmpty — not properly initialised
      store.emplace<number<double>>();
      break;
  }
}

void number<number_any>::constructLateBecauseQtIsAwesome(const number<number_any> &src)
{
  working_assert(store.index()==0);
  working_assert(src.store.index()==0); //was ==0 but doesn't make a lot of sense
  store=src.store;
}

number<number_any>::~number()
{
  nop();
  //delete store;
  //store=nullptr;
}

void number<number_any>::readFrom(void *storage, int index)
{
  visit([](auto &value, void *storage, int index) {
      using T=std::decay_t<decltype(value)>;
      if constexpr (std::is_same_v<T, std::monostate>)
        { }
      else
        value.readFrom(storage, index);
      }, store, storage, index);
  //store->readFrom(storage, index);
}

void number<number_any>::writeTo(void *storage, int index) const
{
  visit([&storage, &index](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
          value.writeTo(storage, index);
  }, store);
  //store->writeTo(storage, index);
}

number<number_any> &number<number_any>::zero(double val)
{
  visit([val](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
          value.zero(val);
  }, store);
  //store->zero(val);
  return *this;
}

number<number_any> &number<number_any>::assign(const number<number_any> &src) noexcept
{
  //working_assert(src.raw_ntype()==typeEmpty);
  working_assert(store.index()!=0);
  working_assert(src.store.index()!=0);
  working_assert(store.index()==src.store.index());
  store=src.store;
  return *this;
}

/*template<>
void number<number_a *>::assign(const number_a &src)
{
  (void)src;
  dbgPoint();
}*/

/*template<> template <>
void number<number_a *>::assign_across<double>(const number<double> *src)
{
  working_assert(store->ntype()==typeDouble);
  working_assert(store->raw_ntype()==typeDouble);
  specific_cast<number<double> *, number_a *>(store)->assign(*src);
}

template<> template <>
void number<number_a *>::assign_across<__float128>(const number<__float128> *src)
{
  specific_cast<number<__float128> *, number_a *>(store)->assign(*src);
}

template<> template <>
void number<number_a *>::assign_across<dd_real>(const number<dd_real> *src)
{
  specific_cast<number<dd_real> *, number_a *>(store)->assign(*src);
}

template<> template <>
void number<number_a *>::assign_across<real642>(const number<real642> *src)
{
  specific_cast<number<real642> *, number_a *>(store)->assign(*src);
}*/

template<typename OTHER_BASE>
number<number_any> &number<number_any>::assign_across(number<OTHER_BASE> const &src)
{
  //working_assert(src->raw_ntype()==typeEmpty);
  //working_assert(store->ntype()==src->store->ntype());
  //working_assert(store->raw_ntype()==src->store->raw_ntype());
  //TODO: convert if store->ntype()!=src->store->ntype()
  class Assigner
  {
  public:
    void operator()(std::monostate &, number<OTHER_BASE> const &) { working_assert(false); }
    void operator()(number<double> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<__float128> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<dd_real> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<dq_real> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<real642> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<dd_heap> &my_store, number<OTHER_BASE> const &other) { my_store.assign_across(other); }
    void operator()(number<number_any> &my_store, number<OTHER_BASE> const &other) {
    //  visit([](auto &my_store, auto const &other_store) { my_store.assign_across(other_store); }, self.store, other.store);
      (void)my_store;
      (void)other;
      working_assert(false);
    }
  } assigner;
  visit(assigner, store, src);

  /*visit([](auto &value, number_any const &src) {
      using T=std::decay_t<decltype(value)>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
        value.assign_across(src);
        / *  visit([](auto const &other_value, decltype(value) &value2) {
            using O=std::decay_t<decltype(other_value)>;
            if constexpr (std::is_same_v<O, std::monostate>)
              { }
            else if constexpr (std::is_same_v<O, T>)
            {
              value2.assign_across(other_value);
            }
            else
              working_assert(false);
          }, src.store, value);* /
      }
    }, store, src);*/
  return *this;
}

/*template<>
void number<number_a *>::assign_across(const number_a *src)
{
  //working_assert(store->ntype()==src->store->ntype());
  //store->assign(*src->store);
  (void)src;
  dbgPoint();
}*/

number<number_any> &number<number_any>::chs()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.chs();
  }, store);
  return *this;
}

number<number_any> &number<number_any>::lshift(int shoft)
{
  visit([shoft](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.lshift(shoft);
  }, store);
  //store->zero(val);
  return *this;
}

number<number_any> &number<number_any>::round()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.round();
  }, store);
  return *this;
}

number<number_any> &number<number_any>::frac()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.frac();
  }, store);
  return *this;
}

number<number_any> &number<number_any>::mod1()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.mod1();
  }, store);
  return *this;
}

number<number_any> &number<number_any>::add_double(double x)
{
  visit([x](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.add_double(x);
  }, store);
  return *this;
}

number<number_any> &number<number_any>::mul_double(double x)
{
  visit([x](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.mul_double(x);
  }, store);
  return *this;
}

number<number_any> &number<number_any>::add(const number<number_any> &other)
{
  visit([&other](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
          visit([&value](auto const &other_value) {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
              { }
              else if constexpr (std::is_same_v<O, T>)
                  value.add(other_value);
              else
              {
                  working_assert(false);
              }
          }, other.store);
      }
  }, store);
  return *this;
}

/*template<>
number<number_any> &number<number_any>::add(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}*/

number<number_any> &number<number_any>::sub(const number<number_any> &other)
{
  visit([&other](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
          visit([&value](auto &other_value) {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
              { }
              else if constexpr (std::is_same_v<O, T>)
                  value.sub(other_value);
              else
              {
                  working_assert(false);
              }
          }, other.store);
      }
  }, store);
  return *this;
}

/*template<>
void number<number_a *>::sub(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}*/

number<number_any> &number<number_any>::rsub(const number<number_any> &other)
{
  visit([&other](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
          visit([&value](auto &other_value) {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
              { }
              else if constexpr (std::is_same_v<O, T>)
              {
                  value.rsub(other_value);
              }
              else
              {
                  working_assert(false);
              }
          }, other.store);
      }
  }, store);
  return *this;
}

/*template<>
void number<number_a *>::rsub(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}*/

number<number_any> &number<number_any>::mul(const number<number_any> &other)
{
  visit([&other](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
          visit([&value](auto &other_value) {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
              { }
              else if constexpr (std::is_same_v<O, T>)
              {
                 value.mul(other_value);
              }
              else
              {
                  working_assert(false);
              }
          }, other.store);
      }
  }, store);
  return *this;
}

/*template<>
void number<number_a *>::mul(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}*/

number<number_any> &number<number_any>::sqr()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
          value.sqr();
  }, store);
  return *this;
}

double number<number_any>::radixfloor() const
{
  return visit([](auto &value)->double {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return 1;
      else
      {
        return value.radixfloor();
      }
  }, store);
}

/*template<>
double number<number_a *>::radixfloor(const number_a &store2) const
{
  (void)store2;
  dbgPoint();
  return 1;
}*/

number<number_any> &number<number_any>::recip()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.recip();
  }, store);
  return *this;
}

number<number_any> &number<number_any>::sqrt()
{
  visit([](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.sqrt();
  }, store);
  return *this;
}

bool number<number_any>::reduce_angle()
{
  return visit([](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return false;
      else
        return value.reduce_angle();
  }, store);
}

number<number_any> &number<number_any>::add_pi(double x)
{
  visit([](auto &value, double x) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
        value.add_pi(x);
  }, store, x);
  return *this;
}

std::strong_ordering number<number_any>::compare(const number<number_any> &other) const
{
  return visit([](auto const &value, const number<number_any> &other)->std::strong_ordering {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return std::strong_ordering::equal;
      else
      {
          return visit([](auto const &other_value, auto const &value)->std::strong_ordering {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
                  return std::strong_ordering::equal;
              else if constexpr (std::is_same_v<O, T>)
                  return value.compare(other_value);
              else
              {
                  working_assert(false);
                  return std::strong_ordering::equal;
              }
          }, other.store, value);
      }
  }, store, other);
}

/*template<>
int number<number_a *>::compare(const number_a &other) const
{
  (void)other;
  dbgPoint();
  return 0;
}*/

bool number<number_any>::isequal(const number<number_any> &other) const
{
  return visit([&other](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return true;
      else
      {
          return visit([&value](auto &other_value)->bool {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
                  return true;
              else if constexpr (std::is_same_v<O, T>)
                  return value.isequal(other_value);
              else
              {
                  working_assert(false);
                  return true;
              }
          }, other.store);
      }
  }, store);
}

/*template<>
bool number<number_a *>::isequal(const number_a &other) const
{
  (void)other;
  dbgPoint();
  return false;
}*/

bool number<number_any>::is0() const
{
  return visit([](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return true;
      else
        return value.is0();
  }, store);
}

bool number<number_any>::isle(const number<number_any> &other) const
{
  return visit([&other](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
          return true;
      else
      {
          return visit([&value](auto &other_value)->bool {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
                  return true;
              else if constexpr (std::is_same_v<O, T>)
                  return value.isle(other_value);
              else
              {
                  working_assert(false);
                  return true;
              }
          }, other.store);
      }
  }, store);
}

/*template<>
bool number<number_a *>::isle(const number_a &other) const
{
  (void)other;
  dbgPoint();
  return true;
}*/

bool number<number_any>::isle0() const
{
  return visit([](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return true;
      else
        return value.isle0();
  }, store);
}

bool number<number_any>::isl0() const
{
  return visit([](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return false;
      else
        return value.isl0();
  }, store);
}

bool number<number_any>::isl1() const
{
  return visit([](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return true;
      else
        return value.isl1();
  }, store);
}

number<number_any> &number<number_any>::min(const number<number_any> &other)
{
  visit([&other](auto &value) {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
      { }
      else
      {
          visit([&value](auto &other_value) {
              using O=std::remove_cv_t<std::remove_reference_t<decltype(other_value)>>;
              if constexpr (std::is_same_v<O, std::monostate>)
              { }
              else if constexpr (std::is_same_v<O, T>)
              {
                  value.min(other_value);
              }
              else
              {
                  working_assert(false);
              }
          }, other.store);
      }
  }, store);
  return *this;
}

/*template<>
void number<number_a *>::min(const number_a &other)
{
  (void)other;
  dbgPoint();
}*/

QString number<number_any>::toString() const
{
  return visit([](auto &value)->QString {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return QString("none");
      else
        return value.toString();
  }, store);
}

bool number<number_any>::fromString(QString const &str)
{
  return visit([&str](auto &value)->bool {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return false;
      else
        return value.fromString(str);
  }, store);
}

int number<number_any>::toRound() const
{
  return visit([](auto &value)->int {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return 0;
      else
        return value.toRound();
  }, store);
}

double number<number_any>::toDouble() const
{
  return visit([](auto &value)->double {
      using T=std::remove_cv_t<std::remove_reference_t<decltype(value)>>;
      if constexpr (std::is_same_v<T, std::monostate>)
        return 0;
      else
        return value.toDouble();
  }, store);
}

template<typename BASE> template<typename OTHER_BASE>
complex<BASE> &complex<BASE>::assign_across(const complex<OTHER_BASE> &src)
{
  re.assign_across(src.re);
  im.assign_across(src.im);
  return *this;
}

template<typename BASE>
QString complex<BASE>::toString() const
{
  return re.toString()+" +i* "+im.toString();
}

void complex_double_sqrt(double *res_re, double *res_im, double in_re, double in_im)
{
  double t1;
  if (in_re>=0)
  {
    t1=sqrt((sqrt(in_re*in_re+in_im*in_im)+in_re)/2);
    *res_re=t1;
    if (t1==0)
      *res_im=0;
    else
      *res_im=in_im/(2*t1);
  }
  else
  {
    t1=sqrt((sqrt(in_re*in_re+in_im*in_im)-in_re)/2);
    if (in_im>=0)
    {
      *res_re=in_im/(2*t1);
      *res_im=t1;
    }
    else
    {
      *res_re=-in_im/(2*t1);
      *res_im=-t1;
    };
  };
}

double radixfloor_double(double x1, double x2)
{
  int ilog1=std::ilogb(x1);
  int ilog2=std::ilogb(x2);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

/*
 finds smaller (in abs) root of ax^2+2bx+c=0
*/
void real_double_quadratic(double *res,
                           double a, double b2, double c)
{
  /*
  a x^2 + 2 b x + c = 0
  x1,x2= -(b +- sqrt(b^2-a*c))/a
  for small a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/(a*c/b^2)*c/b^2
  HELP(x)=(1+-sqrt(1-x))/x =1/(1-+sqrt(1-x))
  x1,x2= -c/b*HELP(a*c/b^2)   or rather, for a<b, c<b
  for b<a
  (F1)  x1,x2= -(b/a +- sqrt(b^2/a^2-c/a))     good for a>b, c>b until ~ c>b^2/a
  (F1') x1,x2= -(1/a)*(b - sqrt(b^2-a*c))      good for b^2<<a*c
  (F2)  x1,x2= -(b/a)*(1 +- sqrt(1-a*c/b^2))   good for a>b, c<b^2/a
  for b>a  b^2/a>b
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a/c*b^2*c/b^2
        x1,x2= -(c/b)*HELP(a*c/b^2)            good for c<b, until c<b^2/a
  (F3)  x1=    -(c/b)/(1+sqrt(1-a*c/b^2))      good for b>0 (cancellation for x2)
  (F3')        -(c)/(b+-sqrt(b^2-a*c))         good except both b,c small e.g. 0

  Muller's method: B=2b   2c/(-B+-sqrt(B^2-4ac))
  x1,x2= c/(-b+-sqrt(b^2-ac))     -c/(b+sqrt(b^2-ac))
  x1,x2= c/b/(-1+-sqrt(1-ac/b^2))
  */
  double bb=b2*b2;
  double bbac=bb-a*c; //b^2-ac
  double d=std::sqrt(bbac); //caller should make sure that c<0 so that it's rather b^2+a*c
  double t1;
  if (b2*d<0)//if (b2<0)
  {
    t1=b2-d;
  }
  else
  {
    t1=b2+d;
  }
  if (std::abs(t1)>std::abs(a)) //do we prefer to divide by t1 or by a? the bigger!
  { //F3'
    *res=-c/t1; //-c/(b+sqrt(b^2-a*c))
  }
  else
  { //when a is large
    t1=2*b2-t1;//b2-d;
    *res=-t1/a; //-(b - sqrt(b^2-a*c))/a
  }
}

/*
 finds smaller (in abs) root of ax^2+2bx+c=0
*/
void complex_double_quadratic(double *res_re, double *res_im,
                              double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im)
{
  /*
  a x^2 + 2 b x + c = 0
  x1,x2= -(b +- sqrt(b^2-a*c))/a
  for small a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/(a*c/b^2)*c/b^2
  HELP(x)=(1+-sqrt(1-x))/x =1/(1-+sqrt(1-x))
  x1,x2= -c/b*HELP(a*c/b^2)   or rather, for a<b, c<b
  for b<a
  (F1)  x1,x2= -(b/a +- sqrt(b^2/a^2-c/a))     good for a>b, c>b until ~ c>b^2/a
  (F1') x1,x2= -(1/a)*(b - sqrt(b^2-a*c))      good for b^2<<a*c
  (F2)  x1,x2= -(b/a)*(1 +- sqrt(1-a*c/b^2))   good for a>b, c<b^2/a
  for b>a  b^2/a>b
  x1,x2= -b*(1 +- sqrt(1-a*c/b^2))/a/c*b^2*c/b^2
        x1,x2= -(c/b)*HELP(a*c/b^2)            good for c<b, until c<b^2/a
  (F3)  x1=    -(c/b)/(1+sqrt(1-a*c/b^2))      good for b>0 (cancellation for x2)
  (F3')        -(c)/(b+-sqrt(b^2-a*c))         good except both b,c small e.g. 0

  Muller's method: B=2b   2c/(-B+-sqrt(B^2-4ac))
  x1,x2= c/(-b+-sqrt(b^2-ac))     -c/(b+sqrt(b^2-ac))
  x1,x2= c/b/(-1+-sqrt(1-ac/b^2))
  */
  double bb_re=b2_re*b2_re-b2_im*b2_im;
  double bb_im=2*b2_re*b2_im;
  double bbac_re=bb_re-a_re*c_re+a_im*c_im;
  double bbac_im=bb_im-a_im*c_re-a_re*c_im; //b^2-ac
  double d_re, d_im;
  complex_double_sqrt(&d_re, &d_im, bbac_re, bbac_im);
  double t1_re, t1_im;
  if (b2_re*d_re+b2_im*d_im<0)//if (b2_re<0)
  {
    t1_re=b2_re-d_re;
    t1_im=b2_im-d_im;
  }
  else
  {
    t1_re=b2_re+d_re;
    t1_im=b2_im+d_im;
  }
  double am=a_re*a_re+a_im*a_im;
  double t1m=t1_re*t1_re+t1_im*t1_im;
  if (t1m>am)
  { //F3'
    *res_re=-(c_re*t1_re+c_im*t1_im)/t1m;
    *res_im=-(c_im*t1_re-c_re*t1_im)/t1m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  { //when? wlog. a=1: b+sqrt(b^2-c)<1, full b/a+sqrt(b^2/a^2-c/a)<1
    t1_re=2*b2_re-t1_re;//b2_re-d_re;
    t1_im=2*b2_im-t1_im;//b2_im-d_im;
    *res_re=-(t1_re*a_re+t1_im*a_im)/am;
    *res_im=-(t1_im*a_re-t1_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
}

/*
 finds both roots of ax^2+2bx+c=0
*/
void complex_double_quadratic2(double *res1_re, double *res1_im,
                               double *res2_re, double *res2_im,
                               double a_re, double a_im, double b2_re, double b2_im, double c_re, double c_im)
{
  //see complex_double_quadratic()
  double bb_re=b2_re*b2_re-b2_im*b2_im;
  double bb_im=2*b2_re*b2_im;
  double bbac_re=bb_re-a_re*c_re+a_im*c_im;
  double bbac_im=bb_im-a_im*c_re-a_re*c_im; //b^2-ac
  double d_re, d_im;
  complex_double_sqrt(&d_re, &d_im, bbac_re, bbac_im);
  double t1_re, t1_im;
  double t2_re, t2_im;
  if (b2_re*d_re+b2_im*d_im<0)//if (b2_re<0)
  {
    t1_re=b2_re-d_re;
    t1_im=b2_im-d_im;
    t2_re=b2_re+d_re;
    t2_im=b2_im+d_im;
  }
  else
  {
    t1_re=b2_re+d_re;
    t1_im=b2_im+d_im;
    t2_re=b2_re-d_re;
    t2_im=b2_im-d_im;
  } //t1=big, t2=small
  double am=a_re*a_re+a_im*a_im;
  double t1m=t1_re*t1_re+t1_im*t1_im;
  double t2m=t2_re*t2_re+t2_im*t2_im;
  if (t1m>am) //the smaller root, one way or the other
  { //F3'
    *res1_re=-(c_re*t1_re+c_im*t1_im)/t1m;
    *res1_im=-(c_im*t1_re-c_re*t1_im)/t1m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  {
    *res1_re=-(t2_re*a_re+t2_im*a_im)/am;
    *res1_im=-(t2_im*a_re-t2_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
  if (t2m>am) //the larger root, one way or the other
  { //F3'
    *res2_re=-(c_re*t2_re+c_im*t2_im)/t2m;
    *res2_im=-(c_im*t2_re-c_re*t2_im)/t2m; //-c/(b+sqrt(b^2-a*c))
  }
  else
  {
    *res2_re=-(t1_re*a_re+t1_im*a_im)/am;
    *res2_im=-(t1_im*a_re-t1_re*a_im)/am; //-(b - sqrt(b^2-a*c))/a
  }
}

// Explicit instantiations for number<T>::assign_across<number_any> (needed by MandelEvaluator, MandelModel, etc.)
template number<double> &number<double>::assign_across<number_any>(const number<number_any> &src);
template number<number_any> &number<number_any>::assign_across<number_any>(const number<number_any> &src);
#if !NUMBER_DOUBLE_ONLY
template number<__float128> &number<__float128>::assign_across<number_any>(const number<number_any> &src);
template number<dd_real> &number<dd_real>::assign_across<number_any>(const number<number_any> &src);
template number<dq_real> &number<dq_real>::assign_across<number_any>(const number<number_any> &src);
template number<real642> &number<real642>::assign_across<number_any>(const number<number_any> &src);
template number<dd_heap> &number<dd_heap>::assign_across<number_any>(const number<number_any> &src);
#endif
template class complex<double>;
template class complex<number_any>;
template complex<double> &complex<double>::assign_across<number_any>(const complex<number_any> &src);
template complex<double> &complex<double>::assign_across<double>(const complex<double> &src);
template complex<number_any> &complex<number_any>::assign_across<number_any>(const complex<number_any> &src);
#if !NUMBER_DOUBLE_ONLY
template class complex<__float128>;
template complex<__float128> &complex<__float128>::assign_across<number_any>(const complex<number_any> &src);
template complex<__float128> &complex<__float128>::assign_across<__float128>(const complex<__float128> &src);
template class complex<dd_real>;
template complex<dd_real> &complex<dd_real>::assign_across<number_any>(const complex<number_any> &src);
template complex<dd_real> &complex<dd_real>::assign_across<dd_real>(const complex<dd_real> &src);
template class complex<dq_real>;
template complex<dq_real> &complex<dq_real>::assign_across<number_any>(const complex<number_any> &src);
template complex<dq_real> &complex<dq_real>::assign_across<dq_real>(const complex<dq_real> &src);
template class complex<real642>;
template complex<real642> &complex<real642>::assign_across<number_any>(const complex<number_any> &src);
template complex<real642> &complex<real642>::assign_across<real642>(const complex<real642> &src);
template class complex<dd_heap>;
template complex<dd_heap> &complex<dd_heap>::assign_across<number_any>(const complex<number_any> &src);
template complex<dd_heap> &complex<dd_heap>::assign_across<dd_heap>(const complex<dd_heap> &src);

#endif
/*template class complex<worker_multi_float128>;
template class complex<worker_multi_ddouble>;
template class complex<worker_multi_qdouble>;
template class complex<worker_multi_real642>;*/

} // namespace MandelMath
