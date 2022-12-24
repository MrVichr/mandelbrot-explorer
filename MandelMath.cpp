//#define working_assert(x) { if (!(x)) dbgPoint(); }
#include "MandelMath.hpp"
#include <math.h>

//#include <cassert>
#include <cmath>

void doNothing(int &x)
{
  x++;
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

number_a::number_a()
{
  CreatedInstancesOfNumber++;
}

number_a::~number_a()
{
  CreatedInstancesOfNumber--;
  if (CreatedInstancesOfNumber<=0)
    nop(); //should be hit exactly once on program end
}

//these 2 have to be in front... maybe move all <number_a *> ?
template<>
NumberType number<number_a *>::ntype() const
{
  return store->ntype();
}

template<>
NumberType number<number_a *>::raw_ntype() const
{
  return NumberType::typeEmpty;
}


template<>
double number<double>::eps2() const
{
  return 1.232596e-32; // 2^-(2*53) rounded up
}

template<>
double number<double>::eps234() const
{
  return 1.17e-24; // eps2^(3/4)
}

template<>
number<double>::number(NumberType ntype): store(0)
{
  working_assert(ntype==NumberType::typeDouble);
}

template<>
number<double>::~number()
{
  store=12.345;
}

template<>
NumberType number<double>::ntype() const
{
  return NumberType::typeDouble;
}

template<>
NumberType number<double>::raw_ntype() const
{
  return NumberType::typeDouble;
}

template<>
void number<double>::readFrom(void *storage, int index)
{
  store=((double *)storage)[index];
}

template<>
void number<double>::writeTo(void *storage, int index) const
{
  ((double *)storage)[index]=store;
}

template<>
double *number<double>::convert_block(NumberType old_type, const void *old_data, int count)
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


template<>
void number<double>::zero(double val)
{
  store=val;
}

#if 1 //1..all explicit, 0..try number_a* + others
template<>
void number<double>::assign(const number<double> &src)
{
  store=src.store;
}

template<>
void number<double>::assign(const number_a &src)
{
  working_assert(src.raw_ntype()==NumberType::typeDouble);
  //const number<double> *src2=(const number<double> *)&src;
  //assign(*src2);
  assign(*(const number<double> *)&src);
}

#else
template<typename BASE>
void number<BASE>::assign(const number<BASE> &src)
{
  store=src.store;
}

#endif
template<> //template <>
void number<double>::assign_across(const number<double> *src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;
  store=src->store;
}

/*template<> template <>
void number<double>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeDouble);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  store=access->store;
}*/

template<>
void number<double>::assign_across(const number_a *src)
{
  //working_assert(src->ntype()==typeDouble);
  //working_assert(src->raw_ntype()==typeDouble);
  switch (src->raw_ntype())
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
  }
}

template<>
void number<double>::chs()
{
  store=-store;
}

template<>
void number<double>::lshift(int shoft)
{
  store=ldexp(store, shoft);
}

template<>
void number<double>::round()
{
  store=std::round(store);
}

template<>
void number<double>::frac()
{
  if (store<0)
    store -= ceil(store);
  else
    store -= floor(store);
}

template<>
void number<double>::mod1()
{
  store -= floor(store);
}

template<>
void number<double>::add_double(double x)
{
  store += x;
}

template<>
void number<double>::mul_double(double x)
{
  store *= x;
}

template<>
void number<double>::add(const number<double> &other)
{
  store += other.store;
}

template<>
void number<double>::add(const number_a &other)
{
  add(*(const number<double> *)&other);
}

template<>
void number<double>::sub(const number<double> &other)
{
  store -= other.store;
}

template<>
void number<double>::sub(const number_a &other)
{
  sub(*(const number<double> *)&other);
}

template<>
void number<double>::rsub(const number<double> &other)
{
  store = other.store-store;
}

template<>
void number<double>::rsub(const number_a &other)
{
  rsub(*(const number<double> *)&other);
}

template<>
void number<double>::mul(const number<double> &other)
{
  store *= other.store;
}

template<>
void number<double>::mul(const number_a &other)
{
  mul(*(const number<double> *)&other);
}

template<>
void number<double>::sqr()
{
  store *= store;
}

template<>
double number<double>::radixfloor(const number<double> &other) const
{
  int ilog1=std::ilogb(store);
  int ilog2=std::ilogb(other.store);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

template<>
double number<double>::radixfloor(const number_a &other) const
{
  return radixfloor(*(const number<double> *)&other);
}

template<>
void number<double>::recip()
{
  store = 1/store;
}

template<>
void number<double>::sqrt()
{
  store = std::sqrt(store);
}

template<>
bool number<double>::reduceAngle()
{
  if (store<-M_PI)
  {
    store+=2*M_PI;
    return true;
  }
  else if (store>=M_PI)
  {
    store-=2*M_PI;
    return true;
  }
  else
    return false;
}

template<>
void number<double>::add_pi(double x)
{
  store+=x*M_PI;
}

template<>
int number<double>::compare(const number<double> &other) const
{
  if (store == other.store)
    return 0;
  else if (store < other.store)
    return -1;
  else
    return +1;
}

template<>
int number<double>::compare(const number_a &other) const
{
  return compare(*(const number<double> *)&other);
}

template<>
bool number<double>::isequal(const number<double> &other) const
{
  return store==other.store;
}

template<>
bool number<double>::isequal(const number_a &other) const
{
  return isequal(*(const number<double> *)&other);
}

template<>
bool number<double>::is0() const
{
  return store==0;
}

template<>
bool number<double>::isle(const number<double> &other)
{
  return store <= other.store;
}

template<>
bool number<double>::isle(const number_a &other)
{
  return isle(*(const number<double> *)&other);
}

template<>
bool number<double>::isle0() const
{
  return store<=0;
}

template<>
bool number<double>::isl0() const
{
  return store<0;
}

template<>
bool number<double>::isl1() const
{
  return store<1;
}

template<>
void number<double>::min(const number<double> &other)
{
  if (store>other.store)
  {
    store=other.store;
  };
}

template<>
void number<double>::min(const number_a &other)
{
  min(*(const number<double> *)&other);
}

template<>
QString number<double>::toString() const
{
  if (store<-1e6 || store>1e6)
    return QString::number(store, 'g', 10);
  else
    return QString::number(store, 'f', 16);
}

template<>
int number<double>::toRound() const
{
  return qRound(store);
}

template<>
double number<double>::toDouble() const
{
  return store;
}















#if !NUMBER_DOUBLE_ONLY

template<>
double number<__float128>::eps2() const
{
  return 9.27302e-69; // 2^-(2*113) rounded up
}

template<>
double number<__float128>::eps234() const
{
  return 9.45e-52; // eps2^(3/4)
}

template<>
number<__float128>::number(NumberType ntype): store(0)
{
  working_assert(ntype==NumberType::typeFloat128);
}

template<>
number<__float128>::~number()
{
  store=12.345;
}

template<>
NumberType number<__float128>::ntype() const
{
  return NumberType::typeFloat128;
}

template<>
NumberType number<__float128>::raw_ntype() const
{
  return NumberType::typeFloat128;
}

template<>
void number<__float128>::readFrom(void *storage, int index)
{
  store=((__float128 *)storage)[index];
}

template<>
void number<__float128>::writeTo(void *storage, int index) const
{
  ((__float128 *)storage)[index]=store;
}

template<>
__float128 *number<__float128>::convert_block(NumberType old_type, const void *old_data, int count)
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

template<>
void number<__float128>::zero(double val)
{
  store=val;
}

template<>
void number<__float128>::assign(const number<__float128> &src)
{
  store=src.store;
}

template<>
void number<__float128>::assign(const number_a &src)
{
  working_assert(src.raw_ntype()==NumberType::typeFloat128);
  //const number<double> *src2=(const number<double> *)&src;
  //assign(*src2);
  assign(*(const number<__float128> *)&src);
}

template<> //template <>
void number<__float128>::assign_across(const number<__float128> *src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;
  store=src->store;
}

/*template<> template <>
void number<__float128>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeFloat128);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<__float128> *access=specific_cast<number<__float128> *, number_a *>(src->store);
  store=access->store;
}*/

template<>
void number<__float128>::assign_across(const number_a *src)
{
  //working_assert(src->ntype()==typeFloat128);
  //working_assert(src->raw_ntype()==typeFloat128);
  switch (src->raw_ntype())
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
  }
}

template<>
void number<__float128>::chs()
{
  store=-store;
}

template<>
void number<__float128>::lshift(int shoft)
{
  //could try ldexp only on upper half but...
  if (store != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)&store)[7]+=shoft;
    //TODO: should check for overflow/underflow but unlikely to happen
}

template<>
void number<__float128>::round()
{
  __float128 remainder=store;
  store=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    store+=some;
    remainder-=some;
  }
}

template<>
void number<__float128>::frac()
{
  if (store<0)
  {
    double some=ceil((double)store);
    while (some)
    {
      store -= some;
      some=ceil((double)store);
    }
  }
  else
  {
    double some=floor((double)store);
    while (some)
    {
      store -= some;
      some=floor((double)store);
    }
  }
}

template<>
void number<__float128>::mod1()
{
  double some=floor((double)store);
  while (some)
  {
    store -= some;
    some=floor((double)store);
  }
}

template<>
void number<__float128>::add_double(double x)
{
  store += x;
}

template<>
void number<__float128>::mul_double(double x)
{
  store *= x;
}

template<>
void number<__float128>::add(const number<__float128> &other)
{
  store += other.store;
}

template<>
void number<__float128>::add(const number_a &other)
{
  add(*(const number<__float128> *)&other);
}

template<>
void number<__float128>::sub(const number<__float128> &other)
{
  store -= other.store;
}

template<>
void number<__float128>::sub(const number_a &other)
{
  sub(*(const number<__float128> *)&other);
}

template<>
void number<__float128>::rsub(const number<__float128> &other)
{
  store = other.store-store;
}

template<>
void number<__float128>::rsub(const number_a &other)
{
  rsub(*(const number<__float128> *)&other);
}

template<>
void number<__float128>::mul(const number<__float128> &other)
{
  store *= other.store;
}

template<>
void number<__float128>::mul(const number_a &other)
{
  mul(*(const number<__float128> *)&other);
}

template<>
void number<__float128>::sqr()
{
  store *= store;
}

template<>
double number<__float128>::radixfloor(const number<__float128> &other) const
{
  int ilog1=std::ilogb((double)store);
  int ilog2=std::ilogb((double)other.store);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

template<>
double number<__float128>::radixfloor(const number_a &other) const
{
  return radixfloor(*(const number<__float128> *)&other);
}

template<>
void number<__float128>::recip()
{
  store = __float128(1)/store;
}

template<>
void number<__float128>::sqrt()
{
  double x=1/std::sqrt((double)store);
  __float128 ax=x*(double)store;
  double tmp=store - ax*ax;
  double tmp2=tmp*x/2;
#if 0
  *store.asf128 = ax + tmp2;//double(tmp)*x/2;
  //error >30 eps because float128 is 112 bits but 2*double is 104 bits (256eps)
#else
  //one more iteration:
  ax+=tmp2;
  tmp=store - ax*ax;
  tmp2=tmp*x/2;
  store = ax + tmp2;
#endif
  //error analysis:
  //x=1/(sqrt + e)
  //ax=a/(sqrt + e)+f
  //result=a/(sqrt + e)+f + (a - (a/(sqrt + e)+f)*(a/(sqrt + e)+f) + g)*x/2
  //result=a/(sqrt + e)+f + (a - (a*a/(sqrt + e)/(sqrt + e)+2*f*a/(sqrt + e)+f*f) + g)*x/2
  //E=sqrt/(sqrt+e), a=sqrt*sqrt
  //result=sqrt*E+f + (a - a*E*E-2*f*sqrt*E-f*f + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + (a*(1-E*E)-2*f*sqrt + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + sqrt*(1-E*E)*E/2-f*E + g/sqrt*E/2
  //1-E*E(1-E)*(1+E)
  //result ~ sqrt*E+f + sqrt*(1+E)*(1-E)*E/2-f*E + g/sqrt*E/2
  //E=1+h
  //result ~ sqrt*(1+h)+f + sqrt*(2+h)*(-h)*(1+h)/2-f*(1+h) + g/sqrt*(1+h)/2
  //result ~ sqrt*(1+h) + sqrt*(-h)*(1+h)+sqrt*(-h)*(1+h)*h/2 + g/sqrt*(1+h)/2 -f*h
  //result ~ sqrt + g/sqrt/2 + second order terms
}

template<>
bool number<__float128>::reduceAngle()
{
  static constexpr __float128 M_PIq=3.141592653589793238462643383279502884197Q;
  if (store<-M_PIq)
  {
    store+=2*M_PIq;
    return true;
  }
  else if (store>=M_PIq)
  {
    store-=2*M_PIq;
    return true;
  }
  else
    return false;
}

template<>
void number<__float128>::add_pi(double x)
{
  static constexpr __float128 M_PIq=3.141592653589793238462643383279502884197Q;
  store+=x*M_PIq;
}

template<>
int number<__float128>::compare(const number<__float128> &other) const
{
  if (store == other.store)
    return 0;
  else if (store < other.store)
    return -1;
  else
    return +1;
}

template<>
int number<__float128>::compare(const number_a &other) const
{
  return compare(*(const number<__float128> *)&other);
}

template<>
bool number<__float128>::isequal(const number<__float128> &other) const
{
  return store==other.store;
}

template<>
bool number<__float128>::isequal(const number_a &other) const
{
  return isequal(*(const number<__float128> *)&other);
}

template<>
bool number<__float128>::is0() const
{
  return store==0;
}

template<>
bool number<__float128>::isle(const number<__float128> &other)
{
  return store <= other.store;
}

template<>
bool number<__float128>::isle(const number_a &other)
{
  return isle(*(const number<__float128> *)&other);
}

template<>
bool number<__float128>::isle0() const
{
  return store<=0;
}

template<>
bool number<__float128>::isl0() const
{
  return store<0;
}

template<>
bool number<__float128>::isl1() const
{
  return store<1;
}

template<>
void number<__float128>::min(const number<__float128> &other)
{
  if (store>other.store)
  {
    store=other.store;
  };
}

template<>
void number<__float128>::min(const number_a &other)
{
  min(*(const number<__float128> *)&other);
}

template<>
QString number<__float128>::toString() const
{
  return QString::number(store, 'f', 16);
}

template<>
int number<__float128>::toRound() const
{
  return qRound((double)store);
}

template<>
double number<__float128>::toDouble() const
{
  return (double)store;
}






template<>
double number<dd_real>::eps2() const
{
  return 6.07717e-64; // 2^-(2*(53+52)) rounded up
}

template<>
double number<dd_real>::eps234() const
{
  return 3.87e-48; // veps2^(3/4)
}

template<>
number<dd_real>::number(NumberType ntype): store()
{
  working_assert(ntype==NumberType::typeDDouble || ntype==NumberType::typeQDouble);
}

template<>
number<dd_real>::~number()
{
  store.hi=12.345;
  store.lo_=1.234e-30;
}

template<>
NumberType number<dd_real>::ntype() const
{
  return NumberType::typeDDouble;
}

template<>
NumberType number<dd_real>::raw_ntype() const
{
  return NumberType::typeDDouble;
}

template<>
void number<dd_real>::readFrom(void *storage, int index)
{
  store=((dd_real *)storage)[index];
}

template<>
void number<dd_real>::writeTo(void *storage, int index) const
{
  ((dd_real *)storage)[index]=store;
}

template<>
dd_real *number<dd_real>::convert_block(NumberType old_type, const void *old_data, int count)
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

template<>
void number<dd_real>::zero(double val)
{
  store.zero(val);
}

template<>
void number<dd_real>::assign(const number<dd_real> &src)
{
  store=src.store;
}

template<>
void number<dd_real>::assign(const number_a &src)
{
  working_assert(src.raw_ntype()==NumberType::typeDDouble);
  //const number<double> *src2=(const number<double> *)&src;
  //assign(*src2);
  assign(*(const number<dd_real> *)&src);
}

template<> //template <>
void number<dd_real>::assign_across(const number<dd_real> *src)
{
  store=src->store;
}

/*template<> template <>
void number<dd_real>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeDDouble);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<dd_real> *access=specific_cast<number<dd_real> *, number_a *>(src->store);
  store=access->store;
}*/

template<>
void number<dd_real>::assign_across(const number_a *src)
{
  //working_assert(src->ntype()==typeDDouble);
  //working_assert(src->raw_ntype()==typeDDouble);
  //store=specific_cast<const number<dd_real> *, const number_a *>(src)->store;
  switch (src->raw_ntype())
  {
    case NumberType::typeEmpty: assign_across(((const number<number_a *> *)src)->store); break;
    case NumberType::typeDouble: store.hi=((const number<double> *)src)->store; store.lo_=0; break;
    case NumberType::typeFloat128:
    {
      const __float128 &access=((const number<__float128> *)src)->store;
      store.hi=access;
      store.lo_=access-store.hi;
    } break;
    case NumberType::typeDDouble:
    {
      const dd_real &access=((const number<dd_real> *)src)->store;
      store=access;
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real &access=((const number<dq_real> *)src)->store;
      store=access;
    } break;
    case NumberType::typeReal642:
    {
      const real642 &access=((const number<real642> *)src)->store;
      store.hi=access.val1;
      store.lo_=0;
    } break;
  }
}

template<>
void number<dd_real>::chs()
{
  store.chs();
}

template<>
void number<dd_real>::lshift(int shoft)
{
  store.lshift(shoft);
}

template<>
void number<dd_real>::round()
{
  store.round();
}

template<>
void number<dd_real>::frac()
{
  store.frac();
}

template<>
void number<dd_real>::mod1()
{
  store.mod1();
}

template<>
void number<dd_real>::add_double(double x)
{
  store.add_double(x);
}

template<>
void number<dd_real>::mul_double(double x)
{
  store.mul_double(x);
}

template<>
void number<dd_real>::add(const number<dd_real> &other)
{
  store.add(other.store.hi, other.store.lo_);
}

template<>
void number<dd_real>::add(const number_a &other)
{
  add(*(const number<dd_real> *)&other);
}

template<>
void number<dd_real>::sub(const number<dd_real> &other)
{
  store.add(-other.store.hi, -other.store.lo_);
}

template<>
void number<dd_real>::sub(const number_a &other)
{
  sub(*(const number<dd_real> *)&other);
}

template<>
void number<dd_real>::rsub(const number<dd_real> &other)
{
  store.chs();
  store.add(other.store.hi, other.store.lo_);
}

template<>
void number<dd_real>::rsub(const number_a &other)
{
  rsub(*(const number<dd_real> *)&other);
}

template<>
void number<dd_real>::mul(const number<dd_real> &other)
{
  store.mul(other.store.hi, other.store.lo_);
}

template<>
void number<dd_real>::mul(const number_a &other)
{
  mul(*(const number<dd_real> *)&other);
}

template<>
void number<dd_real>::sqr()
{
  store.sqr();
}

template<>
double number<dd_real>::radixfloor(const number<dd_real> &other) const
{
  double rf1=store.radixfloor();
  double rf2=other.store.radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

template<>
double number<dd_real>::radixfloor(const number_a &other) const
{
  return radixfloor(*(const number<dd_real> *)&other);
}

template<>
void number<dd_real>::recip()
{
  store.recip();
}

template<>
void number<dd_real>::sqrt()
{
  store.sqrt();
}

template<>
bool number<dd_real>::reduceAngle()
{
  static const dd_real M_PIdd(3.141, 0.0005926);
  if (store.compare(-M_PIdd.hi, -M_PIdd.lo_)<0)
  {
    store.add(2*M_PIdd.hi, 2*M_PIdd.lo_);
    return true;
  }
  else if (store.isle(&M_PIdd))
  {
    store.add(-2*M_PIdd.hi, -2*M_PIdd.lo_);
    return true;
  }
  else
    return false;
}

template<>
void number<dd_real>::add_pi(double x)
{
  dd_real tmppi(3.141, 0.0005926);
  dbgPoint();
  tmppi.hi=3.141592653589793238462643383279502884197Q;
  tmppi.lo_=3.141592653589793238462643383279502884197Q-tmppi.hi;
  tmppi.mul_double(x);
  store.add(tmppi.hi, tmppi.lo_);
}

template<>
int number<dd_real>::compare(const number<dd_real> &other) const
{
  return store.compare(&other.store);
}

template<>
int number<dd_real>::compare(const number_a &other) const
{
  return compare(*(const number<dd_real> *)&other);
}

template<>
bool number<dd_real>::isequal(const number<dd_real> &other) const
{
  return store.isequal(&other.store);
}

template<>
bool number<dd_real>::isequal(const number_a &other) const
{
  return isequal(*(const number<dd_real> *)&other);
}

template<>
bool number<dd_real>::is0() const
{
  return store.is0();
}

template<>
bool number<dd_real>::isle(const number<dd_real> &other)
{
  return store.isle(&other.store);
}

template<>
bool number<dd_real>::isle(const number_a &other)
{
  return isle(*(const number<dd_real> *)&other);
}

template<>
bool number<dd_real>::isle0() const
{
  return store.isle0();
}

template<>
bool number<dd_real>::isl0() const
{
  return store.isl0();
}

template<>
bool number<dd_real>::isl1() const
{
  return store.isl1();
}

template<>
void number<dd_real>::min(const number<dd_real> &other)
{
  if (!store.isle(&other.store))
  {
    store=other.store;
  };
}

template<>
void number<dd_real>::min(const number_a &other)
{
  min(*(const number<dd_real> *)&other);
}

template<>
QString number<dd_real>::toString() const
{
  return QString("dd(%1,%2)").arg(store.hi, 0, 'f', 16).arg(store.lo_, 0, 'g', 16);
}

template<>
int number<dd_real>::toRound() const
{
  return floor(store.hi+0.5)+floor(store.lo_+0.5);
}

template<>
double number<dd_real>::toDouble() const
{
  return store.hi;
}






template<>
double number<real642>::eps2() const
{
  return 1.232596e-32; // 2^-(2*53) rounded up
}

template<>
double number<real642>::eps234() const
{
  return 1.17e-24; // eps2^(3/4)
}

template<>
number<real642>::number(NumberType ntype): store()
{
  working_assert(ntype==NumberType::typeReal642);
}

template<>
number<real642>::~number()
{
  store.val1=12.345;
  store.val2=12.345;
}

template<>
NumberType number<real642>::ntype() const
{
  return NumberType::typeReal642;
}

template<>
NumberType number<real642>::raw_ntype() const
{
  return NumberType::typeReal642;
}

template<>
void number<real642>::readFrom(void *storage, int index)
{
  store=((real642 *)storage)[index];
}

template<>
void number<real642>::writeTo(void *storage, int index) const
{
  ((real642 *)storage)[index]=store;
}

template<>
real642 *number<real642>::convert_block(NumberType old_type, const void *old_data, int count)
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

template<>
void number<real642>::zero(double val)
{
  store.val1=val;
  store.val2=val;
}

template<>
void number<real642>::assign(const number<real642> &src)
{
  store=src.store;
}

template<>
void number<real642>::assign(const number_a &src)
{
  working_assert(src.raw_ntype()==NumberType::typeReal642);
  //const number<double> *src2=(const number<double> *)&src;
  //assign(*src2);
  assign(*(const number<real642> *)&src);
}

template<> //template <>
void number<real642>::assign_across(const number<real642> *src)
{
  //const number<double> *access=specific_cast<number<double> *, number_a *>(src->store);
  //store=access->store;
  store=src->store;
}

/*template<> template <>
void number<real642>::assign_across<number_a *>(const number<number_a *> *src)
{
  working_assert(src->ntype()==typeReal642);
  working_assert(src->raw_ntype()==typeEmpty);
  const number<real642> *access=specific_cast<number<real642> *, number_a *>(src->store);
  store=access->store;
}*/

template<>
void number<real642>::assign_across(const number_a *src)
{
  //working_assert(src->ntype()==typeReal642);
  //working_assert(src->raw_ntype()==typeReal642);
  //store=specific_cast<const number<real642> *, const number_a *>(src)->store;
  switch (src->raw_ntype())
  {
    case NumberType::typeEmpty: assign_across(((const number<number_a *> *)src)->store); break;
    case NumberType::typeDouble: store.val2=store.val1=((const number<double> *)src)->store; break;
    case NumberType::typeFloat128:
    {
      const __float128 &access=((const number<__float128> *)src)->store;
      store.val1=access;
      store.val2=store.val1;
    } break;
    case NumberType::typeDDouble:
    {
      const dd_real &access=((const number<dd_real> *)src)->store;
      store.val1=access.hi;
      store.val2=access.hi;
    } break;
    case NumberType::typeQDouble:
    {
      const dq_real &access=((const number<dq_real> *)src)->store;
      store.val1=access.hi;
      store.val2=access.hi;
    } break;
    case NumberType::typeReal642:
    {
      const real642 &access=((const number<real642> *)src)->store;
      store=access;
    } break;
  }
}

template<>
void number<real642>::chs()
{
  store.val1=-store.val1;
  store.val2=-store.val2;
}

template<>
void number<real642>::lshift(int shoft)
{
  store.val1=ldexp(store.val1, shoft);
  if (store.val2 != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)&store.val2)[7]+=shoft;
    //TODO: should check for overflow/underflow but unlikely to happen
}

template<>
void number<real642>::round()
{
  store.val1=std::round(store.val1);
  __float128 remainder=store.val2;
  store.val2=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    store.val2+=some;
    remainder-=some;
  }
}

template<>
void number<real642>::frac()
{
  if (store.val1<0)
    store.val1 -= ceil(store.val1);
  else
    store.val1 -= floor(store.val1);
  if (store.val2<0)
  {
    double some=ceil((double)store.val2);
    while (some)
    {
      store.val2 -= some;
      some=ceil((double)store.val2);
    }
  }
  else
  {
    double some=floor((double)store.val2);
    while (some)
    {
      store.val2 -= some;
      some=floor((double)store.val2);
    }
  }
}

template<>
void number<real642>::mod1()
{
  store.val1 -= floor(store.val1);
  double some=floor((double)store.val2);
  while (some)
  {
    store.val2 -= some;
    some=floor((double)store.val2);
  }
}

template<>
void number<real642>::add_double(double x)
{
  store.val1 += x;
  store.val2 += x;
}

template<>
void number<real642>::mul_double(double x)
{
  store.val1 *= x;
  store.val2 *= x;
}

template<>
void number<real642>::add(const number<real642> &other)
{
  store.val1 += other.store.val1;
  store.val2 += other.store.val2;
}

template<>
void number<real642>::add(const number_a &other)
{
  add(*(const number<real642> *)&other);
}

template<>
void number<real642>::sub(const number<real642> &other)
{
  store.val1 -= other.store.val1;
  store.val2 -= other.store.val2;
}

template<>
void number<real642>::sub(const number_a &other)
{
  sub(*(const number<real642> *)&other);
}

template<>
void number<real642>::rsub(const number<real642> &other)
{
  store.val1 = other.store.val1-store.val1;
  store.val2 = other.store.val2-store.val2;
}

template<>
void number<real642>::rsub(const number_a &other)
{
  rsub(*(const number<real642> *)&other);
}

template<>
void number<real642>::mul(const number<real642> &other)
{
  store.val1 *= other.store.val1;
  store.val2 *= other.store.val2;
}

template<>
void number<real642>::mul(const number_a &other)
{
  mul(*(const number<real642> *)&other);
}

template<>
void number<real642>::sqr()
{
  store.val1 *= store.val1;
  store.val2 *= store.val2;
}

template<>
double number<real642>::radixfloor(const number<real642> &other) const
{
  int ilog1=std::ilogb(store.val1);
  int ilog2=std::ilogb(other.store.val1);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

template<>
double number<real642>::radixfloor(const number_a &other) const
{
  return radixfloor(*(const number<real642> *)&other);
}

template<>
void number<real642>::recip()
{
  store.val1 = 1.0/store.val1;
  store.val2 = __float128(1)/store.val2;
}

template<>
void number<real642>::sqrt()
{
  store.val1 = std::sqrt(store.val1);
  double x=1/std::sqrt((double)store.val2);
  __float128 ax=x*(double)store.val2;
  ax += double(store.val2 - ax*ax)*x/2;
  store.val2 = ax + double(store.val2 - ax*ax)*x/2;
}

template<>
bool number<real642>::reduceAngle()
{
  if (store.val1<-M_PI)
  {
    store.val1+=2*M_PI;
    store.val2+=2.0Q*3.141592653589793238462643383279502884197Q;
    return true;
  }
  else if (store.val1>=M_PI)
  {
    store.val1-=2*M_PI;
    store.val2-=2.0Q*3.141592653589793238462643383279502884197Q;
    return true;
  }
  else
    return false;
}

template<>
void number<real642>::add_pi(double x)
{
  store.val1+=x*M_PI;
  store.val2+=x*3.141592653589793238462643383279502884197Q;
}

template<>
int number<real642>::compare(const number<real642> &other) const
{
  if (store.val1 == other.store.val1)
    return 0;
  else if (store.val1 < other.store.val1)
    return -1;
  else
    return +1;
}

template<>
int number<real642>::compare(const number_a &other) const
{
  return compare(*(const number<real642> *)&other);
}

template<>
bool number<real642>::isequal(const number<real642> &other) const
{
  return store.val1==other.store.val1;
}

template<>
bool number<real642>::isequal(const number_a &other) const
{
  return isequal(*(const number<real642> *)&other);
}

template<>
bool number<real642>::is0() const
{
  return store.val1==0;
}

template<>
bool number<real642>::isle(const number<real642> &other)
{
  return store.val1 <= other.store.val1;
}

template<>
bool number<real642>::isle(const number_a &other)
{
  return isle(*(const number<real642> *)&other);
}

template<>
bool number<real642>::isle0() const
{
  return store.val1<=0;
}

template<>
bool number<real642>::isl0() const
{
  return store.val1<0;
}

template<>
bool number<real642>::isl1() const
{
  return store.val1<1;
}

template<>
void number<real642>::min(const number<real642> &other)
{
  if (store.val1>other.store.val1)
  {
    store=other.store;
  };
}

template<>
void number<real642>::min(const number_a &other)
{
  min(*(const number<real642> *)&other);
}

template<>
QString number<real642>::toString() const
{
  return QString::number(store.val1, 'f', 16);
}

template<>
int number<real642>::toRound() const
{
  return qRound(store.val1);
}

template<>
double number<real642>::toDouble() const
{
  return store.val1;
}

#endif // !NUMBER_DOUBLE_ONLY

template<>
double number<number_a *>::eps2() const
{
  return store->eps2();
}

template<>
double number<number_a *>::eps234() const
{
  return store->eps234();
}

template<>
number<number_a *>::number(NumberType ntype)
{
  switch (ntype)
  {
    //case NumberType::typeEmpty: dbgPoint(); store=new number<double>(); break;
    case NumberType::typeDouble: store=new number<double>(ntype); break;
#if !NUMBER_DOUBLE_ONLY
    case NumberType::typeFloat128: store=new number<__float128>(); break;
    case NumberType::typeDDouble: store=new number<dd_real>(); break;
    case NumberType::typeQDouble: store=new number<dq_real>(); break;
    case NumberType::typeReal642: store=new number<real642>(); break;
#endif
    default:
      dbgPoint(); store=new number<double>(); break;
  }
}

template<>
void number<number_a *>::constructLateBecauseQtIsAwesome(NumberType ntype)
{
  working_assert(store==nullptr);
  switch (ntype)
  {
    //case NumberType::typeEmpty: dbgPoint(); store=new number<double>(); break;
    case NumberType::typeDouble: store=new number<double>(ntype); break;
#if !NUMBER_DOUBLE_ONLY
    case NumberType::typeFloat128: store=new number<__float128>(); break;
    case NumberType::typeDDouble: store=new number<dd_real>(); break;
    case NumberType::typeQDouble: store=new number<dq_real>(); break;
    case NumberType::typeReal642: store=new number<real642>(); break;
#endif
    default:
      dbgPoint(); store=new number<double>(); break;
  }
}

template<>
number<number_a *>::~number()
{
  nop();
  delete store;
  store=nullptr;
}

template<>
void number<number_a *>::readFrom(void *storage, int index)
{
  store->readFrom(storage, index);
}

template<>
void number<number_a *>::writeTo(void *storage, int index) const
{
  store->writeTo(storage, index);
}

template<>
void number<number_a *>::zero(double val)
{
  store->zero(val);
}

template<>
void number<number_a *>::assign(const number<number_a *> &src)
{
  //working_assert(src.raw_ntype()==typeEmpty);
  working_assert(store);
  working_assert(src.store);
  working_assert(store->raw_ntype()==src.store->raw_ntype());
  store->assign(*src.store);
}

template<>
void number<number_a *>::assign(const number_a &src)
{
  (void)src;
  dbgPoint();
}

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

template<> //template<>
void number<number_a *>::assign_across(const number<number_a *> *src)
{
  working_assert(src->raw_ntype()==typeEmpty);
  //working_assert(store->ntype()==src->store->ntype());
  //working_assert(store->raw_ntype()==src->store->raw_ntype());
  //TODO: convert if store->ntype()!=src->store->ntype()
  store->assign_across(src->store);
}

template<>
void number<number_a *>::assign_across(const number_a *src)
{
  //working_assert(store->ntype()==src->store->ntype());
  //store->assign(*src->store);
  (void)src;
  dbgPoint();
}

template<>
void number<number_a *>::chs()
{
  store->chs();
}

template<>
void number<number_a *>::lshift(int shoft)
{
  store->lshift(shoft);
}

template<>
void number<number_a *>::round()
{
  store->round();
}

template<>
void number<number_a *>::frac()
{
  store->frac();
}

template<>
void number<number_a *>::mod1()
{
  store->mod1();
}

template<>
void number<number_a *>::add_double(double x)
{
  store->add_double(x);
}

template<>
void number<number_a *>::mul_double(double x)
{
  store->mul_double(x);
}

template<>
void number<number_a *>::add(const number<number_a *> &other)
{
  store->add(*other.store);
}

template<>
void number<number_a *>::add(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}

template<>
void number<number_a *>::sub(const number<number_a *> &other)
{
  store->sub(*other.store);
}

template<>
void number<number_a *>::sub(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}

template<>
void number<number_a *>::rsub(const number<number_a *> &other)
{
  store->rsub(*other.store);
}

template<>
void number<number_a *>::rsub(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}

template<>
void number<number_a *>::mul(const number<number_a *> &other)
{
  store->mul(*other.store);
}

template<>
void number<number_a *>::mul(const number_a &other)
{
  (void)other;
  dbgPoint();//store->add(*other.store);
}

template<>
void number<number_a *>::sqr()
{
  store->sqr();
}

template<>
double number<number_a *>::radixfloor(const number<number_a *> &store2) const
{
  return store->radixfloor(*store2.store);
}

template<>
double number<number_a *>::radixfloor(const number_a &store2) const
{
  (void)store2;
  dbgPoint();
  return 1;
}

template<>
void number<number_a *>::recip()
{
  store->recip();
}

template<>
void number<number_a *>::sqrt()
{
  store->sqrt();
}

template<>
bool number<number_a *>::reduceAngle()
{
  return store->reduceAngle();
}

template<>
void number<number_a *>::add_pi(double x)
{
  store->add_pi(x);
}

template<>
int number<number_a *>::compare(const number<number_a *> &other) const
{
  return store->compare(*other.store);
}

template<>
int number<number_a *>::compare(const number_a &other) const
{
  (void)other;
  dbgPoint();
  return 0;
}

template<>
bool number<number_a *>::isequal(const number<number_a *> &other) const
{
  return store->isequal(*other.store);
}

template<>
bool number<number_a *>::isequal(const number_a &other) const
{
  (void)other;
  dbgPoint();
  return false;
}

template<>
bool number<number_a *>::is0() const
{
  return store->is0();
}

template<>
bool number<number_a *>::isle(const number<number_a *> &other)
{
  return store->isle(*other.store);
}

template<>
bool number<number_a *>::isle(const number_a &other)
{
  (void)other;
  dbgPoint();
  return true;
}

template<>
bool number<number_a *>::isle0() const
{
  return store->isle0();
}

template<>
bool number<number_a *>::isl0() const
{
  return store->isl0();
}

template<>
bool number<number_a *>::isl1() const
{
  return store->isl1();
}

template<>
void number<number_a *>::min(const number<number_a *> &other)
{
  store->min(*other.store);
}

template<>
void number<number_a *>::min(const number_a &other)
{
  (void)other;
  dbgPoint();
}

template<>
QString number<number_a *>::toString() const
{
  return store->toString();
}

template<>
int number<number_a *>::toRound() const
{
  return store->toRound();
}

template<>
double number<number_a *>::toDouble() const
{
  return store->toDouble();
}




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
void complex<BASE>::zero(double r, double i)
{
  re.zero(r);
  im.zero(i);
}

template<typename BASE>
void complex<BASE>::assign(const complex *other)
{
  re.assign(other->re);
  im.assign(other->im);
}

template<typename BASE> template<typename OTHER_BASE>
void complex<BASE>::assign_across(const complex<OTHER_BASE> *src)
{
  re.assign_across(&src->re);
  im.assign_across(&src->im);
}

template<typename BASE>
void complex<BASE>::lshift(int shoft)
{
  re.lshift(shoft);
  im.lshift(shoft);
}

template<typename BASE>
void complex<BASE>::mul_double(double m)
{
  re.mul_double(m);
  im.mul_double(m);
}

template<typename BASE>
void complex<BASE>::add(const complex *const other)
{
  re.add(other->re);
  im.add(other->im);
}

template<typename BASE>
void complex<BASE>::chs()
{
  re.chs();
  im.chs();
}

template<typename BASE>
void complex<BASE>::sub(const complex *other)
{
  re.sub(other->re);
  im.sub(other->im);
}

template<typename BASE>
void complex<BASE>::rsub(const complex *other)
{
  re.rsub(other->re);
  im.rsub(other->im);
}

template<typename BASE>
void complex<BASE>::sqr(Scratchpad *tmp)
{
#if 0
  //r:=r*r-i*i
  //i=2*r*i
  worker->assign(tmp1, im);
  worker->sqr(tmp1);
  worker->mul(im, re);
  worker->lshift(im, 1);
  worker->sqr(re);
  worker->sub(re, tmp1);
#else
  //r:=(r+i)*(r-i)
  //i=2*r*i
  //should be more precise (+ fewer mul but sqr is cheaper than mul so not really)
  tmp->tmp1.assign(re);
  tmp->tmp1.sub(im);    //r-i
  tmp->tmp2.assign(im);
  im.mul(re);
  im.lshift(1); //2*r*i
  re.add(tmp->tmp2);    //r+i
  re.mul(tmp->tmp1);
#endif
}

template<typename BASE>
void complex<BASE>::mul(const number<BASE> &other)
{
  re.mul(other);
  im.mul(other);
}

template<typename BASE>
void complex<BASE>::mul(const complex *other, Scratchpad *tmp)
{
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  re.mul(other->im);
  im.mul(other->re);
  tmp->tmp1.mul(other->re); //r1*r2
  tmp->tmp2.mul(other->im); //i1*i2
  tmp->tmp1.sub(tmp->tmp2); //r1*r2-i1*i2
  im.add(re); //i1*r2+r1*i2
  re.assign(tmp->tmp1);
}

template<typename BASE>
void complex<BASE>::recip(Scratchpad *tmp)
{
  getMag_tmp(tmp);
  recip_prepared(tmp);
}

template<typename BASE>
void complex<BASE>::recip_prepared(Scratchpad *tmp)
{ // 1/(re+i*im) = (re-i*im)/((re+i*im)*(re-i*im)) = (re-i*im)/(re*re+im*im)
  tmp->tmp1.recip();
  im.chs();
  re.mul(tmp->tmp1);
  im.mul(tmp->tmp1);
}

template<typename BASE>
void complex<BASE>::sqrt(Scratchpad *tmp)
{
  /*tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  tmp->tmp1.sqr();
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2); //re*re+im*im*/
  getMag_tmp(tmp);
  tmp->tmp1.sqrt();
  if (!re.isle0())
  {
    tmp->tmp1.add(re);
    tmp->tmp1.lshift(-1);
    tmp->tmp1.sqrt(); //t1=sqrt((sqrt(re*re+im*im)+re)/2);
    re.assign(tmp->tmp1); //re=t1
    /*if (t1==0)
      *res_im=0;
    else*/
    tmp->tmp1.lshift(1);
    tmp->tmp1.recip();
    im.mul(tmp->tmp1); //im=im/(2*t1);
  }
  else
  {
    tmp->tmp1.sub(re);
    tmp->tmp1.lshift(-1);
    tmp->tmp1.sqrt(); //t1=sqrt((sqrt(re*re+im*im)-re)/2);
    if (tmp->tmp1.isle0()) //t1==0
    {
      re.zero();
      im.zero();
    }
    else
    {
      re.assign(im);
      im.assign(tmp->tmp1); //new im=t1
      tmp->tmp1.lshift(1);
      tmp->tmp1.recip();
      re.mul(tmp->tmp1); //re=old im/(2*t1);
      if (re.isl0())
      {
        re.chs();
        im.chs();
      };
    }
  };
}

template<typename BASE>
void complex<BASE>::root_approx(int n)
{
  if (is0())
    return;
  double this_re=re.toDouble();
  double this_im=im.toDouble();
  double r_abs=std::exp(std::log(this_re*this_re+this_im*this_im)/(2*n));
  double r_phi=std::atan2(this_im, this_re)/n;
  double r_re=r_abs*cos(r_phi);
  double r_im=r_abs*sin(r_phi);
  zero(r_re, r_im);
}

template<typename BASE>
void complex<BASE>::ln_approx()
{
  double r=re.toDouble();
  double i=im.toDouble();
  zero(std::log(r*r+i*i)/2, std::atan2(i, r));
}

template<typename BASE>
void complex<BASE>::exp_approx()
{
  double d=exp(re.toDouble());
  double p=im.toDouble();
  zero(d*std::cos(p), d*std::sin(p));
}

template<typename BASE>
void complex<BASE>::sign(Scratchpad *tmp)
{
  getMag_tmp(tmp);
  tmp->tmp1.sqrt();
  tmp->tmp1.recip();
  re.mul(tmp->tmp1);
  im.mul(tmp->tmp1);
}

template<typename BASE>
void complex<BASE>::cossin(number<BASE> &angle, Scratchpad *tmp)
{
  zero(cos(angle.toDouble()), sin(angle.toDouble()));
  (void)tmp;
}

template<typename BASE>
void complex<BASE>::arctan2(number<BASE> *result, Scratchpad *tmp) const
{
  if (im.is0() && re.isl0())
    result->zero(std::atan2(-0, 1));
  else
    result->zero(std::atan2(im.toDouble(), re.toDouble()));
  (void)tmp;
}

template<typename BASE>
double complex<BASE>::getMag_double() const
{/*
  re.worker->assign(tmp1, re.getptr());
  re.worker->sqr(tmp1);
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return re.worker->toDouble(tmp1);*/
  return sqr_double(re.toDouble())+sqr_double(im.toDouble());
}

template<typename BASE>
const number<BASE> *complex<BASE>::getMag_tmp(Scratchpad *tmp) const
{
  tmp->tmp1.assign(re);
  tmp->tmp1.sqr();
  tmp->tmp2.assign(im);
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2);
  return &tmp->tmp1;
}

template<typename BASE>
const number<BASE> *complex<BASE>::getMag1_tmp(Scratchpad *tmp) const
{
  //TODO: merge with mag_cmp_1
  tmp->tmp1.assign(re);
  tmp->tmp1.sqr();
  tmp->tmp2.assign(im);
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2);
  tmp->tmp1.add_double(-1);
  return &tmp->tmp1;
}

template<typename BASE>
const number<BASE> *complex<BASE>::getDist1_tmp(Scratchpad *tmp) const
{
  tmp->tmp1.assign(re);
  tmp->tmp1.add_double(-1);
  tmp->tmp1.sqr();
  tmp->tmp2.assign(im);
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2);
  return &tmp->tmp1;
}

template<typename BASE>
const number<BASE> *complex<BASE>::mulreT_tmp(const complex *other, Scratchpad *tmp) const
{
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  tmp->tmp1.mul(other->re);
  tmp->tmp2.mul(other->im);
  tmp->tmp1.add(tmp->tmp2);
  return &tmp->tmp1;
}

template<typename BASE>
const number<BASE> *complex<BASE>::ccw_tmp(const complex *other, Scratchpad *tmp) const
{
  //im(other/this)=im((o.re+i*o.im)/(t.re+i*t.im))=im((o.re+i*o.im)*(t.re-i*t.im)/((t.re+i*t.im)*(t.re-i*t.im)))=
  //im((o.re*t.re+o.im*t.im+i*(o.im*t.re-o.re*t.im))/(t.re*t.re+t.im*t.im))=
  //(t.re*o.im-t.im*o.re)/(t.re*t.re+t.im*t.im)
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  tmp->tmp1.mul(other->im);
  tmp->tmp2.mul(other->re);
  tmp->tmp1.sub(tmp->tmp2); //re*o.im-im*o.re
  return &tmp->tmp1;
}

template<typename BASE>
double complex<BASE>::dist2_double(const complex *other, Scratchpad *tmp) const
{
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  tmp->tmp1.sub(other->re);
  tmp->tmp2.sub(other->im);
  /*tmp->tmp1.sqr();
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2);
  return tmp->tmp1.toDouble();*/
  return sqr_double(tmp->tmp1.toDouble())+sqr_double(tmp->tmp2.toDouble());
}

template<typename BASE>
const number<BASE> *complex<BASE>::dist2_tmp(const complex *other, Scratchpad *tmp) const
{
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  tmp->tmp1.sub(other->re);
  tmp->tmp2.sub(other->im);
  tmp->tmp1.sqr();
  tmp->tmp2.sqr();
  tmp->tmp1.add(tmp->tmp2);
  return &tmp->tmp1;
}

template<typename BASE>
bool complex<BASE>::isequal(const complex *other) const
{
  return re.isequal(other->re) &&
      im.isequal(other->im);
}

template<typename BASE>
bool complex<BASE>::is0() const
{
  return re.is0() && im.is0();
}

template<typename BASE>
bool complex<BASE>::isNegative() const
{
  //im<0 || im==0 && re<0
  return im.isl0() || (im.is0() && re.isl0());
}

template<typename BASE>
int complex<BASE>::mag_cmp_1(Scratchpad *tmp) const
{
  //first reduce fz_re, fz_im to 2nd octant:
  //  fz.re>=0, fz.im>=0, fz.im>=fz.re
  tmp->tmp1.assign(re);
  tmp->tmp2.assign(im);
  if (tmp->tmp1.isl0()) tmp->tmp1.chs();
  if (tmp->tmp2.isl0()) tmp->tmp2.chs();
  if (!tmp->tmp1.isle(tmp->tmp2))
  {
    tmp->tmp3.assign(tmp->tmp1);
    tmp->tmp1.assign(tmp->tmp2);
    tmp->tmp2.assign(tmp->tmp3);
  };
  //im>=re, now check re^2+im^2>1
  if ((tmp->tmp1.toDouble()>=0.71) || //im>=re>=0.71 -> mag>=1.0082
      (tmp->tmp2.toDouble()>=1.01))   //im>=1.01 -> mag>=1.01
  {
    return 1;
  }
  else if (tmp->tmp2.toDouble()<=0.70) //0.70>=im>=re -> mag<=0.98
    return -1;
  else
  { //   re*re+im*im>1
    //   gets inaccurate for re~0, im~1
    //   im*im>1-re*re ; im>sqrt(1-re*re) ; im-1>sqrt(1-re*re)-1
    //>> x=Sqrt(1-r*r)-1
    //   could try x:=((-r2/16-1/8)*r2-1/2)*r2; //-r^6/16-r^4/8-r^2/2
    //   but one cycle of Newton should work with any precision
    //   (x+1)*(x+1)=(1-r*r)
    //   x*x+2*x+r*r=0  f'=2*x+2
    //   x2=x-(x*x+2*x+r*r)/(2*x+2)
    //>> x2=(x*x-r*r)/(2*x+2)
    //   x2=(x-1)(r*r-x*x)/(1-x*x)/2
    tmp->tmp1.sqr(); //r^2
    tmp->tmp1.chs();
    tmp->tmp3.assign(tmp->tmp1);
    tmp->tmp3.add_double(1);
    tmp->tmp3.sqrt();              //sqrt(1-r*r) = x+1
    tmp->tmp4.assign(tmp->tmp3);
    tmp->tmp4.add_double(-1);        //x
    tmp->tmp4.sqr();
    tmp->tmp1.add(tmp->tmp4);    //x*x-r*r
    tmp->tmp3.recip(); //-1..-0.51 -> -1..-1.96
    tmp->tmp1.mul(tmp->tmp3);
    tmp->tmp1.lshift(-1);           //(x*x-r*r)/(x+1)/2 = better x
    tmp->tmp2.add_double(-1);
    return tmp->tmp2.compare(tmp->tmp1); //im-1 vs limit
  };
}

template<typename BASE>
QString complex<BASE>::toString() const
{
  return re.toString()+" +i* "+im.toString();
}
































#if 0

template<> worker_multi::Allocator<worker_multi>::Allocator(int capacity):
  lowest(0), allocUpTo(0), capacity(capacity), above(nullptr), worker(nullptr)
{
}

template<class WORKER_MULTI>
MandelMath::worker_multi::Allocator<WORKER_MULTI>::Allocator(MandelMath::worker_multi *worker, int capacity):
  lowest(0), allocUpTo(0), capacity(capacity), above(nullptr), worker((WORKER_MULTI *)worker)
{
}

template<> worker_multi::Allocator<worker_multi>::Allocator(Allocator *allo, int capacity):
  lowest(allo->allocUpTo), allocUpTo(allo->allocUpTo), capacity(capacity), above(allo), worker(allo->worker)
{
  allo->allocUpTo+=capacity;
  working_assert(allo->allocUpTo<=allo->lowest+allo->capacity);
}

template<> worker_multi::Allocator<worker_multi>::Allocator(Allocator *allo, int lowest, int count):
  lowest(allo->lowest+lowest), allocUpTo(allo->lowest+lowest), capacity(count), above(nullptr), worker(allo->worker)
{
  working_assert(lowest>=allo->lowest);
  working_assert(lowest+capacity<=allo->allocUpTo);
}

template<class WORKER_MULTI>
MandelMath::worker_multi::Allocator<WORKER_MULTI>::~Allocator()
{
   working_assert(allocUpTo==lowest);
   if (above!=nullptr)
   {
     working_assert(above->allocUpTo==lowest+capacity);
     working_assert(above->allocUpTo-above->lowest>=capacity);
     above->allocUpTo-=capacity;
   };
}

template<class WORKER_MULTI>
number_index worker_multi::Allocator<WORKER_MULTI>::alloc()
{
  if (worker==nullptr) //support for Qt that requires parameter-less constructors for ShareableViewInfo
  {
    dbgPoint();
    return number_index(-1, nullptr);
  };
  working_assert(allocUpTo<lowest+capacity);
  allocUpTo++;
  return worker->getNumber(allocUpTo-1);
}

template<class WORKER_MULTI>
void worker_multi::Allocator<WORKER_MULTI>::dealloc(number_index ptr)
{
  if (worker==nullptr) //support for Qt that requires parameter-less constructors for ShareableViewInfo
    return;
  working_assert(ptr.dbg_worker==worker);
  working_assert(allocUpTo>lowest);
  allocUpTo--;
  number_index expected=worker->getNumber(allocUpTo);
  working_assert((int)expected==(int)ptr);
}


worker_multi *worker_multi::create(Type ntype, int capacity)
{
  switch (ntype)
  {
    default://case Type::typeEmpty:
      dbgPoint();
      goto lolwut;
#if NUMBER_DOUBLE_EXISTS
    case Type::typeDouble: lolwut:
      return new worker_multi_double(capacity);
#endif
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
      return new worker_multi_float128(capacity);
    case Type::typeDDouble:
      return new worker_multi_ddouble(capacity);
    case Type::typeQDouble:
      return new worker_multi_qdouble(capacity);
    case Type::typeReal642:
      return new worker_multi_real642(capacity);
#endif
  }
}

worker_multi *worker_multi::create(Type ntype, Allocator<worker_multi> *allocator)
{
  switch (ntype)
  {
    default://case Type::typeEmpty:
      dbgPoint();
      goto lolwut;
#if NUMBER_DOUBLE_EXISTS
    case Type::typeDouble: lolwut:
      return new worker_multi_double(allocator);
#endif
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
      return new worker_multi_float128(allocator);
    case Type::typeDDouble:
      return new worker_multi_ddouble(allocator);
    case Type::typeQDouble:
      return new worker_multi_qdouble(allocator);
    case Type::typeReal642:
      return new worker_multi_real642(allocator);
#endif
  }
}


#if NUMBER_DOUBLE_EXISTS
//worker_multi_double::worker_multi_double(worker_multi *source):
worker_multi_double::worker_multi_double(Allocator<worker_multi> *source):
  worker_multi(Type::typeDouble, specific_cast<Allocator<worker_multi_double> *, Allocator<worker_multi> *>(source)->capacity)
{
  int oldfirst, capacity;
  source->_getFirstCapac(oldfirst, capacity);
  storage_=new double[iiw__COUNT+capacity];
  for (int i=0; i<iiw__COUNT; i++)
    storage_[i]=0;
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage_[iiw__COUNT+i]=0;
    } break;
    case Type::typeDouble:
    {
      double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->storage_+iiw__COUNT;
      src_storage+=oldfirst;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=*src_storage;
        src_storage++;
      }
    } break;
#if !ONLY_DOUBLE_WORKER
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=*src_storage;
        src_storage++;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->hi;
        src_storage++;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->hi;
        src_storage++;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->val1;
        src_storage++;
      }
    } break;
#endif
  }
}

worker_multi_double::~worker_multi_double()
{
  delete[] storage_;
  storage_=nullptr;
}

number_index worker_multi_double::getNumber(int index)
{
  working_assert(allocator.checkIndex(index));
  //return number_pointer(&storage_[COUNT_OF_TMPS+index]);
  return number_index(iiw__COUNT+index, this);
}

/*void worker_multi_double::getTmp12(number_index &t1, number_index &t2)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
}

void worker_multi_double::getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  working_assert(t3.dbg_worker==this);
  working_assert(t4.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
  t3.index_into_worker=2;
  t4.index_into_worker=3;
}*/

void worker_multi_double::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  working_assert(src_worker->ntype()==Type::typeDouble);
  worker_multi_double *access=(worker_multi_double *)src_worker;
  working_assert(allocator.checkIndex(dst));
  working_assert(allocator.checkIndex(dst+len-1));
  working_assert(access->allocator.checkIndex(src));
  working_assert(access->allocator.checkIndex(src+len-1));
  //working_assert(dst>=0 && dst+len<=capacity);
  //working_assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage_[iiw__COUNT+dst+i]=access->storage_[iiw__COUNT+src+i];
}

void worker_multi_double::zero(const number_index store, double val)
{
  working_assert(store.dbg_worker==this);
  storage_[store]=val;
}

void worker_multi_double::assign(const number_index store, const number_index src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src.dbg_worker==this);
  working_assert(allocator.checkIndex(store));
  working_assert(allocator.checkIndex(src));
  storage_[store]=storage_[src];
}

void worker_multi_double::assign_across(const number_index store, const number_place *src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src->ptr.dbg_worker->ntype()==Type::typeDouble);
  const worker_multi_double *access=specific_cast<const worker_multi_double *, const worker_multi *>(src->allocator->worker);
  working_assert((src->ptr.dbg_worker==access));
  working_assert(allocator.checkIndex(store));
  working_assert(src->allocator->checkIndex(src->ptr));
  storage_[store]=access->storage_[src->ptr];
}

void worker_multi_double::chs(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] = -storage_[store];
}

void worker_multi_double::lshift(const number_index store, int shoft)
{
  working_assert(store.dbg_worker==this);
  storage_[store]=ldexp(storage_[store], shoft);
}

void worker_multi_double::round(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store]=std::round(storage_[store]);
}

void worker_multi_double::frac(const number_index store)
{
  working_assert(store.dbg_worker==this);
  if (storage_[store]<0)
    storage_[store] -= ceil(storage_[store]);
  else
    storage_[store] -= floor(storage_[store]);
}

void worker_multi_double::mod1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] -= floor(storage_[store]);
}

void worker_multi_double::add_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store] += x;
}

void worker_multi_double::mul_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store] *= x;
}

void worker_multi_double::add(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] += storage_[other];
}

void worker_multi_double::sub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] -= storage_[other];
}

void worker_multi_double::rsub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] = storage_[other]-storage_[store];
}

void worker_multi_double::mul(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] *= storage_[other];
}

void worker_multi_double::sqr(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] *= storage_[store];
}

double worker_multi_double::radixfloor(const number_index store1, const number_index store2)
{
  working_assert(store1.dbg_worker==this);
  working_assert(store2.dbg_worker==this);
  int ilog1=std::ilogb(storage_[store1]);
  int ilog2=std::ilogb(storage_[store2]);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_double::recip(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] = 1/(storage_[store]);
}

void worker_multi_double::sqrt(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] = std::sqrt(storage_[store]);
}

int worker_multi_double::compare(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store] == storage_[other])
    return 0;
  else if (storage_[store] < storage_[other])
    return -1;
  else
    return +1;
}

bool worker_multi_double::isequal(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store] == storage_[other];
}

bool worker_multi_double::is0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] == 0;
}

bool worker_multi_double::isle(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store] <= storage_[other];
}

bool worker_multi_double::isle0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] <= 0;
}

bool worker_multi_double::isl0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] < 0;
}

bool worker_multi_double::isl1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] < 1;
}

void worker_multi_double::min(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store]>storage_[other])
  {
    storage_[store]=storage_[other];
  };
}

QString worker_multi_double::toString(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return QString::number(storage_[store], 'f', 16);
}

int worker_multi_double::toRound(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return qRound(storage_[store]);
}

double worker_multi_double::toDouble(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store];
}
#endif //NUMBER_DOUBLE_EXISTS


#if !ONLY_DOUBLE_WORKER

worker_multi_float128::worker_multi_float128(Allocator<worker_multi> *source):
  worker_multi(Type::typeFloat128, specific_cast<Allocator<worker_multi_float128> *, Allocator<worker_multi> *>(source)->capacity)
{
  int oldfirst, capacity;
  source->_getFirstCapac(oldfirst, capacity);
  storage_=new __float128[iiw__COUNT+capacity];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage_[iiw__COUNT+i]=0;
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=*src_storage;
        src_storage++;
      }
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=*src_storage;
        src_storage++;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i]=src_storage->val1;
        src_storage++;
      }
    } break;
  }
}

worker_multi_float128::~worker_multi_float128()
{
  delete[] storage_;
  storage_=nullptr;
}

number_index worker_multi_float128::getNumber(int index)
{
  working_assert(allocator.checkIndex(index));
  //working_assert(index>=0 && index<capacity);
  return number_index(iiw__COUNT+index, this);
}

/*void worker_multi_float128::getTmp12(number_index &t1, number_index &t2)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
}

void worker_multi_float128::getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  working_assert(t3.dbg_worker==this);
  working_assert(t4.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
  t3.index_into_worker=2;
  t4.index_into_worker=3;
}*/

void worker_multi_float128::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  working_assert(src_worker->ntype()==Type::typeFloat128);
  worker_multi_float128 *access=(worker_multi_float128 *)src_worker;
  working_assert(allocator.checkIndex(dst));
  working_assert(allocator.checkIndex(dst+len-1));
  working_assert(access->allocator.checkIndex(src));
  working_assert(access->allocator.checkIndex(src+len-1));
  //working_assert(dst>=0 && dst+len<=capacity);
  //working_assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage_[iiw__COUNT+dst+i]=access->storage_[iiw__COUNT+src+i];
}

void worker_multi_float128::zero(const number_index store, double val)
{
  working_assert(store.dbg_worker==this);
  storage_[store]=val;
}

void worker_multi_float128::assign(const number_index store, const number_index src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src.dbg_worker==this);
  working_assert(allocator.checkIndex(store));
  working_assert(allocator.checkIndex(src));
  storage_[store]=storage_[src];
}

void worker_multi_float128::assign_across(const number_index store, const number_place *src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src->ptr.dbg_worker->ntype()==Type::typeFloat128);
  const worker_multi_float128 *access=specific_cast<const worker_multi_float128 *, const worker_multi *>(src->allocator->worker);
  working_assert((src->ptr.dbg_worker==access));
  working_assert(allocator.checkIndex(store));
  working_assert(src->allocator->checkIndex(src->ptr));
  storage_[store]=access->storage_[src->ptr];
}

void worker_multi_float128::chs(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] = -storage_[store];
}

void worker_multi_float128::lshift(const number_index store, int shoft)
{
  working_assert(store.dbg_worker==this);
  //*store.asf128=ldexp(*store.asf128, shoft);
  //could try ldexp only on upper half but...
  if (storage_[store] != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)&storage_[store])[7]+=shoft;
}

void worker_multi_float128::round(const number_index store)
{
  working_assert(store.dbg_worker==this);
  __float128 remainder=storage_[store];
  storage_[store]=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    storage_[store]+=some;
    remainder-=some;
  }
}

void worker_multi_float128::frac(const number_index store)
{
  working_assert(store.dbg_worker==this);
  if (storage_[store]<0)
  {
    double some=ceil((double)storage_[store]);
    while (some)
    {
      storage_[store] -= some;
      some=ceil((double)storage_[store]);
    }
  }
  else
  {
    double some=floor((double)storage_[store]);
    while (some)
    {
      storage_[store] -= some;
      some=floor((double)storage_[store]);
    }
  }
}

void worker_multi_float128::mod1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  double some=floor((double)storage_[store]);
  while (some)
  {
    storage_[store] -= some;
    some=floor((double)storage_[store]);
  }
}

void worker_multi_float128::add_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store] += x;
}

void worker_multi_float128::mul_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store] *= x;
}

void worker_multi_float128::add(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] += storage_[other];
}

void worker_multi_float128::sub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] -= storage_[other];
}

void worker_multi_float128::rsub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] = storage_[other]-storage_[store];
}

void worker_multi_float128::mul(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store] *= storage_[other];
}

void worker_multi_float128::sqr(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] *= storage_[store];
}

double worker_multi_float128::radixfloor(const number_index store1, const number_index store2)
{
  working_assert(store1.dbg_worker==this);
  working_assert(store2.dbg_worker==this);
  int ilog1=std::ilogb((double)storage_[store1]);
  int ilog2=std::ilogb((double)storage_[store2]);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_float128::recip(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store] = __float128(1)/(storage_[store]);
}

void worker_multi_float128::sqrt(const number_index store)
{
  working_assert(store.dbg_worker==this);
  double x=1/std::sqrt((double)storage_[store]);
  __float128 ax=x*(double)storage_[store];
  double tmp=storage_[store] - ax*ax;
  double tmp2=tmp*x/2;
#if 0
  *store.asf128 = ax + tmp2;//double(tmp)*x/2;
  //error >30 eps because float128 is 112 bits but 2*double is 104 bits (256eps)
#else
  //one more iteration:
  ax+=tmp2;
  tmp=storage_[store] - ax*ax;
  tmp2=tmp*x/2;
  storage_[store] = ax + tmp2;
#endif
  //x=1/(sqrt + e)
  //ax=a/(sqrt + e)+f
  //result=a/(sqrt + e)+f + (a - (a/(sqrt + e)+f)*(a/(sqrt + e)+f) + g)*x/2
  //result=a/(sqrt + e)+f + (a - (a*a/(sqrt + e)/(sqrt + e)+2*f*a/(sqrt + e)+f*f) + g)*x/2
  //E=sqrt/(sqrt+e), a=sqrt*sqrt
  //result=sqrt*E+f + (a - a*E*E-2*f*sqrt*E-f*f + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + (a*(1-E*E)-2*f*sqrt + g)*1/(sqrt + e)/2
  //result ~ sqrt*E+f + sqrt*(1-E*E)*E/2-f*E + g/sqrt*E/2
  //1-E*E(1-E)*(1+E)
  //result ~ sqrt*E+f + sqrt*(1+E)*(1-E)*E/2-f*E + g/sqrt*E/2
  //E=1+h
  //result ~ sqrt*(1+h)+f + sqrt*(2+h)*(-h)*(1+h)/2-f*(1+h) + g/sqrt*(1+h)/2
  //result ~ sqrt*(1+h) + sqrt*(-h)*(1+h)+sqrt*(-h)*(1+h)*h/2 + g/sqrt*(1+h)/2 -f*h
  //result ~ sqrt + g/sqrt/2 + second order terms
}

int worker_multi_float128::compare(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store] == storage_[other])
    return 0;
  else if (storage_[store] < storage_[other])
    return -1;
  else
    return +1;
}

bool worker_multi_float128::isequal(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store] == storage_[other];
}

bool worker_multi_float128::is0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] == 0;
}

bool worker_multi_float128::isle(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store] <= storage_[other];
}

bool worker_multi_float128::isle0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] <= 0;
}

bool worker_multi_float128::isl0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] < 0;
}

bool worker_multi_float128::isl1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store] < 1; //enough to check top 16-17 bits?
}

void worker_multi_float128::min(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store]>storage_[other])
  {
    storage_[store]=storage_[other];
  };
}

QString worker_multi_float128::toString(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return QString::number(storage_[store], 'f', 16);
}

int worker_multi_float128::toRound(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return qRound((double)storage_[store]);
}

double worker_multi_float128::toDouble(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return (double)storage_[store];
}







number_index worker_multi_ddouble::getNumber(int index)
{
  working_assert(allocator.checkIndex(index));
  //working_assert(index>=0 && index<capacity);
  return number_index(iiw__COUNT+index, this);
}

/*void worker_multi_ddouble::getTmp12(number_index &t1, number_index &t2)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
}

void worker_multi_ddouble::getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  working_assert(t3.dbg_worker==this);
  working_assert(t4.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
  t3.index_into_worker=2;
  t4.index_into_worker=3;
}*/

worker_multi_ddouble::worker_multi_ddouble(Allocator<worker_multi> *source):
  worker_multi(Type::typeDDouble, specific_cast<Allocator<worker_multi_ddouble> *, Allocator<worker_multi> *>(source)->capacity)
{
  int oldfirst, capacity;
  source->_getFirstCapac(oldfirst, capacity);
  storage_=new dd_real[iiw__COUNT+capacity];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage_[iiw__COUNT+i].zero(0);
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].zero(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].hi=(double)*src_storage;
        storage_[iiw__COUNT+i].lo_=*src_storage-storage_[iiw__COUNT+i].hi;
        src_storage++;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].hi=src_storage->val1;
        storage_[iiw__COUNT+i].lo_=0;
        src_storage++;
      }
    } break;
  }
}

worker_multi_ddouble::~worker_multi_ddouble()
{
  delete[] storage_;
  storage_=nullptr;
}

void worker_multi_ddouble::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  working_assert(src_worker->ntype()==Type::typeDDouble);
  worker_multi_ddouble *access=(worker_multi_ddouble *)src_worker;
  working_assert(allocator.checkIndex(dst));
  working_assert(allocator.checkIndex(dst+len-1));
  working_assert(access->allocator.checkIndex(src));
  working_assert(access->allocator.checkIndex(src+len-1));
  //working_assert(dst>=0 && dst+len<=capacity);
  //working_assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage_[iiw__COUNT+dst+i].assign(access->storage_[iiw__COUNT+src+i]);
}


void worker_multi_ddouble::zero(const number_index store, double val)
{
  working_assert(store.dbg_worker==this);
  storage_[store].zero(val);
}

void worker_multi_ddouble::assign(const number_index store, const number_index src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src.dbg_worker==this);
  working_assert(allocator.checkIndex(store));
  working_assert(allocator.checkIndex(src));
  storage_[store] = storage_[src];
}

void worker_multi_ddouble::assign_across(const number_index store, const number_place *src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src->ptr.dbg_worker->ntype()==Type::typeDDouble);
  const worker_multi_ddouble *access=specific_cast<const worker_multi_ddouble *, const worker_multi *>(src->allocator->worker);
  working_assert((src->ptr.dbg_worker==access));
  working_assert(allocator.checkIndex(store));
  working_assert(src->allocator->checkIndex(src->ptr));
  storage_[store]=access->storage_[src->ptr];
}

void worker_multi_ddouble::chs(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].chs();
}

void worker_multi_ddouble::lshift(const number_index store, int shoft)
{
  working_assert(store.dbg_worker==this);
  storage_[store].lshift(shoft);
}

void worker_multi_ddouble::round(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].round();
}

void worker_multi_ddouble::frac(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].frac();
}

void worker_multi_ddouble::mod1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].mod1();
}

void worker_multi_ddouble::add_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].add_double(x);
}

void worker_multi_ddouble::mul_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].mul_double(x);
}

void worker_multi_ddouble::add(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].add(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_ddouble::sub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].add(-storage_[other].hi, -storage_[other].lo_);
}
void worker_multi_ddouble::rsub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].chs();
  storage_[store].add(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_ddouble::mul(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].mul(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_ddouble::sqr(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].sqr();
}

double worker_multi_ddouble::radixfloor(const number_index store1, const number_index store2)
{
  working_assert(store1.dbg_worker==this);
  working_assert(store2.dbg_worker==this);
  double rf1=storage_[store1].radixfloor();
  double rf2=storage_[store2].radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void worker_multi_ddouble::recip(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].recip();
}

void worker_multi_ddouble::sqrt(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].sqrt();
}

int worker_multi_ddouble::compare(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].compare(&storage_[other]);
}

bool worker_multi_ddouble::isequal(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].isequal(&storage_[other]);
}

bool worker_multi_ddouble::is0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].is0();
}

bool worker_multi_ddouble::isle(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].isle(&storage_[other]);
}

bool worker_multi_ddouble::isle0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isle0();
}

bool worker_multi_ddouble::isl0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isl0();
}

bool worker_multi_ddouble::isl1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isl1();
}

void worker_multi_ddouble::min(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (!storage_[store].isle(&storage_[other]))
  {
    storage_[store]=storage_[other];
  };
}

QString worker_multi_ddouble::toString(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return QString("dd(%1,%2)").arg(storage_[store].hi, 0, 'f', 16).arg(storage_[store].lo_, 0, 'g', 16);
}

int worker_multi_ddouble::toRound(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return floor(storage_[store].hi+0.5)+floor(storage_[store].lo_+0.5);
}

double worker_multi_ddouble::toDouble(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].hi;
}


number_index worker_multi_qdouble::getNumber(int index)
{
  working_assert(allocator.checkIndex(index));
  //working_assert(index>=0 && index<capacity);
  return number_index(iiw__COUNT+index, this);
}

/*void worker_multi_qdouble::getTmp12(number_index &t1, number_index &t2)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
}

void worker_multi_qdouble::getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  working_assert(t3.dbg_worker==this);
  working_assert(t4.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
  t3.index_into_worker=2;
  t4.index_into_worker=3;
}*/

worker_multi_qdouble::worker_multi_qdouble(Allocator<worker_multi> *source):
  worker_multi(Type::typeDDouble, specific_cast<Allocator<worker_multi_qdouble> *, Allocator<worker_multi> *>(source)->capacity)
{
  int oldfirst, capacity;
  source->_getFirstCapac(oldfirst, capacity);
  storage_=new dq_real[iiw__COUNT+capacity];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
        storage_[i].zero(0);
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].zero(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].hi=(double)*src_storage;
        storage_[iiw__COUNT+i].lo_=*src_storage-storage_[iiw__COUNT+i].hi;
        src_storage++;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].assign(*src_storage);
        src_storage++;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].hi=src_storage->val1;
        storage_[iiw__COUNT+i].lo_=0;
        src_storage++;
      }
    } break;
  }
}

worker_multi_qdouble::~worker_multi_qdouble()
{
  delete[] storage_;
  storage_=nullptr;
}

void worker_multi_qdouble::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  working_assert(src_worker->ntype()==Type::typeQDouble);
  worker_multi_qdouble *access=(worker_multi_qdouble *)src_worker;
  working_assert(allocator.checkIndex(dst));
  working_assert(allocator.checkIndex(dst+len-1));
  working_assert(access->allocator.checkIndex(src));
  working_assert(access->allocator.checkIndex(src+len-1));
  //working_assert(dst>=0 && dst+len<=capacity);
  //working_assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage_[iiw__COUNT+dst+i].assign(access->storage_[iiw__COUNT+src+i]);
}

void worker_multi_qdouble::zero(const number_index store, double val)
{
  working_assert(store.dbg_worker==this);
  storage_[store].zero(val);
}

void worker_multi_qdouble::assign(const number_index store, const number_index src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src.dbg_worker==this);
  working_assert(allocator.checkIndex(store));
  working_assert(allocator.checkIndex(src));
  storage_[store] = storage_[src];
}

void worker_multi_qdouble::assign_across(const number_index store, const number_place *src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src->ptr.dbg_worker->ntype()==Type::typeQDouble);
  const worker_multi_qdouble *access=specific_cast<const worker_multi_qdouble *, const worker_multi *>(src->allocator->worker);
  working_assert((src->ptr.dbg_worker==access));
  working_assert(allocator.checkIndex(store));
  working_assert(src->allocator->checkIndex(src->ptr));
  storage_[store]=access->storage_[src->ptr];
}

void worker_multi_qdouble::chs(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].chs();
}

void worker_multi_qdouble::lshift(const number_index store, int shoft)
{
  working_assert(store.dbg_worker==this);
  storage_[store].lshift(shoft);
}

void worker_multi_qdouble::round(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].round();
}

void worker_multi_qdouble::frac(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].frac();
}

void worker_multi_qdouble::mod1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].mod1();
}

void worker_multi_qdouble::add_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].add_double(x);
}

void worker_multi_qdouble::mul_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].mul_double(x);
}

void worker_multi_qdouble::add(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].add(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_qdouble::sub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].add(-storage_[other].hi, -storage_[other].lo_);
}

void worker_multi_qdouble::rsub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].chs();
  storage_[store].add(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_qdouble::mul(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].mul(storage_[other].hi, storage_[other].lo_);
}

void worker_multi_qdouble::sqr(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].sqr();
}

double worker_multi_qdouble::radixfloor(const number_index store1, const number_index store2)
{
  working_assert(store1.dbg_worker==this);
  working_assert(store2.dbg_worker==this);
  double rf1=storage_[store1].radixfloor();
  double rf2=storage_[store2].radixfloor();
  if (rf1<rf2)
    return rf2;
  return rf1;
}

void worker_multi_qdouble::recip(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].recip();
}

void worker_multi_qdouble::sqrt(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].sqrt();
}

int worker_multi_qdouble::compare(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].compare(&storage_[other]);
}

bool worker_multi_qdouble::isequal(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].isequal(&storage_[other]);
}

bool worker_multi_qdouble::is0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].is0();
}

bool worker_multi_qdouble::isle(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].isle(&storage_[other]);
}

bool worker_multi_qdouble::isle0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isle0();
}

bool worker_multi_qdouble::isl0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isl0();
}

bool worker_multi_qdouble::isl1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].isl1();
}

void worker_multi_qdouble::min(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (!storage_[store].isle(&storage_[other]))
  {
    storage_[store]=storage_[other];
  };
}

QString worker_multi_qdouble::toString(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return QString("dd(%1,%2)").arg(storage_[store].hi, 0, 'f', 16).arg(storage_[store].lo_, 0, 'g', 16);
}

int worker_multi_qdouble::toRound(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return floor(storage_[store].hi+0.5)+floor(storage_[store].lo_+0.5);
}

double worker_multi_qdouble::toDouble(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].hi;
}





worker_multi_real642::worker_multi_real642(Allocator<worker_multi> *source):
  worker_multi(Type::typeReal642, specific_cast<Allocator<worker_multi_real642> *, Allocator<worker_multi> *>(source)->capacity)
{
  int oldfirst, capacity;
  source->_getFirstCapac(oldfirst, capacity);
  storage_=new real642[iiw__COUNT+capacity];
  switch (source->worker->ntype())
  {
    case Type::typeEmpty:
    {
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=0;
        storage_[iiw__COUNT+i].val2=0;
      }
    } break;
    case Type::typeDouble:
    {
      const double *src_storage=specific_cast<worker_multi_double *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      //TODO: memmove but we'll see later
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=*src_storage;
        storage_[iiw__COUNT+i].val2=*src_storage;
        src_storage++;
      }
    } break;
    case Type::typeFloat128:
    {
      const __float128 *src_storage=specific_cast<worker_multi_float128 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=*src_storage;
        storage_[iiw__COUNT+i].val2=*src_storage;
        src_storage++;
      }
    } break;
    case Type::typeDDouble:
    {
      const dd_real *src_storage=specific_cast<worker_multi_ddouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=src_storage->hi;
        storage_[iiw__COUNT+i].val2=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case Type::typeQDouble:
    {
      const dq_real *src_storage=specific_cast<worker_multi_qdouble *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=src_storage->hi;
        storage_[iiw__COUNT+i].val2=src_storage->hi+(__float128)src_storage->lo_;
        src_storage++;
      }
    } break;
    case Type::typeReal642:
    {
      const real642 *src_storage=specific_cast<worker_multi_real642 *, worker_multi *>(source->worker)->_getStorage_()+iiw__COUNT;
      src_storage+=oldfirst;
      for (int i=0; i<capacity; i++)
      {
        storage_[iiw__COUNT+i].val1=src_storage->val1;
        storage_[iiw__COUNT+i].val2=src_storage->val2;
        src_storage++;
      }
    } break;
  }
}

worker_multi_real642::~worker_multi_real642()
{
  delete[] storage_;
  storage_=nullptr;
}

number_index worker_multi_real642::getNumber(int index)
{
  working_assert(allocator.checkIndex(index));
  //working_assert(index>=0 && index<capacity);
  return number_index(this, iiw__COUNT+index);
}

/*void worker_multi_real642::getTmp12(number_index &t1, number_index &t2)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
}

void worker_multi_real642::getTmp1234(number_index &t1, number_index &t2, number_index &t3, number_index &t4)
{
  working_assert(t1.dbg_worker==this);
  working_assert(t2.dbg_worker==this);
  working_assert(t3.dbg_worker==this);
  working_assert(t4.dbg_worker==this);
  t1.index_into_worker=0;
  t2.index_into_worker=1;
  t3.index_into_worker=2;
  t4.index_into_worker=3;
}*/

void worker_multi_real642::assign_block(int dst, worker_multi *src_worker, int src, int len)
{
  working_assert(src_worker->ntype()==Type::typeReal642);
  worker_multi_real642 *access=(worker_multi_real642 *)src_worker;
  working_assert(allocator.checkIndex(dst));
  working_assert(allocator.checkIndex(dst+len-1));
  working_assert(access->allocator.checkIndex(src));
  working_assert(access->allocator.checkIndex(src+len-1));
  //working_assert(dst>=0 && dst+len<=capacity);
  //working_assert(src>=0 && src+len<=access->capacity);

  for (int i=0; i<len; i++)
    storage_[iiw__COUNT+dst+i]=access->storage_[iiw__COUNT+src+i];
}

void worker_multi_real642::zero(const number_index store, double val)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1=val;
  storage_[store].val2=val;
}

void worker_multi_real642::assign(const number_index store, const number_index src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src.dbg_worker==this);
  working_assert(allocator.checkIndex(store));
  working_assert(allocator.checkIndex(src));
  storage_[store]=storage_[src];
}

void worker_multi_real642::assign_across(const number_index store, const number_place *src)
{
  working_assert(store.dbg_worker==this);
  working_assert(src->ptr.dbg_worker->ntype()==Type::typeReal642);
  const worker_multi_real642 *access=specific_cast<const worker_multi_real642 *, const worker_multi *>(src->allocator->worker);
  working_assert((src->ptr.dbg_worker==access));
  working_assert(allocator.checkIndex(store));
  working_assert(src->allocator->checkIndex(src->ptr));
  storage_[store]=access->storage_[src->ptr];
}

void worker_multi_real642::chs(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 = -storage_[store].val1;
  storage_[store].val2 = -storage_[store].val2;
}

void worker_multi_real642::lshift(const number_index store, int shoft)
{
  working_assert(store.dbg_worker==this);
  //*store.asf128=ldexp(*store.asf128, shoft);
  //could try ldexp only on upper half but...
  storage_[store].val1=ldexp(storage_[store].val1, shoft);
  if (storage_[store].val2 != 0) //don't play with 0's exponent or you get NaN
    ((uint16_t *)&storage_[store].val2)[7]+=shoft;
}

void worker_multi_real642::round(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1=std::round(storage_[store].val1);
  __float128 remainder=storage_[store].val2;
  storage_[store].val2=0;
  for (;;)
  {
    double some=std::round((double)remainder);
    if (some==0)
      break;
    storage_[store].val2+=some;
    remainder-=some;
  }
}

void worker_multi_real642::frac(const number_index store)
{
  working_assert(store.dbg_worker==this);
  if (storage_[store].val1<0)
    storage_[store].val1 -= ceil(storage_[store].val1);
  else
    storage_[store].val1 -= floor(storage_[store].val1);
  if (storage_[store].val2<0)
  {
    double some=ceil((double)storage_[store].val2);
    while (some)
    {
      storage_[store].val2 -= some;
      some=ceil((double)storage_[store].val2);
    }
  }
  else
  {
    double some=floor((double)storage_[store].val2);
    while (some)
    {
      storage_[store].val2 -= some;
      some=floor((double)storage_[store].val2);
    }
  }
}

void worker_multi_real642::mod1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 -= floor(storage_[store].val1);
  double some=floor((double)storage_[store].val2);
  while (some)
  {
    storage_[store].val2 -= some;
    some=floor((double)storage_[store].val2);
  }
}

void worker_multi_real642::add_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 += x;
  storage_[store].val2 += x;
}

void worker_multi_real642::mul_double(const number_index store, double x)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 *= x;
  storage_[store].val2 *= x;
}

void worker_multi_real642::add(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].val1 += storage_[other].val1;
  storage_[store].val2 += storage_[other].val2;
}

void worker_multi_real642::sub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].val1 -= storage_[other].val1;
  storage_[store].val2 -= storage_[other].val2;
}

void worker_multi_real642::rsub(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].val1 = storage_[other].val1-storage_[store].val1;
  storage_[store].val2 = storage_[other].val2-storage_[store].val2;
}

void worker_multi_real642::mul(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  storage_[store].val1 *= storage_[other].val1;
  storage_[store].val2 *= storage_[other].val2;
}

void worker_multi_real642::sqr(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 *= storage_[store].val1;
  storage_[store].val2 *= storage_[store].val2;
}

double worker_multi_real642::radixfloor(const number_index store1, const number_index store2)
{
  working_assert(store1.dbg_worker==this);
  working_assert(store2.dbg_worker==this);
  int ilog1=std::ilogb(storage_[store1].val1);
  int ilog2=std::ilogb(storage_[store2].val1);
  if (ilog1<ilog2)
    ilog1=ilog2;
  return ldexp(1, ilog1);
}

void worker_multi_real642::recip(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 = 1/storage_[store].val1;
  storage_[store].val2 = 1/storage_[store].val2;
}

void worker_multi_real642::sqrt(const number_index store)
{
  working_assert(store.dbg_worker==this);
  storage_[store].val1 = std::sqrt(storage_[store].val1);
  double x=1/std::sqrt((double)storage_[store].val2);
  double ax=x*(double)storage_[store].val2;
  storage_[store].val2 = ax + double(storage_[store].val2 - __float128(ax)*__float128(ax))*x/2;
}

int worker_multi_real642::compare(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store].val1 == storage_[other].val1)
    return 0;
  else if (storage_[store].val1 < storage_[other].val1)
    return -1;
  else
    return +1;
}

bool worker_multi_real642::isequal(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].val1 == storage_[other].val1;
}

bool worker_multi_real642::is0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].val1 == 0;
}

bool worker_multi_real642::isle(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  return storage_[store].val1 <= storage_[other].val1;
}

bool worker_multi_real642::isle0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].val1 <= 0;
}

bool worker_multi_real642::isl0(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].val1 < 0;
}

bool worker_multi_real642::isl1(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store].val1 < 1;
}

void worker_multi_real642::min(const number_index store, const number_index other)
{
  working_assert(store.dbg_worker==this);
  working_assert(other.dbg_worker==this);
  if (storage_[store].val1>storage_[other].val1)
  {
    storage_[store]=storage_[other];
  }
}

QString worker_multi_real642::toString(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return QString::number(storage_[store.index_into_worker].val1, 'f', 16);
}

int worker_multi_real642::toRound(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return qRound(storage_[store.index_into_worker].val1);
}

double worker_multi_real642::toDouble(const number_index store)
{
  working_assert(store.dbg_worker==this);
  return storage_[store.index_into_worker].val1;
}

#endif


/*
number::number(number_pointer *src):
  worker(src->worker), allocator(src->allocator),
  free_storage_on_destroy(src->free_storage_on_destroy),
  ptr(src->ptr)
{
  src->ptr.asf64=nullptr;
}*/

#if 0
template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::zero(double v)
{
  worker->zero(_ptr, v);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::assign(const number_index src)
{
  worker->assign(_ptr, src);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::assign_across(const number_place<WORKER_MULTI> *src)
{
  worker->assign_across(_ptr, src);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::assign_across(const number_place<WORKER_MULTI> &&src)
{
  worker->assign_across(_ptr, &src);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::lshift(int shoft)
{
  worker->lshift(_ptr, shoft);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::add(const number_index other)
{
  worker->add(_ptr, other);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::chs()
{
  worker->chs(_ptr);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::sub(const number_index other)
{
  worker->sub(_ptr, other);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::sqr()
{
  worker->sqr(_ptr);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::mul(const number_index other)
{
  worker->mul(_ptr, other);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::recip()
{
  worker->recip(_ptr);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::sqrt()
{
  worker->sqrt(_ptr);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::add_double(double x)
{
  worker->add_double(_ptr, x);
}

template<class WORKER_MULTI>
bool number_instance<WORKER_MULTI>::is0() const
{
  return worker->is0(_ptr);
}

template <class WORKER_MULTI>
void number_instance<WORKER_MULTI>::min(const number_index other)
{
  worker->min(_ptr, other);
}

template <class WORKER_MULTI>
int number_instance<WORKER_MULTI>::toRound() const
{
  return worker->toRound(_ptr);
}

template <class WORKER_MULTI>
double number_instance<WORKER_MULTI>::toDouble() const
{
  return worker->toDouble(_ptr);
}


template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::zero(double r, double i)
{
  re.zero(r);
  im.zero(i);
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::assign(const complex_index other)
{
  //working_assert(allocator!=nullptr);
  //working_assert(allocator->worker!=nullptr);
  //working_assert(re.dbg_worker!=nullptr);
  //working_assert(im.dbg_worker!=nullptr);
  working_assert(other.dbg_worker!=nullptr);
  re.assign(other.getre());
  im.assign(other.getim());
}

template<class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::assign_across(const complex_place<WORKER_MULTI> *other)
{
  //working_assert(allocator!=nullptr);
  //working_assert(allocator->worker!=nullptr);
  working_assert(this!=nullptr);
  working_assert(other!=nullptr);
  if (re.worker==other->re.dbg_worker)
    nop();
  re.assign_across(other->re_place());
  im.assign_across(other->im_place());
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::lshift(int shoft)
{
  re.lshift(shoft);
  im.lshift(shoft);
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::add(const complex_index other)
{
  re.add(other.getre());
  im.add(other.getim());
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::chs()
{
  re.chs();
  im.chs();
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::sub(const complex_index other)
{
  re.sub(other.getre());
  im.sub(other.getim());
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::rsub(const complex_index other)
{
  re.worker->rsub(re.getptr(), other.getre());
  im.worker->rsub(im.getptr(), other.getim());
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::mul(const complex_index other)
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //working_assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //r:=r1*r2-i1*i2
  //i:=r1*i2+i1*r2
  re.worker->assign(tmp1, re.getptr());
  re.worker->mul(tmp1, other.getre());
  re.worker->assign(tmp2, im.getptr());
  re.worker->mul(tmp2, other.getim());
  re.worker->sub(tmp1, tmp2); //r1*r2-i1*i2
  re.worker->mul(re.getptr(), other.getim());
  re.worker->mul(im.getptr(), other.getre());
  re.worker->add(im.getptr(), re.getptr()); //i1*r2+r1*i2
  re.worker->assign(re.getptr(), tmp1);
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::sqr()
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  //if ((tmp1.store==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //working_assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
#if 0
  //r:=r*r-i*i
  //i=2*r*i
  worker->assign(tmp1, im);
  worker->sqr(tmp1);
  worker->mul(im, re);
  worker->lshift(im, 1);
  worker->sqr(re);
  worker->sub(re, tmp1);
#else
  //r:=(r+i)*(r-i)
  //i=2*r*i
  //should be more precise + fewer mul
  re.worker->assign(tmp1, re.getptr());
  re.worker->sub(tmp1, im.getptr());    //r-i
  re.worker->assign(tmp2, im.getptr());
  re.worker->mul(im.getptr(), re.getptr());
  re.worker->lshift(im.getptr(), 1);    //2*r*i
  re.worker->add(re.getptr(), tmp2);    //r+i
  re.worker->mul(re.getptr(), tmp1);
#endif
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::recip()
{
  getMag_tmp();
  recip_prepared();
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::recip_prepared()
{ // 1/(re+i*im) = (re-i*im)/((re+i*im)*(re-i*im)) = (re-i*im)/(re*re+im*im)
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->recip(tmp1);
  re.worker->mul(re.getptr(), tmp1);
  re.worker->chs(im.getptr());
  re.worker->mul(im.getptr(), tmp1);
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::sqrt()
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp1);
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2); //re*re+im*im
  re.worker->sqrt(tmp1);
  if (!re.worker->isle0(re.getptr()))
  {
    re.worker->add(tmp1, re.getptr());
    re.worker->lshift(tmp1, -1);
    re.worker->sqrt(tmp1); //t1=sqrt((sqrt(re*re+im*im)+re)/2);
    re.worker->assign(re.getptr(), tmp1); //re=t1
    /*if (t1==0)
      *res_im=0;
    else*/
    re.worker->lshift(tmp1, 1);
    re.worker->recip(tmp1);
    re.worker->mul(im.getptr(), tmp1); //im=im/(2*t1);
  }
  else
  {
    re.worker->sub(tmp1, re.getptr());
    re.worker->lshift(tmp1, -1);
    re.worker->sqrt(tmp1); //t1=sqrt((sqrt(re*re+im*im)-re)/2);
    if (re.worker->isle0(tmp1)) //t1==0
    {
      re.worker->zero(re.getptr());
      re.worker->zero(im.getptr());
    }
    else
    {
      re.worker->assign(re.getptr(), im.getptr());
      re.worker->assign(im.getptr(), tmp1); //new im=t1
      re.worker->lshift(tmp1, 1);
      re.worker->recip(tmp1);
      re.worker->mul(re.getptr(), tmp1); //re=old im/(2*t1);
      if (re.worker->isl0(re.getptr()))
      {
        re.worker->chs(re.getptr());
        re.worker->chs(im.getptr());
      };
    }
  };
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::root_approx(int n)
{
  if (is0())
    return;
  double this_re=re.toDouble();
  double this_im=im.toDouble();
  double r_abs=std::exp(std::log(this_re*this_re+this_im*this_im)/(2*n));
  double r_phi=std::atan2(this_im, this_re)/n;
  double r_re=r_abs*cos(r_phi);
  double r_im=r_abs*sin(r_phi);
  zero(r_re, r_im);
}

template <class WORKER_MULTI>
void complex_instance<WORKER_MULTI>::sign()
{
  number_index donotuse(re.worker), tmp2(re.worker);
  re.worker->getTmp12(donotuse, tmp2);
  re.worker->assign(tmp2, getMag_tmp());
  re.worker->sqrt(tmp2);
  re.worker->recip(tmp2);
  re.worker->mul(re.getptr(), tmp2);
  re.worker->mul(im.getptr(), tmp2);
}

template <class WORKER_MULTI>
double complex_instance<WORKER_MULTI>::getMag_double() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //working_assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->sqr(tmp1);
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return re.worker->toDouble(tmp1);
}

template <class WORKER_MULTI>
const number_index complex_instance<WORKER_MULTI>::getMag_tmp() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //working_assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->sqr(tmp1);
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return tmp1;
}

template <class WORKER_MULTI>
const number_index complex_instance<WORKER_MULTI>::getMag1_tmp() const
{
  //if ((tmp1_s==nullptr) || (tmp2.store==nullptr))
  //  dbgPoint();
  //working_assert((tmp1.store!=nullptr) && (tmp2.store!=nullptr));
  //TODO: merge with mag_cmp_1
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->sqr(tmp1);
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  re.worker->add_double(tmp1, -1);
  return tmp1;
}

template <class WORKER_MULTI>
const number_index complex_instance<WORKER_MULTI>::getDist1_tmp() const
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->add_double(tmp1, -1);
  re.worker->sqr(tmp1);
  re.worker->assign(tmp2, im.getptr());
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return tmp1;
}

template <class WORKER_MULTI>
const number_index complex_instance<WORKER_MULTI>::mulreT_tmp(const complex_index other) const //Re(this*conjugate(other))
{ //re*ore+im*oim
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->assign(tmp2, im.getptr());
  re.worker->mul(tmp1, other.getre());
  re.worker->mul(tmp2, other.getim());
  re.worker->add(tmp1, tmp2);
  return tmp1;
}

template <class WORKER_MULTI>
double complex_instance<WORKER_MULTI>::dist2_double(const complex_index other) const
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->assign(tmp2, im.getptr());
  re.worker->sub(tmp1, other.getre());
  re.worker->sub(tmp2, other.getim()->im.getptr());
  re.worker->sqr(tmp1);
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return re.worker->toDouble(tmp1);
}

template <class WORKER_MULTI>
const number_index complex_instance<WORKER_MULTI>::dist2_tmp(const complex_instance<WORKER_MULTI> *other) const
{
  number_index tmp1(re.worker), tmp2(re.worker);
  re.worker->getTmp12(tmp1, tmp2);
  re.worker->assign(tmp1, re.getptr());
  re.worker->assign(tmp2, im.getptr());
  re.worker->sub(tmp1, other->re.getptr());
  re.worker->sub(tmp2, other->im.getptr());
  re.worker->sqr(tmp1);
  re.worker->sqr(tmp2);
  re.worker->add(tmp1, tmp2);
  return tmp1;
}

template <class WORKER_MULTI>
bool complex_instance<WORKER_MULTI>::isequal(const complex_instance<WORKER_MULTI> *other) const
{
  return re.worker->isequal(re.getptr(), other->re.getptr()) &&
         re.worker->isequal(im.getptr(), other->im.getptr());
}

template <class WORKER_MULTI>
bool complex_instance<WORKER_MULTI>::is0() const
{
  return re.is0() && im.is0();
}

template <class WORKER_MULTI>
int complex_instance<WORKER_MULTI>::mag_cmp_1()
{
  number_index tmp1(re.worker), tmp2(re.worker), tmp3(re.worker), tmp4(re.worker);
  re.worker->getTmp1234(tmp1, tmp2, tmp3, tmp4);
  //first reduce fz_re, fz_im to 2nd octant:
  //  fz.re>=0, fz.im>=0, fz.im>=fz.re
  re.worker->assign(tmp1, re.getptr());
  re.worker->assign(tmp2, im.getptr());
  if (re.worker->isle0(tmp1)) re.worker->chs(tmp1);
  if (re.worker->isle0(tmp2)) re.worker->chs(tmp2);
  if (!re.worker->isle(tmp1, tmp2))
  {
    re.worker->assign(tmp3, tmp1);
    re.worker->assign(tmp1, tmp2);
    re.worker->assign(tmp2, tmp3);
  };
  //im>=re, now check re^2+im^2>1
  if ((re.worker->toDouble(tmp1)>=0.71) || //im>=re>=0.71 -> mag>=1.0082
      (re.worker->toDouble(tmp2)>=1.01))   //im>=1.01 -> mag>=1.01
  {
    return 1;
  }
  else if (re.worker->toDouble(tmp2)<=0.70) //0.70>=im>=re -> mag<=0.98
    return -1;
  else
  { //   re*re+im*im>1
    //   gets inaccurate for re~0, im~1
    //   im*im>1-re*re ; im>sqrt(1-re*re) ; im-1>sqrt(1-re*re)-1
    //>> x=Sqrt(1-r*r)-1
    //   could try x:=((-r2/16-1/8)*r2-1/2)*r2; //-r^6/16-r^4/8-r^2/2
    //   but one cycle of Newton should work with any precision
    //   (x+1)*(x+1)=(1-r*r)
    //   x*x+2*x+r*r=0  f'=2*x+2
    //   x2=x-(x*x+2*x+r*r)/(2*x+2)
    //>> x2=(x*x-r*r)/(2*x+2)
    //   x2=(x-1)(r*r-x*x)/(1-x*x)/2
    re.worker->sqr(tmp1); //r^2
    re.worker->chs(tmp1);
    re.worker->assign(tmp3, tmp1);
    re.worker->add_double(tmp3, 1);
    re.worker->sqrt(tmp3);              //sqrt(1-r*r) = x+1
    re.worker->assign(tmp4, tmp3);
    re.worker->add_double(tmp4, -1);        //x
    re.worker->sqr(tmp4);
    re.worker->add(tmp1, tmp4);    //x*x-r*r
    re.worker->recip(tmp3); //-1..-0.51 -> -1..-1.96
    re.worker->mul(tmp1, tmp3);
    re.worker->lshift(tmp1, -1);           //(x*x-r*r)/(x+1)/2 = better x
    re.worker->add_double(tmp2, -1);
    return re.worker->compare(tmp2, tmp1); //im-1 vs limit
  };
}

template <class WORKER_MULTI>
QString complex_instance<WORKER_MULTI>::toString()
{
  return re.toString()+" +i* "+im.toString();
}
#endif //non-double workers removed

#endif // all workers removed

double sqr_double(double x)
{
  return x*x;
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

template class number<double>;
//template void number<double>::assign_across<double>(const number<double> *src);
template class complex<double>;
template class complex<number_a *>;
template void complex<double>::assign_across<number_a *>(const complex<number_a *> *src);
template void complex<double>::assign_across<double>(const complex<double> *src);
template void complex<number_a *>::assign_across<number_a *>(const complex<number_a *> *src);
#if !NUMBER_DOUBLE_ONLY
template class complex<__float128>;
template void complex<__float128>::assign_across<number_a *>(const complex<number_a *> *src);
template void complex<__float128>::assign_across<__float128>(const complex<__float128> *src);
template class complex<dd_real>;
template void complex<dd_real>::assign_across<number_a *>(const complex<number_a *> *src);
template void complex<dd_real>::assign_across<dd_real>(const complex<dd_real> *src);
template class complex<real642>;
template void complex<real642>::assign_across<number_a *>(const complex<number_a *> *src);
template void complex<real642>::assign_across<real642>(const complex<real642> *src);

#endif
/*template class complex<worker_multi_float128>;
template class complex<worker_multi_ddouble>;
template class complex<worker_multi_qdouble>;
template class complex<worker_multi_real642>;*/

} // namespace MandelMath
