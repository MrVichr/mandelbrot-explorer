#include "quad_double.hpp"
#include "double_double.hpp"
#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

static void __nop()
{

}

void qd_dbgPoint()
{
  __nop();
}

#if defined(_WIN32) && !defined(__GNUC__)
#ifdef __BORLANDC__
#else
  /* Win 32 MSVC */
  #include <float.h>
#endif
#endif

namespace MandelMath {

// https://www.davidhbailey.com/dhbsoftware/
// https://www.davidhbailey.com/dhbsoftware/qd-2.3.23.tar.gz
// https://www.davidhbailey.com/dhbsoftware/dqfun-v01.tar.gz

#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

//log((2^53-1)/pi)/log(2) = 51.xxx
//    round(pi*2^51)=7074237752028440
//log((2^53-1)/(pi-7074237752028440/2^51))/log(2) = 105.xxx
//    round((pi-7074237752028440/2^51)*2^105)=4967757600021511
//log((2^53-1)/(pi-7074237752028440/2^51-4967757600021511/2^105))/log(2) = 161.xxx+xxx*I
//    round((pi-7074237752028440/2^51-4967757600021511/2^105)*2^161)=-8753721960665020
//log((2^53-1)/(pi-7074237752028440/2^51-4967757600021511/2^105+8753721960665020/2^161))/log(2) = 215.xxx
//    round((pi-7074237752028440/2^51-4967757600021511/2^105+8753721960665020/2^161)*2^215)=5857755168774013
static const qd_real M_PIqd(7074237752028440/2251799813685248.0, // 3.141592653589793116e+00, *2^51
                            4967757600021511/40564819207303340847894502572032.0, // 1.224646799147353207e-16, *2^105
                            -2188430490166255/730750818665451459101842416358141509827966271488.0, //-2.994769809718339666e-33, *2^159
                            5857755168774013/52656145834278593348959013841835216159447547700274555627155488768.0); //1.112454220863365282e-49); *2^215

/*********** Basic Functions ************/
/* Computes fl(a+b) and err(a+b).  Assumes |a| >= |b|. */
inline void qd_real::quick_nth_sum01()//double a, double b)
{
  double a=x[0];
  x[0] = x[0] + x[1];
  x[1] = x[1] + (a - x[0]);
  //lo_ = b - (hi - a);

}

inline void qd_real::quick_nth_sum12()
{
  double a=x[1];
  x[1] = x[1] + x[2];
  x[2] = x[2] + (a - x[1]);
}

inline void qd_real::quick_nth_sum23()
{
  double a=x[2];
  x[2] = x[2] + x[3];
  x[3] = x[3] + (a - x[2]);
}

inline void qd_real::quick_nth_sum3x(double &x4)
{
  double a=x[3];
  x[3] = x[3] + x4;
  x4 = x4 + (a - x[3]);
}

/* Computes fl(a-b) and err(a-b).  Assumes |a| >= |b| */
inline double quick_two_diff(double a, double b, double &err) {
  double s = a - b;
  err = (a - s) - b;
  return s;
}

/* Computes fl(a+b) and err(a+b).  */
inline void qd_real::two_sum_(double &a, double &b)
{
  double oa=a;
  a = oa + b;
  double bb = a - oa;
  b = (oa - (a - bb)) + (b - bb);
}

inline void qd_real::two_sum0(double &b)
{
  two_sum_(x[0], b);
}

inline void qd_real::two_sum1(double &b)
{
  two_sum_(x[1], b);
}

inline void qd_real::two_sum2(double &b)
{
  two_sum_(x[2], b);
}

inline void qd_real::two_sum3(double &b)
{
  two_sum_(x[3], b);
}

inline void qd_real::three_sum(double &a, double &b, double &c)
{
  qd_real::two_sum_(a, b);
  qd_real::two_sum_(a, c);
  qd_real::two_sum_(b, c);
}

inline void qd_real::three_sum2(double &a, double &b, double &c)
{
  two_sum_(a, b);
  two_sum_(a, c);
  b = b + c;
}

/* Computes fl(a-b) and err(a-b).  */
inline double two_diff(double a, double b, double &err) {
  double s = a - b;
  double bb = s - a;
  err = (a - (s - bb)) - (b + bb);
  return s;
}

/* Computes high word and lo word of a */
/*inline*/ void qd_real::split(double a, double &hi, double &lo)
{
  double temp;
  if (a > _QD_SPLIT_THRESH || a < -_QD_SPLIT_THRESH) {
    a *= 3.7252902984619140625e-09;  // 2^-28
    temp = _QD_SPLITTER * a;
    hi = temp - (temp - a);
    lo = a - hi;
    hi *= 268435456.0;          // 2^28
    lo *= 268435456.0;          // 2^28
  } else {
    temp = _QD_SPLITTER * a; //a=a1*2^27+a2*1  temp=a1*2^54+(a1+a2)*2^27 +(a2)*1 -a2
    hi = temp - (temp - a);  //temp-a=a1*2^54+a2*2^27 -a2    temp-(temp-a)=a1*2^27
    lo = a - hi;
  }
}

/* Computes fl(a*b) and err(a*b). */
inline double qd_real::two_prod(double &a, double b)
{
  double a_hi, a_lo;
  double b_hi, b_lo;
  split(a, a_hi, a_lo);
  split(b, b_hi, b_lo);
  a*=b;
  return ((a_hi * b_hi - a) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
}

/* Computes fl(a*a) and err(a*a).  Faster than the above method. */
inline double qd_real::two_sqr(double &a)
{
  double a_hi, a_lo;
  split(a, a_hi, a_lo);
  a *= a;
  return ((a_hi * a_hi - a) + 2.0 * a_hi * a_lo) + a_lo * a_lo;
}

//hack if you don't want to mess with FPU control word... but is a nightmare
#define MAKE_DOUBLE( tmp, val ) ( tmp=(val), *(volatile double *)&tmp )

double qd_real::eps2() const
{
  return 3.60664e-130; // 2^-(2*(53+54+54+54)) rounded up
}

double qd_real::eps234() const
{
  return 8.276113e-98; // eps2^(3/4)
}

qd_real &qd_real::lshift(int exp)
{
  //lolwut "On many implementations, std::ldexp is less efficient than multiplication or division
  //        by a power of two using arithmetic operators"
  //should I just adjust exponent in *(int *)&hi or what
  x[0]=ldexp(x[0], exp);
  x[1]=ldexp(x[1], exp);
  x[2]=ldexp(x[2], exp);
  x[3]=ldexp(x[3], exp);
  return *this;
}

/*void qd_real::quick_renorm(double &x4)
{
  double t0, t1, t2, t3;
  double s;
  s  = qd::quick_two_sum(c3, c4, t3);
  s  = qd::quick_two_sum(c2, s , t2);
  s  = qd::quick_two_sum(c1, s , t1);
  c0 = qd::quick_two_sum(c0, s , t0);

  s  = qd::quick_two_sum(t2, t3, t2);
  s  = qd::quick_two_sum(t1, s , t1);
  c1 = qd::quick_two_sum(t0, s , t0);

  s  = qd::quick_two_sum(t1, t2, t1);
  c2 = qd::quick_two_sum(t0, s , t0);

  c3 = t0 + t1;
}*/

inline void qd_real::renorm()
{
  //if (QD_ISINF(x[0])) return;

  quick_nth_sum23();
  quick_nth_sum12();
  quick_nth_sum01();

  if (x[1] != 0.0) {
    quick_nth_sum12();
    if (x[2] != 0.0)
        quick_nth_sum23();
    else
    {
        x[2]=x[3];
        quick_nth_sum12();
        x[3]=0;
    }
  } else {
    x[1]=x[2];
    quick_nth_sum01();
    if (x[1] != 0.0)
    {
        x[2]=x[3];
        quick_nth_sum12();
        x[3]=0;
    }
    else
    {
        x[1]=x[3];
        quick_nth_sum01();
        x[2]=0;
        x[3]=0;
    }
  }
}

inline void qd_real::renorm(double &x4)
{
  //double s0, s1, s2 = 0.0, s3 = 0.0;

  //if (QD_ISINF(x[0])) return;

  quick_nth_sum3x(x4);
  quick_nth_sum23();
  quick_nth_sum12();
  quick_nth_sum01();

  //s0 = x[0];
  //s1 = x[1];

  if (x[1] != 0.0) {
    quick_nth_sum12();
    if (x[2] != 0.0) {
        quick_nth_sum23();
        if (x[3] != 0.0)
            x[3] += x4;
        else
        {
            x[3]=x4;
            quick_nth_sum23();
        }
    } else {
        x[2]=x[3];
        quick_nth_sum12();
        if (x[2] != 0.0)
        {
            x[3]=x4;
            quick_nth_sum23();
        }
        else
        {
            x[2]=x4;
            quick_nth_sum12();
            x[3]=0;
        }
    }
  } else {
    x[1]=x[2];
    quick_nth_sum01();
    if (x[1] != 0.0) {
        x[2]=x[3];
        quick_nth_sum12();
        if (x[2] != 0.0)
        {
            x[3]=x4;
            quick_nth_sum23();
        }
        else
        {
            x[2]=x4;
            quick_nth_sum12();
            x[3]=0;
        }
    } else {
        x[1]=x[3];
        quick_nth_sum01();
        if (x[1] != 0.0)
        {
            x[2]=x4;
            quick_nth_sum12();
            x[3]=0;
        }
        else
        {
            x[1]=x4;
            quick_nth_sum01();
            x[2]=0;
            x[3]=0;
        }
    }
  }
}

qd_real &qd_real::add_double(double b)
{
  two_sum0(b);
  two_sum1(b);
  two_sum2(b);
  two_sum3(b);

  renorm(b);
  return *this;
}

qd_real &qd_real::mul_double(double b)
{
  double s4;

  double x5=x[3]*b;
  s4 = two_prod(x[2], b);
  double x4=x[2];
  x[3] = two_prod(x[1], b);
  x[2] = two_prod(x[0], b);

  two_sum_(x[1], x[2]);

  three_sum(x[2], x[3], x4);

  three_sum2(x[3], s4, x5);

  s4 += x4;

  renorm(s4);
  return *this;
}

qd_real &qd_real::add(qd_real const &other)
{
  //sloppy add
  qd_real t(other);

  two_sum0(t.x[0]);
  two_sum1(t.x[1]);
  two_sum2(t.x[2]);
  two_sum3(t.x[3]);

  two_sum1(t.x[0]);
  three_sum(x[2], t.x[0], t.x[1]);
  three_sum2(x[3], t.x[0], t.x[2]);
  t.x[0] = t.x[0] + t.x[1] + t.x[3];

  renorm(t.x[0]);
  return *this;
}

qd_real &qd_real::sub(qd_real const &other)
{
  qd_real t(-other.x[0], -other.x[1], -other.x[2], -other.x[3]);

  two_sum0(t.x[0]);
  two_sum1(t.x[1]);
  two_sum2(t.x[2]);
  two_sum3(t.x[3]);

  two_sum1(t.x[0]);
  three_sum(x[2], t.x[0], t.x[1]);
  three_sum2(x[3], t.x[0], t.x[2]);
  t.x[0] = t.x[0] + t.x[1] + t.x[3];

  renorm(t.x[0]);
  return *this;
}

#if 0 //this won't be faster than Bailey's
void split3(double x, double &h, double &m, double &l)
{
  double temp = 0x400000001 * x;
  double r = temp - (temp - x);  //temp-a=a1*2^54+a2*2^27 -a2    temp-(temp-a)=a1*2^27
  l = (x - r)*0x400000000;
  temp = 0x20001 * r;
  h=temp-(temp-r);
  m=(r-h)*0x20000;
}

void mul_fft(qd_real &result, qd_real const &other)
{
  /*
convolution 12+12->24:
  [ a (4) ]  [ b (4) ]  [ c (4) ]
  [ d (4) ]  [ e (4) ]  [ f (4) ]
->
  [ g (8) ]  [ h (8) ]  [ i (8) ]  [ j (8) ]  [ k (8) ]
  ghijk=[5x5 (int)] ( [3x5 int] abc . [3x5 int] def )
    probably don't need k, maybe j

multiplication [4] x [4] -> [7]
classic: 16 mul + renorm      abcd*efgh=ae,af+be,ag+bf+ce,ah+bg+cf+de,bh+cg+df,ch+dg,dh
karats 1x: [ab] x [cd] =

at 7 points: inf,0,1,-1,i,-1,2
{{0,  0,  0,  0,  0,  0,  1},
 {1,  0,  0,  0,  0,  0,  0},
 {1,  1,  1,  1,  1,  1,  1},
 {1, -1,  1, -1,  1, -1,  1},
 {1,  i, -1, -i,  1,  i, -1},
 {1, -i, -1,  i,  1, -i, -1},
 {1,  2,  4,  8, 16, 32, 64}}

{{0,0,0,0,0,0,1},{1,0,0,0,0,0,0},{1,1,1,1,1,1,1},{1,-1,1,-1,1,-1,1},{1,i,-1,-i,1,i,-1},{1,-i,-1,i,1,-i,-1},{1,2,4,8,16,32,64}}^-1


{{a^0,a^0,a^0,a^0,a^0,a^0,a^0,a^0},
 {a^1,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^2,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^3,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^4,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^5,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^6,a^2,a^3,a^4,a^5,a^6,a^7,a^8},
 {a^7,a^2,a^3,a^4,a^5,a^6,a^7,a^8},

*/
  //might as well do it all in int64 instead of doubles...
  qd_real a, b, c, d, e, f;
  split3(result.x[0], a.x[0], a.x[1], a.x[2]);
  split3(result.x[1], a.x[3], b.x[0], b.x[1]);
  split3(result.x[2], b.x[2], b.x[3], c.x[0]);
  split3(result.x[3], c.x[1], c.x[2], c.x[3]);
  split3(other.x[0], d.x[0], d.x[1], d.x[2]);
  split3(other.x[1], d.x[3], e.x[0], e.x[1]);
  split3(other.x[2], e.x[2], e.x[3], f.x[0]);
  split3(other.x[3], f.x[1], f.x[2], f.x[3]);
  //eval at 0, 1, -1, 2, inf
  qd_real p0(a);
  qd_real p1(a.x[0]+  b.x[0]+  c.x[0], a.x[1]+  b.x[1]+  c.x[1], a.x[2]+  b.x[2]+  c.x[2], a.x[3]+  b.x[3]+  c.x[3]);
  qd_real pn(a.x[0]-  b.x[0]+  c.x[0], a.x[1]-  b.x[1]+  c.x[1], a.x[2]-  b.x[2]+  c.x[2], a.x[3]-  b.x[3]+  c.x[3]);
  qd_real p2(a.x[0]+2*b.x[0]+4*c.x[0], a.x[1]+2*b.x[1]+4*c.x[1], a.x[2]+2*b.x[2]+4*c.x[2], a.x[3]+2*b.x[3]+4*c.x[3]);
  qd_real pi(c);
  qd_real q0(d);
  qd_real q1(d.x[0]+  e.x[0]+  f.x[0], d.x[1]+  e.x[1]+  f.x[1], d.x[2]+  e.x[2]+  f.x[2], d.x[3]+  e.x[3]+  f.x[3]);
  qd_real qn(d.x[0]-  e.x[0]+  f.x[0], d.x[1]-  e.x[1]+  f.x[1], d.x[2]-  e.x[2]+  f.x[2], d.x[3]-  e.x[3]+  f.x[3]);
  qd_real q2(d.x[0]+2*e.x[0]+4*f.x[0], d.x[1]+2*e.x[1]+4*f.x[1], d.x[2]+2*e.x[2]+4*f.x[2], d.x[3]+2*e.x[3]+4*f.x[3]);
  qd_real qi(f);
  mul447(p0, q0);
  mul447(p1, q1);
  mul447(pn, qn);
  mul447(p2, q2);
  mul447(pi, qi);
  //eval back
  // {{1,0,0,0,0},{1,1,1,1,1},{1,-1,1,-1,1},{1,2,4,8,16},{0,0,0,0,1}}^-1
  //   (  1    0    0    0   0)
  //   (-1/2   1  -1/3 -1/6  2)
  //   ( -1   1/2  1/2   0  -1)
  //   ( 1/2 -1/2 -1/6  1/6 -2)
  //   (  0    0    0    0   1)
  qd_real r0h(p0), r0l(q0);
  qd_real r2h(p0.x[0]*0.5-p1.x[0]*0.5+(p2.x[0]-pn.x[0])*(1/6.0)-2*pi.x[0],
              p0.x[1]*0.5-p1.x[1]*0.5+(p2.x[1]-pn.x[1])*(1/6.0)-2*pi.x[1],
              p0.x[2]*0.5-p1.x[2]*0.5+(p2.x[2]-pn.x[2])*(1/6.0)-2*pi.x[2],
              p0.x[3]*0.5-p1.x[3]*0.5+(p2.x[3]-pn.x[3])*(1/6.0)-2*pi.x[3]);
  qd_real r2l(q0.x[0]*0.5-q1.x[0]*0.5+(q2.x[0]-qn.x[0])*(1/6.0)-2*qi.x[0],
              q0.x[1]*0.5-q1.x[1]*0.5+(q2.x[1]-qn.x[1])*(1/6.0)-2*qi.x[1],
              q0.x[2]*0.5-q1.x[2]*0.5+(q2.x[2]-qn.x[2])*(1/6.0)-2*qi.x[2],
              q0.x[3]*0.5-q1.x[3]*0.5+(q2.x[3]-qn.x[3])*(1/6.0)-2*qi.x[3]);
  qd_real r1h(-r2h.x[0]+(p1.x[0]-pn.x[0])*0.5,
              -r2h.x[1]+(p1.x[1]-pn.x[1])*0.5,
              -r2h.x[2]+(p1.x[2]-pn.x[2])*0.5,
              -r2h.x[3]+(p1.x[3]-pn.x[3])*0.5);
  qd_real r1l(-r2l.x[0]+(q1.x[0]-qn.x[0])*0.5,
              -r2l.x[1]+(q1.x[1]-qn.x[1])*0.5,
              -r2l.x[2]+(q1.x[2]-qn.x[2])*0.5,
              -r2l.x[3]+(q1.x[3]-qn.x[3])*0.5);
  qd_real rnh(-p0.x[0]+(p1.x[0]+pn.x[0])*0.5-pi.x[0],
              -p0.x[1]+(p1.x[1]+pn.x[1])*0.5-pi.x[1],
              -p0.x[2]+(p1.x[2]+pn.x[2])*0.5-pi.x[2],
              -p0.x[3]+(p1.x[3]+pn.x[3])*0.5-pi.x[3]);
  qd_real rnl(-q0.x[0]+(q1.x[0]+qn.x[0])*0.5-qi.x[0],
              -q0.x[1]+(q1.x[1]+qn.x[1])*0.5-qi.x[1],
              -q0.x[2]+(q1.x[2]+qn.x[2])*0.5-qi.x[2],
              -q0.x[3]+(q1.x[3]+qn.x[3])*0.5-qi.x[3]); there is no q.x[3]
  qd_real rih(pi), ril(qi);
}
#endif

qd_real &qd_real::mul(qd_real const &other)
{
  //qd_real fromfft(*this);
  //fromfft.mul_fft(other);

  double p3, p4, p5;
  double q0;
  double x4;
  double x99=x[0]*other.x[3] + x[1]*other.x[2] + x[2]*other.x[1] + x[3]*other.x[0];

  p3=x[0];
  p4=x[1];
  p5=x[2];
  x99+= two_prod(p3, other.x[2])+
        two_prod(p4, other.x[1])+
        two_prod(p5, other.x[0]);

  x[2] = x[0]; //p2=x[2]
  x[3] = two_prod(x[2], other.x[1]); //q1=x[3]


  x4 = two_prod(x[1], other.x[0]); //p1=x[1]
  q0 = two_prod(x[0], other.x[0]); //p0=x[0]

  /* Start Accumulation */
  three_sum(x[1], x[2], q0);

  /* Six-Three Sum  of p2, q1, q2, p3, p4, p5. */
  three_sum(x[2], x[3], x4);
  three_sum(p3, p4, p5);
  /* compute (s0, s1, s2) = (p2, q1, q2) + (p3, p4, p5). */
  two_sum_(x[2], p3); //s0=p2=x[2]  t0=p3
  two_sum_(x[3], p4); //s1=q1 t1=p4
  two_sum_(x[3], p3);
  x4 += (p3 + p4 + p5);

  /* O(eps^3) order terms */
  x[3] += x99 + q0;
  renorm(x4);

  //if (x[0]!=fromfft.x[0] || x[1]!=fromfft.x[1] || x[2]!=fromfft.x[2] || x[3]!=fromfft.x[3])
  //  qd_dbgPoint();
  return *this;
}

qd_real &qd_real::sqr()
{
  double p4, p5;
  double q0, q1, q2, q3;
  //double s0, s1;
  //double t0, t1;

  p4 = 2.0 * x[0] * x[3];
  p5 = 2.0 * x[1] * x[2];

  x[3]=x[1];
  q3 = two_sqr(x[3]); //p=x
  q1 = two_prod(x[1], 2.0 * x[0]);
  q2 = two_prod(x[2], 2.0 * x[0]);
  q0 = two_sqr(x[0]);

  two_sum_(x[1], q0);

  two_sum_(q0, q1);
  two_sum_(x[2], x[3]);

  two_sum_(x[2], q0); //s0=p2  t0=q0
  two_sum_(x[3], q1); //s1=p3 t1=q1

  two_sum_(x[3], q0);
  q0 += q1;

  quick_nth_sum3x(q0);
  quick_nth_sum23(); //t1=p3
  quick_nth_sum3x(q0);

  two_sum_(p4, p5);
  two_sum_(q2, q3);

  two_sum_(p4, q2); //t0=p4 t1=q2
  q2 = q2 + p5 + q3;

  two_sum_(x[3], p4);
  p4 = p4 + q0 + q2;

  renorm(p4);
  return *this;
}

double qd_real::radixfloor() const
{
  int ilog1=std::ilogb(x[0]);
  return ldexp(1, ilog1);
}

qd_real &qd_real::recip()
{
  qd_real a(1, 0, 0, 0);
  qd_real b(*this);

  qd_real t;

  x[0] = a.x[0] / b.x[0];
  t=b;
  t.mul_double(-x[0]);
  a.add(t);

  x[1] = a.x[0] / b.x[0];
  t=b;
  t.mul_double(-x[1]);
  a.add(t);

  x[2] = a.x[0] / b.x[0];
  t=b;
  t.mul_double(-x[2]);
  a.add(t);

  x[3] = a.x[0] / b.x[0];

  renorm();
  return *this;
}

qd_real &qd_real::sqrt()
{
  if (x[0]<=0)
  {
    x[3]=x[2]=x[1]=x[0]=0;
    return *this;
  };

  /* Strategy:

     Perform the following Newton iteration:

       x' = x + (1 - a * x^2) * x / 2;

     which converges to 1/sqrt(a), starting with the
     double precision approximation to 1/sqrt(a).
     Since Newton's iteration more or less doubles the
     number of correct digits, we only need to perform it
     twice.
  */

  qd_real h = *this;
  h.lshift(-1);
  h.chs();
  qd_real r(1.0 / std::sqrt(x[0]), 0, 0, 0);
  qd_real t;

  for (int cykl=0; cykl<3; cykl++)
  {
    //r += ((0.5 - h * sqr(r)) * r);  with h->-h
    t=r;
    t.sqr();
    t.mul(h);
    t.add_double(0.5);
    t.mul(r);
    r.add(t);
  }

  mul(r);
  return *this;
}

/* Computes the nearest integer to d. */
static inline double nint(double d) {
  if (d == std::floor(d))
    return d;
  return std::floor(d + 0.5);
}

qd_real &qd_real::round()
{
#if 1 //Bailey's method
  qd_real a=*this;

  x[0] = nint(a.x[0]);
  x[3]=x[2]=x[1]=0.0;

  if (x[0] != a.x[0])
  {
    /* First double is not an integer. */
    if (std::abs(x[0] - a.x[0]) == 0.5 && a.x[1] < 0.0) {
        x[0] -= 1.0;
    }
  }
  else
  {
    /* First double is already an integer. */
    x[1] = nint(a.x[1]);

    if (x[1] != a.x[1]) {
        if (std::abs(x[1] - a.x[1]) == 0.5 && a.x[2] < 0.0) {
            x[1] -= 1.0;
        }

    } else {
        /* Second double is already an integer. */
        x[2] = nint(a.x[2]);

        if (x[2] != a.x[2]) {
            if (std::abs(x[2] - a.x[2]) == 0.5 && a.x[3] < 0.0) {
                x[2] -= 1.0;
            }
        } else {
            /* Third double is already an integer. */
            x[3] = nint(a.x[3]);
        }
    }
  }

  renorm();
#else
  my version of dd_real
  //if lo==+-0.5 exactly, we should decide based on hi whether to round to +-1 or 0 but who cares
  //or could do this+=2^103-=2^103 or something like that
  if ((lo_<=-0.5) || (lo_>=0.5))
  {
    //hi=round(hi);
    lo_=std::round(lo_);
    int discard;
    double testme=std::frexp(lo_, &discard);
    if (testme==-1 || testme==-0.5 || testme==0.5 || testme==1) //spec says it returns +-0.5 instead of 1 but that's just stupid
      dd_dbgPoint(); //may need renormalization
  }
  else
  {
    hi=std::round(hi);
    lo_=0;
    //should be good even at hi= +-2^52
  }
#endif
  return *this;
}

qd_real &qd_real::frac()
{
  bool neg=x[0]<0;
  if (x[3]<=-0.5 || x[3]>=0.5)
  {
    x[0]=x[3];
    x[3]=x[2]=x[1]=0;
  }
  else if (x[2]<=-0.5 || x[2]>=0.5)
  {
    x[0]=x[2];
    x[1]=x[3];
    x[3]=x[2]=0;
  }
  else if (x[1]<=-0.5 || x[1]>=0.5)
  {
    x[0]=x[1];
    x[1]=x[2];
    x[2]=x[3];
    x[3]=0;
  }
  else
  {
  }
  if (neg)
    add_double(-std::ceil(x[0]));
  else
    add_double(-std::floor(x[0]));
  return *this;
}

qd_real &qd_real::mod1()
{
  if (x[3]<=-0.5 || x[3]>=0.5)
  {
    x[0]=x[3];
    x[3]=x[2]=x[1]=0;
  }
  else if (x[2]<=-0.5 || x[2]>=0.5)
  {
    x[0]=x[2];
    x[1]=x[3];
    x[3]=x[2]=0;
  }
  else if (x[1]<=-0.5 || x[1]>=0.5)
  {
    x[0]=x[1];
    x[1]=x[2];
    x[2]=x[3];
    x[3]=0;
  }
  else
  {
  }
  add_double(-std::floor(x[0]));
  return *this;
}

bool qd_real::reduce_angle()
{
  //static const dd_real M_PIqd(7074237752028440.0/2251799813685248.0, 2483878800010755.0/20282409603651670423947251286016.0);
  static const qd_real M_PIqd_chs(-M_PIqd.x[0], -M_PIqd.x[1], -M_PIqd.x[2], -M_PIqd.x[3]);
  static const qd_real M_PIqd_2(2*M_PIqd.x[0], 2*M_PIqd.x[1], 2*M_PIqd.x[2], 2*M_PIqd.x[3]);
  //static const dd_real M_PIdd_2chs(-2*3.141, -2*0.0005926);
  if (!M_PIqd_chs.isle(*this)) // this < -pi
  {
    add(M_PIqd_2);
    return true;
  }
  else if (M_PIqd.isle(*this))
  {
    sub(M_PIqd_2);
    return true;
  }
  else
    return false;
}

qd_real &qd_real::add_pi(double x)
{
  qd_real pix(M_PIqd);
  pix.mul_double(x);
  add(pix);
  return *this;
}

std::strong_ordering qd_real::compare(qd_real const &other) const
{
  for (int i=0; i<4; i++)
  {
    if (x[i]<other.x[i])
      return std::strong_ordering::less;//-1;
    else if (x[i]>other.x[i])
      return std::strong_ordering::greater;
  }
  return std::strong_ordering::equal;
}

bool qd_real::isequal(qd_real const &other) const
{
  return (x[0]==other.x[0]) && (x[1]==other.x[1]) && (x[2]==other.x[2]) && (x[3]==other.x[3]);
}

bool qd_real::is0() const
{
  return x[0]==0;
}

bool qd_real::isle(qd_real const &other) const
{
  for (int i=0; i<4; i++)
    if (x[i]!=other.x[i])
      return x[i]<other.x[i];
  return true;
}

bool qd_real::isle0() const
{
  return x[0]<=0;
}

bool qd_real::isl0() const
{
  return x[0]<0;
}

bool qd_real::isl1() const
{
  return x[0]<1;
}

qd_real &qd_real::min(qd_real const &other)
{
  if (!isle(other))
  {
    assign(other);
  };
  return *this;
}

namespace {
// TODO: Didn't review thoroughly; we'll need actual BigNum for arbitrary precision anyway
// Unsigned big integer (little-endian 32-bit limbs), just enough for exact decimal I/O.
struct BigU
{
  std::vector<uint32_t> w;
  BigU() { }
  explicit BigU(uint64_t v) { while (v) { w.push_back(uint32_t(v)); v>>=32; } }
  bool is0() const { return w.empty(); }
  void trim() { while (!w.empty() && w.back()==0) w.pop_back(); }
  int bitlen() const { return w.empty()?0:32*(int(w.size())-1)+32-std::countl_zero(w.back()); }
  bool bit(int i) const { size_t k=size_t(i)/32; return i>=0 && k<w.size() && ((w[k]>>(i%32))&1); }
  bool anyBelow(int n) const // any of the lowest n bits set
  {
    for (int k=0; k<n/32 && k<int(w.size()); k++)
      if (w[k]) return true;
    return n%32 && size_t(n/32)<w.size() && (w[n/32] & ((1u<<(n%32))-1));
  }
  int ctz() const { int k=0; while (w[k]==0) k++; return 32*k+std::countr_zero(w[k]); }
  uint64_t low64() const { return (w.size()>0?w[0]:0) | (w.size()>1?uint64_t(w[1])<<32:0); }
  void mul_add(uint32_t m, uint32_t a)
  {
    uint64_t c=a;
    for (uint32_t &x: w) { c+=uint64_t(x)*m; x=uint32_t(c); c>>=32; }
    if (c) w.push_back(uint32_t(c));
  }
  uint32_t divmod(uint32_t d)
  {
    uint64_t r=0;
    for (size_t i=w.size(); i-->0;) { r=(r<<32)|w[i]; w[i]=uint32_t(r/d); r%=d; }
    trim();
    return uint32_t(r);
  }
  void shl(int n)
  {
    if (w.empty() || n<=0) return;
    int bits=n%32;
    if (bits)
    {
      uint32_t c=0;
      for (uint32_t &x: w) { uint32_t nc=x>>(32-bits); x=(x<<bits)|c; c=nc; }
      if (c) w.push_back(c);
    }
    w.insert(w.begin(), size_t(n/32), 0u);
  }
  void shr(int n)
  {
    if (n<=0) return;
    if (size_t(n/32)>=w.size()) { w.clear(); return; }
    w.erase(w.begin(), w.begin()+n/32);
    int bits=n%32;
    if (bits)
    {
      for (size_t i=0; i<w.size(); i++)
        w[i]=(w[i]>>bits) | (i+1<w.size()?w[i+1]<<(32-bits):0);
    }
    trim();
  }
  static int cmp(BigU const &a, BigU const &b)
  {
    if (a.w.size()!=b.w.size()) return a.w.size()<b.w.size()?-1:1;
    for (size_t i=a.w.size(); i-->0;)
      if (a.w[i]!=b.w[i]) return a.w[i]<b.w[i]?-1:1;
    return 0;
  }
  void add(BigU const &b)
  {
    if (w.size()<b.w.size()) w.resize(b.w.size(), 0);
    uint64_t c=0;
    for (size_t i=0; i<w.size(); i++) { c+=uint64_t(w[i])+(i<b.w.size()?b.w[i]:0); w[i]=uint32_t(c); c>>=32; }
    if (c) w.push_back(uint32_t(c));
  }
  void sub(BigU const &b) // requires *this>=b
  {
    int64_t c=0;
    for (size_t i=0; i<w.size(); i++) { c+=int64_t(w[i])-(i<b.w.size()?b.w[i]:0); w[i]=uint32_t(c); c>>=32; }
    trim();
  }
  std::string toDecimal() const
  {
    if (w.empty()) return "0";
    BigU t(*this);
    std::string s;
    while (!t.is0())
    {
      uint32_t r=t.divmod(1000000000);
      for (int i=0; i<9; i++) { s.push_back(char('0'+r%10)); r/=10; }
    }
    while (s.size()>1 && s.back()=='0') s.pop_back();
    return std::string(s.rbegin(), s.rend());
  }
};

// |q| = mag*2^e2 exactly, mag odd (or 0); returns the sign
bool qd_exact(qd_real const &q, BigU &mag, int &e2)
{
  int64_t m[4];
  int e[4];
  int n=0;
  for (double xi: q.x)
    if (xi!=0)
    {
      int ex;
      double f=std::frexp(xi, &ex);
      m[n]=(int64_t)std::ldexp(f, 53);
      e[n]=ex-53;
      n++;
    }
  mag.w.clear();
  e2=0;
  if (n==0)
    return std::signbit(q.x[0]);
  e2=*std::min_element(e, e+n);
  BigU pos, neg;
  for (int i=0; i<n; i++)
  {
    BigU t(uint64_t(m[i]<0?-m[i]:m[i]));
    t.shl(e[i]-e2);
    (m[i]<0?neg:pos).add(t);
  }
  bool negative=BigU::cmp(pos, neg)<0;
  if (negative) { neg.sub(pos); mag=neg; } else { pos.sub(neg); mag=pos; }
  if (!mag.is0())
  {
    int tz=mag.ctz();
    mag.shr(tz);
    e2+=tz;
  }
  return negative;
}

// Greedy round-to-nearest split of (Q+sticky)*2^e into up to 4 doubles
qd_real qd_from_big(BigU Q, int e, bool sticky, bool negative)
{
  qd_real r;
  for (int i=0; i<4 && !Q.is0(); i++)
  {
    int sh=std::max(0, Q.bitlen()-53);
    BigU top(Q);
    top.shr(sh);
    bool up=sh>0 && Q.bit(sh-1) && (Q.anyBelow(sh-1) || sticky || (top.low64()&1));
    uint64_t t=top.low64()+(up?1:0);
    r.x[i]=std::ldexp(negative?-double(t):double(t), e+sh);
    BigU tb(t);
    tb.shl(sh);
    if (up) { tb.sub(Q); Q=tb; negative=!negative; sticky=false; }
    else Q.sub(tb);
  }
  return r;
}

// value = digits*10^dexp; %g-like formatting with all of the given digits
std::string formatDecimal(bool negative, std::string const &digits, int dexp)
{
  int prec=int(digits.size());
  int x=prec-1+dexp; // decimal exponent of the first digit
  std::string s=negative?"-":"";
  if (x<-4 || x>=std::max(prec, 1))
  {
    s+=digits[0];
    if (prec>1) { s+='.'; s+=digits.substr(1); }
    char buf[16];
    std::snprintf(buf, sizeof(buf), "e%c%02d", x<0?'-':'+', x<0?-x:x);
    s+=buf;
  }
  else if (x<0)
    s+="0."+std::string(-x-1, '0')+digits;
  else if (x+1>=prec)
    s+=digits+std::string(x+1-prec, '0');
  else
    s+=digits.substr(0, x+1)+"."+digits.substr(x+1);
  return s;
}

} // namespace

bool qd_real::fromString(char const *str)
{
  char const *p=str;
  while (std::isspace((unsigned char)*p)) p++;
  bool negative=false;
  if (*p=='+' || *p=='-') negative=*p++=='-';
  BigU D;
  int ndig=0, dexp=0;
  bool seen=false;
  for (bool frac=false; ; p++)
  {
    if (*p>='0' && *p<='9')
    {
      seen=true;
      if (D.is0() && *p=='0')
      {
        if (frac) dexp--;
        continue;
      }
      if (ndig<120) { D.mul_add(10, uint32_t(*p-'0')); ndig++; if (frac) dexp--; }
      else if (!frac) dexp++; // digits beyond 120 are below qd precision
    }
    else if (*p=='.' && !frac)
      frac=true;
    else
      break;
  }
  bool ok=seen;
  if (ok && (*p=='e' || *p=='E'))
  {
    char *end;
    long ex=std::strtol(p+1, &end, 10);
    if (end==p+1 || std::isspace((unsigned char)p[1]))
      ok=false;
    else
    {
      p=end;
      dexp+=int(std::clamp(ex, -100000L, 100000L));
    }
  }
  while (ok && std::isspace((unsigned char)*p)) p++;
  if (!ok || *p)
  { // inf, nan, hex floats
    __float128 q;
    if (!float128_fromString(str, q))
      return false;
    double h=(double)q;
    if (!std::isfinite(h))
    {
      zero(h);
      return true;
    }
    __float128 rest=q-h;
    double m=(double)rest;
    *this=qd_real(h, m, (double)(rest-m), 0);
    return true;
  }
  if (D.is0() || ndig+dexp<-400)
  {
    zero(negative?-0.0:0.0);
    return true;
  }
  if (ndig+dexp>400)
  {
    zero(negative?-INFINITY:INFINITY);
    return true;
  }
  if (dexp>=0)
  {
    for (int i=0; i<dexp; i++) D.mul_add(10, 0);
    *this=qd_from_big(D, 0, false, negative);
    return true;
  }
  BigU B(1);
  for (int i=0; i<-dexp; i++) B.mul_add(10, 0); // TODO: 400 times?
  // Q=floor(D*2^s/B) with ~256 bits
  int s=256+B.bitlen()-D.bitlen();
  if (s>=0) D.shl(s); else B.shl(-s);
  BigU Q;
  int shift=D.bitlen()-B.bitlen();
  if (shift>=0)
  {
    B.shl(shift);
    Q.w.assign(size_t(shift/32+1), 0u);
    for (int i=shift; i>=0; i--)
    {
      if (BigU::cmp(D, B)>=0)
      {
        D.sub(B);
        Q.w[size_t(i/32)]|=1u<<(i%32);
      }
      B.shr(1);
    }
    Q.trim();
  }
  *this=qd_from_big(Q, -s, !D.is0(), negative);
  return true;
}

QString qd_real::toString() const
{
  if (!std::isfinite(x[0]))
    return std::isnan(x[0])?"nan":x[0]<0?"-inf":"inf";
  BigU mag;
  int e2;
  bool negative=qd_exact(*this, mag, e2);
  if (mag.is0())
    return negative?"-0":"0";
  BigU const exact=mag;
  int dexp=0;
  if (e2>=0)
    mag.shl(e2);
  else
  {
    for (int i=0; i<-e2; i++) mag.mul_add(5, 0);
    dexp=e2;
  }
  std::string all=mag.toDecimal();
  while (all.size()>1 && all.back()=='0') { all.pop_back(); dexp++; }
  for (size_t prec=1; prec<all.size(); prec++)
  {
    std::string d=all.substr(0, prec);
    char next=all[prec];
    bool rest=all.find_first_not_of('0', prec+1)!=std::string::npos;
    int de=dexp+int(all.size()-prec);
    if (next>'5' || (next=='5' && (rest || (d.back()-'0')%2)))
    {
      int i=int(prec)-1;
      while (i>=0 && d[i]=='9') d[i--]='0';
      if (i<0) { d.insert(d.begin(), '1'); d.pop_back(); de++; }
      else d[i]++;
    }
    while (d.size()>1 && d.back()=='0') { d.pop_back(); de++; }
    std::string s=formatDecimal(negative, d, de);
    qd_real back;
    BigU bmag;
    int be2;
    if (back.fromString(s.c_str()) && qd_exact(back, bmag, be2)==negative && be2==e2 && BigU::cmp(bmag, exact)==0)
      return QString::fromStdString(s);
  }
  return QString::fromStdString(formatDecimal(negative, all, dexp));
}

int qd_real::toRound() const
{
  return std::floor(x[0]+0.5)+std::floor(x[1]+0.5)+std::floor(x[2]+0.5)+std::floor(x[3]+0.5);
}

double qd_real::toDouble() const
{
  return x[0];
}

#undef MAKE_DOUBLE

} // namespace MandelMath
