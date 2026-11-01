#ifndef MANDELMATH_QUAD_DOUBLE_HPP
#define MANDELMATH_QUAD_DOUBLE_HPP

#include <compare>
#include <QString>

namespace MandelMath {

void fpu_fix_start(unsigned int *old_cw);
void fpu_fix_end(unsigned int *old_cw);

class qd_real {
protected:
  //inline void quick_two_sum(double a, double b);
  inline void quick_nth_sum01();//double a, double b)
  inline void quick_nth_sum12();
  inline void quick_nth_sum23();
  inline void quick_nth_sum3x(double &x4);
  inline void renorm();
  inline void renorm(double &x4);
  inline static void two_sum_(double &a, double &b);
  inline void two_sum0(double &b);
  inline void two_sum1(double &b);
  inline void two_sum2(double &b);
  inline void two_sum3(double &b);
  inline static void three_sum(double &a, double &b, double &c);
  inline static void three_sum2(double &a, double &b, double &c);
  //made public for debugging inline void split(double a);
  inline static double two_prod(double &a, double b); //a*b->[a,result]
  inline static double two_sqr(double &a);
public:
  /*inline*/ static void split(double a, double &hi, double &lo);
  qd_real(): x{0,0,0,0} { }
  qd_real(double h, double x1, double x2, double l): x{h, x1, x2, l} { }
  double x[4]; // !!can have different signs
  double eps2() const;
  double eps234() const;
  qd_real &assign(const qd_real &src) noexcept { *this=src; return *this; }
  qd_real &zero(double v=0) { x[0]=v; x[1]=x[2]=x[3]=0; return *this; }
  qd_real &chs() { for (double &xi: x) {xi=-xi;}; return *this; }
  qd_real &lshift(int exp); //*=2^exp
  qd_real &add_double(double b);
  qd_real &mul_double(double h2);
  //qd_real &add(double h2, double l2); //TODO: -> add(dd_real)
  qd_real &add(qd_real const &other);
  qd_real &sub(qd_real const &other);
  qd_real &rsub(qd_real const &other) { chs(); return add(other); }
  //qd_real &mul(double h2, double l2);
  qd_real &mul(qd_real const &other);
  qd_real &sqr();
  double radixfloor() const; //nearest smaller power of 2 (1.5->1->1)
  qd_real &recip();
  qd_real &sqrt();
  qd_real &round();
  qd_real &frac();
  qd_real &mod1();
  bool reduce_angle();
  qd_real &add_pi(double x);
  //int compare(const dd_real *other) const; //return -1 if <, 0 if =, +1 if >
  std::strong_ordering compare(qd_real const &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  //int compare(double other_h, double other_l) const;
  bool isequal(qd_real const &other) const;
  bool is0() const;
  bool isle(const qd_real &other) const;
  bool isle0() const;
  bool isl0() const;
  bool isl1() const;
  qd_real &min(qd_real const &other);

  //explicit operator double();
  QString toString() const;
  int toRound() const;
  double toDouble() const;
};

} // namespace MandelMath

#endif // MANDELMATH_QUAD_DOUBLE_HPP
