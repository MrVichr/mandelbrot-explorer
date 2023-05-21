#ifndef MANDELMATH_DOUBLE_DOUBLE_HPP
#define MANDELMATH_DOUBLE_DOUBLE_HPP

#include <compare>
#include <QString>

namespace MandelMath {

void fpu_fix_start(unsigned int *old_cw);
void fpu_fix_end(unsigned int *old_cw);

class dd_real {
protected:
  inline void quick_two_sum(double a, double b);
  inline void two_sum(double a, double b);
  //made public for debugging inline void split(double a);
  inline void two_prod(double a, double b);
  inline void two_sqr(double a);
public:
  /*inline*/ void split(double a);
  dd_real(): hi(0), lo_(0) { }
  dd_real(double h, double l): hi(h), lo_(l) { }
  double hi;
  double lo_; // !!can have other sign than hi!!
  double eps2() const;
  double eps234() const;
  dd_real &assign(const dd_real &src) noexcept { hi=src.hi; lo_=src.lo_; return *this; }
  dd_real &zero(double v) { hi=v; lo_=0; return *this; }
  dd_real &chs() { hi=-hi; lo_=-lo_; return *this; }
  dd_real &lshift(int exp); //*=2^exp
  dd_real &add_double(double h2);
  dd_real &mul_double(double h2);
  //dd_real &add(double h2, double l2); //TODO: -> add(dd_real)
  dd_real &add(dd_real const &other);
  dd_real &sub(dd_real const &other);
  dd_real &rsub(dd_real const &other) { chs(); return add(other); }
  //dd_real &mul(double h2, double l2);
  dd_real &mul(dd_real const &other);
  dd_real &sqr();
  double radixfloor() const; //nearest smaller power of 2 (1.5->1->1)
  dd_real &recip();
  dd_real &sqrt();
  dd_real &round();
  dd_real &frac();
  dd_real &mod1();
  bool reduce_angle();
  dd_real &add_pi(double x);
  //int compare(const dd_real *other) const; //return -1 if <, 0 if =, +1 if >
  std::strong_ordering compare(dd_real const &other) const; //return -1 if <, 0 if =, +1 if >; std::strong_ordering not in my compiler yet
  //int compare(double other_h, double other_l) const;
  bool isequal(dd_real const &other) const;
  bool is0() const;
  bool isle(const dd_real &other) const;
  bool isle0() const;
  bool isl0() const;
  bool isl1() const;
  dd_real &min(dd_real const &other);

  //explicit operator double();
  QString toString() const;
  int toRound() const;
  double toDouble() const;
};

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_DOUBLE_HPP
