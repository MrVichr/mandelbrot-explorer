#ifndef MANDELMATH_DOUBLE_DOUBLE_HPP
#define MANDELMATH_DOUBLE_DOUBLE_HPP

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
  void assign(const dd_real &src) { hi=src.hi; lo_=src.lo_; }
  void zero(double v) { hi=v; lo_=0; }
  void chs() { hi=-hi; lo_=-lo_; }
  void lshift(int exp); //*=2^exp
  void add_double(double h2);
  void mul_double(double h2);
  void add(double h2, double l2); //TODO: -> add(dd_real)
  void mul(double h2, double l2);
  void sqr();
  double radixfloor() const; //nearest smaller power of 2 (1.5->1->1)
  void recip();
  void sqrt();
  void round();
  void frac();
  void mod1();
  int compare(const dd_real *other) const; //return -1 if <, 0 if =, +1 if >
  int compare(double other_h, double other_l) const;
  bool isequal(const dd_real *other) const;
  bool is0() const;
  bool isle(const dd_real *other) const;
  bool isle0() const;
  bool isl0() const;
  bool isl1() const;

  //explicit operator double();
};

} // namespace MandelMath

#endif // MANDELMATH_NUMBER_DOUBLE_HPP
