#ifndef MANDELMATH_MULTIPREC_HPP
#define MANDELMATH_MULTIPREC_HPP


namespace MandelMath {

class multiprec
{
public:
  multiprec();
  void set(double val);
  void assign(const multiprec &src) { (void)src; }
  void chs();
  void lshift(int shoft);
  void round();
  void frac();
  void mod1();
  void add_double(double x);
  void mul_double(double x);
  void add(multiprec *other);
  void mul(multiprec *other);
  void sqr();
  double radixfloor(); //nearest smaller power of 2 (1.5->1->1)
  void recip();
  void sqrt();
  int toround();
  int compare(const multiprec *other); //return -1 if <, 0 if =, +1 if >
  bool isequal(const multiprec *other);
  bool is0();
  bool isle(const multiprec *other);
  bool isle0();
  bool isl0();
  bool isl1();
  double toDouble();
};


} // namespace MandelMath

#endif // MANDELMATH_MULTIPREC_HPP
