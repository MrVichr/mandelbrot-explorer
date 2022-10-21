#include "multiprec.hpp"
#include "MandelMath.hpp" //dbgPoint

namespace MandelMath {

multiprec::multiprec()
{
  dbgPoint();
}

void multiprec::set(double val)
{
  (void)val;
}

void multiprec::chs()
{

}

void multiprec::lshift(int shoft)
{
  (void)shoft;
}

void multiprec::round()
{

}

void multiprec::frac()
{

}

void multiprec::mod1()
{

}

void multiprec::add_double(double x)
{
  (void)x;
}

void multiprec::mul_double(double x)
{
  (void)x;
}

void multiprec::add(multiprec *other)
{
  (void)other;
}

void multiprec::mul(multiprec *other)
{
  (void)other;
}

void multiprec::sqr()
{
  mul(this);
}

double multiprec::radixfloor()
{
  return 1;
}

void multiprec::recip()
{
}

void multiprec::sqrt()
{
}

int multiprec::toround()
{
  return 0;
}

int multiprec::compare(const multiprec *other)
{
  (void)other;
  return 0;
}

bool multiprec::isequal(const multiprec *other)
{
  (void)other;
  return true;
}

bool multiprec::is0()
{
  return false;
}

bool multiprec::isle(const multiprec *other)
{
  (void)other;
  return false;
}

bool multiprec::isle0()
{
  return false;
}

bool multiprec::isl0()
{
  return false;
}

bool multiprec::isl1()
{
  return false;
}

double multiprec::toDouble()
{
  return 0;
}

} // namespace MandelMath
