#ifndef DD_HEAP_HPP
#define DD_HEAP_HPP

#include <cassert>
#include <utility>

#include <QDebug>

#include "double_double.hpp"

namespace MandelMath {

// dd_real kept on the heap, the pattern for heap-backed types such as MPFR.
// Buffers are recycled through a per-thread free list, so temporaries rarely reach new/delete.
// Point arrays and convert_block use plain dd_real as the storage form.
class dd_heap
{
  union Block
  {
    dd_real value;
    Block *next;
    Block(): value() { }
  };
  struct Pool
  {
    Block *head=nullptr;
    ~Pool()
    {
      while (head)
      {
        Block *n=head->next;
        delete head;
        head=n;
      }
      dead()=true;
    }
  };
  static bool &dead()
  {
    thread_local bool d=false; // trivially destructible: still valid while the thread's Pool is torn down
    return d;
  }
  static Pool &pool()
  {
    thread_local Pool p;
    return p;
  }
  static Block *acquire()
  {
    if (dead())
      return new Block();
    Pool &p=pool();
    Block *b=p.head;
    if (!b)
      return new Block();
    p.head=b->next;
    return b;
  }
  static void release(Block *b)
  {
    if (dead())
    {
      delete b;
      return;
    }
    Pool &p=pool();
    b->next=p.head;
    p.head=b;
  }

  Block *p; // nullptr only in a moved-from object, which should only be destroyed
  void ensure_block() noexcept
  {
    if (!p) [[unlikely]]
    {
      assert(!"assignment to a moved-from dd_heap");
      qDebug() << "assignment to a moved-from dd_heap";
      p=acquire();
    }
  }
public:
  dd_heap(): p(acquire()) { p->value=dd_real(); }
  dd_heap(dd_heap const &o) noexcept: p(acquire()) { p->value=o.p->value; }
  dd_heap(dd_heap &&o) noexcept: p(o.p) { o.p=nullptr; }
  // assignment copies the value, so a long-lived object keeps its own buffer
  dd_heap &operator=(dd_heap const &o) noexcept
  {
    ensure_block();
    p->value=o.p->value;
    return *this;
  }
  dd_heap &operator=(dd_heap &&o) noexcept { return *this=o; }
  dd_heap &operator=(dd_real const &v) noexcept
  {
    ensure_block();
    p->value=v;
    return *this;
  }
  explicit operator dd_real() const { return p->value; }
  ~dd_heap() { if (p) release(p); }
  friend void swap(dd_heap &a, dd_heap &b) noexcept { std::swap(a.p, b.p); }
  dd_real &value() { return p->value; }
  dd_real const &value() const { return p->value; }
};

} // namespace MandelMath

#endif // DD_HEAP_HPP
