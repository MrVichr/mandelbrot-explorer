
#include "atomic_update_rect.hpp"

void atomic_min(std::atomic<int> &a, int val) noexcept
{
    int prev_value = a;
    while(prev_value > val &&
           !a.compare_exchange_weak(prev_value, val))
    {}
}

void atomic_max(std::atomic<int> &a, int val) noexcept
{
    int prev_value = a;
    while(prev_value < val &&
           !a.compare_exchange_weak(prev_value, val))
    {}
}

AtomicUpdateRect::AtomicUpdateRect()
{

}

void AtomicUpdateRect::addPoint(int x, int y)
{
  atomic_min(left, x);
  atomic_max(right, x);
  atomic_min(top, y);
  atomic_max(bottom, y);
}

void AtomicUpdateRect::fetch()
{
  //need some std::atomic expert here...
  fleft=bright;
  fright=bleft-1;
  ftop=bbottom;
  fbottom=btop-1;
  do
  {
    int left=this->left.exchange(bright);
    int right=this->right.exchange(bleft-1);
    int top=this->top.exchange(bbottom);
    int bottom=this->bottom.exchange(btop-1);
    if (left==bright && right==bleft-1 && top==bbottom && bottom==btop-1)
      break;
    if (fleft>left)
      fleft=left;
    if (fright<right)
      fright=right;
    if (ftop>top)
      ftop=top;
    if (fbottom<bottom)
      fbottom=bottom;
  }
  while (true);
}

