#ifndef ATOMICUPDATERECT_HPP
#define ATOMICUPDATERECT_HPP

#include <atomic>



class AtomicUpdateRect
{
protected:
    int bleft=0, bright=0, btop=0, bbottom=0;
    std::atomic<int> left, right, top, bottom; //left<=right top<=bottom
public:
    AtomicUpdateRect();
    void setBase(int l, int r, int t, int b)
    {
        bleft=l; bright=r; btop=t; bbottom=b;
    }
    void invalidate()
    {
        left=bleft;
        right=bright-1;
        top=btop;
        bottom=bbottom-1;
    }
    void addPoint(int x, int y);
    int fleft=0, fright=0, ftop=0, fbottom=0;
    void fetch(); //ideally would return QRect<int>
};

#endif // ATOMICUPDATERECT_HPP
