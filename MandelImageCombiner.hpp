#ifndef MANDELIMAGECOMBINER_H
#define MANDELIMAGECOMBINER_H

//#include <QObject>
//#include <QQuickItem>
#include <QImage>
#include <QQuickPaintedItem>

#include "ShareableImageWrapper.hpp"

class MandelImageCombiner: public QQuickPaintedItem
{
  Q_OBJECT
public:
  MandelImageCombiner(QQuickItem *parent = nullptr);
  Q_INVOKABLE void resetBgImage(int width, int height, unsigned int color);
  //Q_INVOKABLE void setBgImage(const QImage &image);
  //Q_INVOKABLE void setBaseImage(const QImage &image);
  Q_INVOKABLE ShareableImageWrapper getBaseImage();
  Q_INVOKABLE ShareableImageWrapper getOverlayImage();
  //Q_INVOKABLE void setFgImage(QImage &image);
  void paint(QPainter *painter) override;

private:
    QImage m_bgImage;
    QImage m_baseImage;
    QImage m_fgImage;
};

#endif // MANDELIMAGECOMBINER_H
