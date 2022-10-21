#include "MandelImageCombiner.hpp"
#include <QPainter>

MandelImageCombiner::MandelImageCombiner(QQuickItem *parent):
  QQuickPaintedItem(parent)
{

}

void MandelImageCombiner::resetBgImage(int width, int height, unsigned int color)
{
  if ((width<=0) || (height<=0)) //Qt begins with width=0, height=-13
    return;
  {
    QImage newBg(width, height, QImage::Format::Format_ARGB32);
    QPainter bgPainter(&newBg);
    bgPainter.fillRect(0+0, 0+0, width-0, height-0, QColor(color));
    m_bgImage.swap(newBg);
  }

  {
    QImage newBase(width, height, QImage::Format::Format_ARGB32);
    QPainter basePainter(&newBase);
    basePainter.fillRect(0, 0, width, height, Qt::GlobalColor::transparent);
    m_baseImage.swap(newBase);
  }
}

/*void MandelImageCombiner::setBgImage(const QImage &image)
{
  m_bgImage = image;
  update();
}

void MandelImageCombiner::setBaseImage(const QImage &image)
{
  m_baseImage = image;
  update();
}*/
ShareableImageWrapper MandelImageCombiner::getBaseImage()
{
  return ShareableImageWrapper(&m_baseImage);
}

ShareableImageWrapper MandelImageCombiner::getOverlayImage()
{
  return ShareableImageWrapper(&m_fgImage);
}

/*void MandelImageCombiner::setFgImage(QImage &image)
{
  m_fgImage.swap(image);
  update();
}*/

void MandelImageCombiner::paint(QPainter *painter)
{
  painter->setPen(QColor("red"));
  painter->drawEllipse(50, 50, 100, 50);
  painter->drawImage(0, 0, m_bgImage);
  //painter->setPen(QColor(0x0000ff));
  //painter->drawEllipse(50, 100, 100, 50);
  painter->drawImage(0, 0, m_baseImage);
  painter->drawImage(0, 0, m_fgImage);
  //painter->setPen(QColor("green"));
  //painter->drawEllipse(50, 150, 100, 50);
}
