#ifndef SHAREABLEIMAGEWRAPPER_H
#define SHAREABLEIMAGEWRAPPER_H

#include <QImage>

class ShareableImageWrapper
{
public:
  ShareableImageWrapper(): image(nullptr) { }
  ShareableImageWrapper(QImage *image);
  QImage *image;
};

Q_DECLARE_METATYPE(ShareableImageWrapper);

#endif // SHAREABLEIMAGEWRAPPER_H
