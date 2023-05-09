QT += quick

CONFIG += c++20
CONFIG += force_debug_info

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

#helps with templated version but not much with virtual version
#removes dbgPoint :-(
#QMAKE_CXXFLAGS_DEBUG += -O3
#QMAKE_CFLAGS_DEBUG += -O3
QMAKE_CXXFLAGS_DEBUG += -fconcepts-diagnostics-depth=9
#QMAKE_CXXFLAGS_RELEASE += -O3
#QMAKE_CXXFLAGS_RELEASE_WITH_DEBUGINFO += -O3
#QMAKE_CFLAGS_RELEASE += -O3
#QMAKE_CFLAGS_RELEASE_WITH_DEBUGINFO += -O3

SOURCES += \
        JuliaModel.cpp \
        LaguerreModel.cpp \
        MandelEvaluator.cpp \
        MandelImageCombiner.cpp \
        MandelMath.cpp \
        MandelModel.cpp \
        ShareableImageWrapper.cpp \
        double_double.cpp \
        main.cpp \
        multiprec.cpp

RESOURCES += qml.qrc

# Additional import path used to resolve QML modules in Qt Creator's code model
QML_IMPORT_PATH =

# Additional import path used to resolve QML modules just for Qt Quick Designer
QML_DESIGNER_IMPORT_PATH =

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
  JuliaModel.hpp \
  LaguerreModel.hpp \
  MandelEvaluator.hpp \
  MandelImageCombiner.hpp \
  MandelMath.hpp \
  MandelModel.hpp \
  ShareableImageWrapper.hpp \
  double_double.hpp \
  multiprec.hpp
