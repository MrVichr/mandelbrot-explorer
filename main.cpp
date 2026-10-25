#include <QGuiApplication>
#include <QQmlApplicationEngine>

#include "MandelImageCombiner.hpp"
#include "MandelModel.hpp"
#include "JuliaModel.hpp"
#include "LaguerreModel.hpp"
#include "ShareableImageWrapper.hpp"

#define TRACKING_NEW_DELETE 0
  //add breakpoint in qtentrypoint_win.cpp, func qtEntryPoint(), last line (return exitCode)
  //memory_stats should be all 0 there
#if TRACKING_NEW_DELETE
struct MemoryStats
{
    int bytes;
    int blocks;
    int abytes;
    int ablocks;
} memory_stats;

constexpr size_t MEMORY_STAMP_BEGIN=0x3141592653589793;
constexpr size_t MEMORY_STAMP_END=0x2718281828459045;

void *operator new(size_t size)
{
    if (size==108)
        asm ( "nop" );
    memory_stats.bytes+=size;
    memory_stats.blocks++;
    size_t *result=(size_t *)calloc(size+3*sizeof(size_t), 1);
    result[0]=size;
    result[1]=MEMORY_STAMP_BEGIN;
    *(size_t *)((char *)result+size+2*sizeof(size_t))=MEMORY_STAMP_END;
    return result+2;
}

void operator delete(void *data) noexcept
{
    size_t *actual1=(size_t *)data;
    size_t size=actual1[-2];
    size_t stamp_begin=actual1[-1];
    if (size==0xfeeefeeefeeefeee && stamp_begin!=MEMORY_STAMP_BEGIN)
    {
      free(data);
      return;
    };
    size_t stamp_end=*(size_t *)((char *)data+size+0*sizeof(size_t));
    if (stamp_begin!=MEMORY_STAMP_BEGIN)
        asm ( "nop" );
    if (stamp_end!=MEMORY_STAMP_END)
        asm ( "nop" );
    memory_stats.bytes-=size;
    memory_stats.blocks--;
    free(actual1-2);
}

void operator delete(void *data, size_t size) noexcept
{
    size_t *actual1=(size_t *)data;
    if (actual1[-2]!=size)
        asm ( "nop" );
    size_t stamp_begin=actual1[-1];
    size_t stamp_end=*(size_t *)((char *)data+size+0*sizeof(size_t));
    if (stamp_begin!=MEMORY_STAMP_BEGIN)
        asm ( "nop" );
    if (stamp_end!=MEMORY_STAMP_END)
        asm ( "nop" );
    memory_stats.bytes-=size;
    memory_stats.blocks--;
    free(actual1-2);
}

void *operator new[](size_t size)
{
    //if (size==108)    Qt does not free argv[1] because it's nulled in QCoreApplicationPrivate::processCommandLineArguments
    //  asm ( "nop" );
    memory_stats.abytes+=size;
    memory_stats.ablocks++;
    size_t *result=(size_t *)calloc(size+3*sizeof(size_t), 1);
    result[0]=size;
    result[1]=MEMORY_STAMP_BEGIN;
    *(size_t *)((char *)result+size+2*sizeof(size_t))=MEMORY_STAMP_END;
    return result+2;
}

void operator delete[](void *data) noexcept
{
    size_t *actual1=(size_t *)data;
    size_t size=actual1[-2];
    size_t stamp_begin=actual1[-1];
    size_t stamp_end=*(size_t *)((char *)data+size+0*sizeof(size_t));
    if (stamp_begin!=MEMORY_STAMP_BEGIN)
        asm ( "nop" );
    if (stamp_end!=MEMORY_STAMP_END)
        asm ( "nop" );
    memory_stats.abytes-=size;
    memory_stats.ablocks--;
    free(actual1-2);
}

void operator delete[](void *data, size_t size) noexcept
{
    size_t *actual1=(size_t *)data;
    if (actual1[-2]!=size)
        asm ( "nop" );
    size_t stamp_begin=actual1[-1];
    size_t stamp_end=*(size_t *)((char *)data+size+0*sizeof(size_t));
    if (stamp_begin!=MEMORY_STAMP_BEGIN)
        asm ( "nop" );
    if (stamp_end!=MEMORY_STAMP_END)
        asm ( "nop" );
    memory_stats.abytes-=size;
    memory_stats.ablocks--;
    free(actual1-2);
}
#endif

int main(int argc, char *argv[])
{
  //memset(&memory_stats, 0x00, sizeof(memory_stats));
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
  QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

#if TRACKING_NEW_DELETE
  int original_argc=argc;
  char **original_argv=new char *[argc+1];
  for (int i=0; i<=argc; i++) //ends with nullptr after end, copy that too justin case
    original_argv[i]=argv[i];
#endif
  QGuiApplication app(argc, argv);

  qmlRegisterType<MandelImageCombiner>("com.example.MandelImageCombiner", 1, 0, "MandelImageCombiner");
  qmlRegisterType<MandelModel>("com.example.MandelModel", 1, 0, "MandelModel");
  qmlRegisterType<JuliaModel>("com.example.JuliaModel", 1, 0, "JuliaModel");
  qmlRegisterType<LaguerreModel>("com.example.LaguerreModel", 1, 0, "LaguerreModel");
  qRegisterMetaType<ShareableImageWrapper>();
  qRegisterMetaType<ShareableViewInfo>();

  QQmlApplicationEngine engine;
  const QUrl url(QStringLiteral("qrc:/main.qml"));
  QObject::connect(&engine, &QQmlApplicationEngine::objectCreated,
                   &app, [url](QObject *obj, const QUrl &objUrl) {
    if (!obj && url == objUrl)
      QCoreApplication::exit(-1);
  }, Qt::QueuedConnection);
  engine.load(url);

  int result=app.exec();
#if TRACKING_NEW_DELETE
  for (int i=0; i<=original_argc; i++) //ends with nullptr after end, copy that too justin case
    argv[i]=original_argv[i];
  delete[] original_argv;
#endif
  return result;
}
