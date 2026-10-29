#!/bin/sh

if [ "$1" = "run" ]; then
  cd build-gui
  make -j$(nproc)
  if [ $? -ne 0 ]; then
    echo "Failed to build GUI"
    cd ..
    exit 1
  fi
  ./Mandelbrot
  cd ..
  exit 0
fi

if [ "$1" = "runp" ]; then
  make -f Makefile-headless -j$(nproc)
  if [ $? -ne 0 ]; then
    echo "Failed to build headless"
    exit 1
  fi
  ./mandel_runner --point --fractal jul --precision double --cx -0.09991404543315241 --cy 0.9084110242318193
  exit 0
fi

if [ "$1" = "rund" ]; then
  make -B -f Makefile-headless -j$(nproc) OPTFLAGS=-O0
  if [ $? -ne 0 ]; then
    echo "Failed to build headless"
    exit 1
  fi
  ./mandel_runner --point --fractal jul --precision double --cx -0.09991404543315241 --cy 0.9084110242318193
  exit 0
fi

if [ "$1" = "check" ]; then
  make -f Makefile-headless check-all
  exit $?
fi

if [ "$1" = "bless" ]; then
  make -f Makefile-headless bless-all
  exit $?
fi

echo "Usage: $0 [run|runp|check|bless]"
exit 1
