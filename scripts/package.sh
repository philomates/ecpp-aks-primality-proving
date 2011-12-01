#!/bin/sh

mkdir -p final/program
cp README final
cp aks.h aks.cpp aks-driver.cpp ecpp.cpp gprime.cpp miller-rabin.h miller-rabin.cpp miller-rabin-driver.cpp final/program
cp Makefile-final final/program/Makefile

cd final
tar czvf ../final.tgz *
cd ..
rm -rf final

