#!/bin/sh
INC="../inc/"
SRC="../src/*"
LIB="-lm"

#build library with test code
g++ -o testmtx -I $INC $SRC $LIB -O3
echo "Build complete"
