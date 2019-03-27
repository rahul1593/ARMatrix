#!/bin/sh
INC="../inc/"
SRC="../src/"
LIB="-lm"
CMODE="-c"
COMPILER="gcc"
DEF="MODE_CPP"

if [ "x"$1 = "x"$CMODE ]
then
    SRC=$SRC'*.c'
    DEF="MODE_C"
else
    SRC=$SRC'*.cpp'
    COMPILER="g++"
fi

#build library with test code
$COMPILER -o testmtx -I $INC $SRC $LIB -O3 -D$DEF
echo "Build complete"
