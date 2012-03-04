#!/bin/sh

#autoscan
#cp configure.scan configure.ac
# modify configure.ac by hand

echo "Remove *.o."
rm -f *.o

echo "Execute autoheader."
autoheader

echo "Execute aclocal."
aclocal

echo "Execute automake --add-missing."
automake --add-missing

echo "Execute automake."
automake

echo "Execute autoconf."
autoconf

./configure
make
