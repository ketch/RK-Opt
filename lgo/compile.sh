#!/bin/sh

gfortran -c -Wall LGOMAIN.F90 > messages.txt
gfortran -c -Wall LGO.FOR >> messages.txt
gfortran LGOMAIN.O LGO.O -o LGO.EXE >> messages.txt
