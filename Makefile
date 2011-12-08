    # Copyright (C) 2011  Akshay Srinivasan 
    # <akshay@ncbs.res.in>

    # This program is free software: you can redistribute it and/or modify
    # it under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with this program.  If not, see <http://www.gnu.org/licenses/>.

DEBUG = #-g
FCFLAGS = $(DEBUG) -fPIC -std=legacy -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
FC = gfortran $(FCFLAGS)
LINFLAGS = -lm -lsnopt_c -lf2c -lblas -lsnprint_c -lgfortran -lsnfwrap_c -ladolc -lgsl -lgslcblas -l:libclmech.a -l:libopkd.a -lglut -lnlopt_cxx -lGL -lGLU
INCFLAGS = -I optcon -I /home/neptune/include -L /home/neptune/lib -I solvers/snopt -I solvers/ipopt -I. -I optcon/orthopoly
CCFLAGS = $(DEBUG) -fPIC -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
CPPFLAGS = $(DEBUG) -fPIC -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
CC = gcc $(CCFLAGS) $(INCFLAGS)
CPP = g++ $(CPPFLAGS) $(INCFLAGS)

OBJECTS = dlsode.o utility.o main.o eom_double.o

PROG = dpend_dlsode

.f.o:
	$(FC) -c $< -o $@

.c.o:
	$(CC) -c $< -o $@

.cpp.o:
	$(CPP) -c $< -o $@

main: $(OBJECTS)
	$(CC) $(OBJECTS) -o $(PROG) $(LINFLAGS)

clean:
	find -iname "*.o" -exec rm '{}' \;