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

FCFLAGS = -fPIC -std=legacy -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
FC = gfortran $(FCFLAGS)
LINFLAGS = -lm -lgfortran
INCFLAGS =
CCFLAGS = -fPIC -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
CC = gcc $(CCFLAGS) $(INCFLAGS)

OBJECTS = dlsode.o utility.o

LIBNAME = cdlsode

.f.o:
	$(FC) -c $< -o $@

.c.o:
	$(CC) -c $< -o $@

static: $(OBJECTS)
	ar -cvq lib$(LIBNAME).a $(OBJECTS)

dynamic: $(OBJECTS)
	$(CC) $(LINFLAGS) -shared -Wl,-soname,lib$(LIBNAME).so -o lib$(LIBNAME).so $(OBJECTS)