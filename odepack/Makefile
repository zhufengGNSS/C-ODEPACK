FCFLAGS = -fPIC -std=legacy -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
FC = gfortran $(FCFLAGS)
LINFLAGS = -lm -lgfortran
INCFLAGS =
CCFLAGS = -fPIC -O3 -ftree-vectorize -msse2 -ftree-vectorizer-verbose=0 -ffast-math
CC = gcc $(CCFLAGS) $(INCFLAGS)

DOBJECTS = opkdmain.o opkda1.o opkda2.o
SOBJECTS = opksmain.o opksa1.o opksa2.o

DLIBNAME = opkd
SLIBNAME = opks

.f.o:
	$(FC) -c $< -o $@

.c.o:
	$(CC) -c $< -o $@

all:
	make double-dynamic
	make double-static
	make single-dynamic
	make single-static

double-dynamic: $(DOBJECTS)
	$(CC) $(LINFLAGS) -shared -Wl,-soname,lib$(DLIBNAME).so -o lib$(DLIBNAME).so $(DOBJECTS)

double-static: $(DOBJECTS)
	ar -cvq lib$(DLIBNAME).a $(DOBJECTS)

single-dynamic: $(SOBJECTS)
	$(CC) $(LINFLAGS) -shared -Wl,-soname,lib$(SLIBNAME).so -o lib$(SLIBNAME).so $(SOBJECTS)

single-static: $(SOBJECTS)
	ar -cvq lib$(SLIBNAME).a $(SOBJECTS)

clean:
	rm -f *.o *~