#CC = icpc
CC = g++
#CFLAGS = -m64 -openmp
CFLAGS = -O0 -g
#CFLAGS = -O3

#LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -ldl -lpthread -lfftw3 -lm
#LIBS = -lmkl_intel_lp64 -lfftw3 -lm -ltr
#LIBS = -lmkl_intel_lp64 -lmkl_core -ldl -lfftw3 -lm

#LIBDIR =/lib64
#LIBS=-L${LIBDIR} -lblas -lgfortran -lfftw3 -lm
LIBS=-lblas -lgfortran -lfftw3 -lm

all: aim

aim: AIM.C timer.C
	$(CC)  $(CFLAGS) $(CPPFLAGS) $(LFLAGS) -o aim *.C $(LIBS)

clean:
	$(RM) aim *.info

