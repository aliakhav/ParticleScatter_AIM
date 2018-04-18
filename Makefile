CC = g++
#CFLAGS = -O0 -g
CFLAGS = -O3

#LIBDIR =/lib64
#LIBS=-L${LIBDIR} -lblas -lgfortran -lfftw3 -lm
LIBS=-lgfortran -lfftw3 -lm

all: aim

aim: AIM.C timer.C
	$(CC)  $(CFLAGS) $(CPPFLAGS) $(LFLAGS) -o aim *.C $(LIBS)

clean:
	$(RM) aim *.info

