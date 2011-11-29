CC=gcc
CXX=g++

INCLUDE=-I/opt/local/include

CC_LIBS=-L/opt/local/lib -lgmp -lmpfr
CXX_LIBS=-L/opt/local/lib -lgmpxx -lgmp


all: ecpp atkin aks miller-rabin gprime

# primality tests
ecpp: ecpp.cpp aks.c aks.h miller-rabin.c miller-rabin.h
	$(CC) $(INCLUDE) $? -o ecpp-aks $(CC_LIBS)

miller-rabin: miller-rabin.c miller-rabin.h miller-rabin-driver.c
	$(CC) $(INCLUDE) $? -o miller-rabin $(CC_LIBS)

aks: aks.c aks.h aks-driver.c
	$(CC) $(INCLUDE) $? -o aks $(CC_LIBS)

# helpers
gprime: miller-rabin.c miller-rabin.h gprime.c
	$(CC) $(INCLUDE) $? -o gprime $(CC_LIBS)

# reference implementation
atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $? -o atkin $(CXX_LIBS)

clean:
	@echo Cleaning Primes Project
	rm -f *.o ecpp-aks atkin aks miller-rabin gprime

