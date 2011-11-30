CC=gcc
CXX=g++

INCLUDE=-I/opt/local/include -O3

CC_LIBS=-L/opt/local/lib -lgmp -lmpfr
CXX_LIBS=-L/opt/local/lib -lgmpxx -lgmp

all: run atkin aks miller-rabin gprime

# primality tests
run: ecpp.cpp aks.o miller-rabin.o
	$(CXX) $(INCLUDE) $? -o run $(CC_LIBS)

miller-rabin:  miller-rabin-driver.c miller-rabin.o
	$(CC) $(INCLUDE) $? -o miller-rabin $(CC_LIBS)

aks: aks-driver.c aks.o
	$(CC) $(INCLUDE) $? -o aks $(CC_LIBS)

# helpers
gprime: gprime.c miller-rabin.o
	$(CC) $(INCLUDE) $? -o gprime $(CC_LIBS)

# libs
aks.o: aks.c aks.h
	$(CC) $(INCLUDE) $? -c

miller-rabin.o: miller-rabin.c miller-rabin.c
	$(CC) $(INCLUDE) $? -c

# reference implementation
atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $? -o atkin $(CXX_LIBS)

clean:
	@echo Cleaning Primes Project
	rm -f *.o run atkin aks miller-rabin gprime

