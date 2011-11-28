CC=gcc
CXX=g++
INCLUDE=-I/opt/local/include
CC_LIBS=-L/opt/local/lib -lgmp
CXX_LIBS=-L/opt/local/lib -lgmpxx -lgmp

all: ecpp atkin aks miller-rabin gprime

ecpp: ecpp.cpp
	$(CC) $(INCLUDE) $(CC_LIBS) -o ecpp ecpp.cpp

atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $(CXX_LIBS) -o atkin atkin237.cpp

debug: ecpp.cpp
	$(CC) $(INCLUDE) $(CC_LIBS) -g -o ecpp ecpp.cpp

aks: aks.c
	$(CC) $(INCLUDE) $(CC_LIBS) -lmpfr -o aks aks.c

mrlib.o: mrlib.c
	$(CC) $(INCLUDE) -c mrlib.c

miller-rabin: mrlib.o miller-rabin.c
	$(CC) $(INCLUDE) $(CC_LIBS) -o miller-rabin mrlib.o miller-rabin.c

gprime: mrlib.o gprime.c
	$(CC) $(INCLUDE) $(CC_LIBS) -o gprime mrlib.o gprime.c


clean:
	@echo Cleaning Primes Project
	rm -f *.o ecpp atkin aks miller-rabin gprime
