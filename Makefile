CC=gcc
CXX=g++

INCLUDE=-I/opt/local/include

CC_LIBS=-L/opt/local/lib -lmpfr -lgmp
CXX_LIBS=-L/opt/local/lib -lgmpxx -lgmp


all: ecpp atkin aks miller-rabin gprime

# primality tests
ecpp: ecpp.o aks.o
	$(CC) -o ecpp aks.o ecpp.o $(INCLUDE) $(CC_LIBS)

ecpp.o: ecpp.cpp aks.h
	$(CC) -c ecpp.cpp

miller-rabin: mrlib.o miller-rabin.c
	$(CC) -o miller-rabin mrlib.o miller-rabin.c $(INCLUDE) $(CC_LIBS)

aks: aks.o aks_driver.o
	$(CC) -o aks aks_driver.o aks.o $(INCLUDE) $(CC_LIBS)

aks.o: aks.c aks.h
	$(CC) $(CC_LIBS) -c $<

aks_driver.o: aks_driver.c aks.h
	$(CC) -c $<

# helpers
gprime: mrlib.o gprime.c
	$(CC) $(INCLUDE) $(CC_LIBS) -o gprime mrlib.o gprime.c

mrlib.o: mrlib.c
	$(CC) $(INCLUDE) -c mrlib.c

# reference implementation
atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $(CXX_LIBS) -o atkin atkin237.cpp

clean:
	@echo Cleaning Primes Project
	rm -f *.o ecpp atkin aks miller-rabin gprime

