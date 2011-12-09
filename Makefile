CXX=g++
INCLUDE=-I/opt/local/include -O3
CXX_LIBS=-L/opt/local/lib -lgmpxx -lmpfr -lgmp

all: run aks miller-rabin gprime hybrid

# primality tests
run: ecpp.cpp aks.o miller-rabin.o
	$(CXX) $(INCLUDE) $^ -o run $(CXX_LIBS)

hybrid: ecpp-to-aks.cpp aks.o miller-rabin.o
	$(CXX) $(INCLUDE) $^ -o hybrid $(CXX_LIBS)

miller-rabin: miller-rabin-driver.cpp miller-rabin.o
	$(CXX) $(INCLUDE) $^ -o miller-rabin $(CXX_LIBS)

aks: aks-driver.cpp aks.o
	$(CXX) $(INCLUDE) $^ -o aks $(CXX_LIBS)

# helpers
gprime: gprime.cpp miller-rabin.o
	$(CXX) $(INCLUDE) $^ -o gprime $(CXX_LIBS)

# libs
aks.o: aks.cpp aks.h
	$(CXX) $(INCLUDE) $^ -c

miller-rabin.o: miller-rabin.cpp miller-rabin.h
	$(CXX) $(INCLUDE) $^ -c

tarball: run
	mkdir program
	cp Makefile *.cpp *.h run program
	tar czvf primes.tgz README program
	rm -rf program

clean:
	@echo Cleaning Primes Project
	rm -f *.o *.gch run aks miller-rabin gprime

