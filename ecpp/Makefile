CXX=gcc
INCLUDE=-I/opt/local/include
LIBS=-L/opt/local/lib -lgmp

all: atkin ecpp

.PHONY: clean

utils.o: utils.cpp utils.h
	$(CXX) -c $<

atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $(LIBS) -o atkin atkin237.cpp

ecpp: ecpp.cpp utils.o
	$(CXX) $(INCLUDE) $(LIBS) -o ecpp ecpp.cpp utils.o

clean:
	@echo Cleaning ECPP Project
	rm -f *.o
