CXX=gcc
INCLUDE=-I/opt/local/include
LIBS=-L/opt/local/lib -lgmp

all: ecpp atkin

.PHONY: clean

utils.o: utils.cpp utils.h
	$(CXX) -c $<

ecpp: ecpp.cpp utils.o
	$(CXX) $(INCLUDE) $(LIBS) -o ecpp ecpp.cpp utils.o

atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $(LIBS) -o atkin atkin237.cpp

debug: ecpp.cpp utils.o
	$(CXX) $(INCLUDE) $(LIBS) -g -o ecpp ecpp.cpp utils.o

clean:
	@echo Cleaning ECPP Project
	rm -f *.o
