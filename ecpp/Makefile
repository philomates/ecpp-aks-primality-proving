CC=gcc
CXX=g++
INCLUDE=-I/opt/local/include
CC_LIBS=-L/opt/local/lib -lgmp
CXX_LIBS=-L/opt/local/lib -lgmpxx -lgmp

all: ecpp atkin

.PHONY: clean

ecpp: ecpp.cpp
	$(CC) $(INCLUDE) $(CC_LIBS) -o ecpp ecpp.cpp

atkin: atkin237.cpp
	$(CXX) $(INCLUDE) $(CXX_LIBS) -o atkin atkin237.cpp

debug: ecpp.cpp
	$(CC) $(INCLUDE) $(CC_LIBS) -g -o ecpp ecpp.cpp

clean:
	@echo Cleaning ECPP Project
	rm -f *.o
