
all: aks miller-rabin gprime

aks: aks.c
	gcc -I/opt/local/include -L/opt/local/lib -lmpfr -lgmp -o aks aks.c

mrlib.o: mrlib.c
	gcc -I/opt/local/include -c mrlib.c

miller-rabin: mrlib.o miller-rabin.c
	gcc -I/opt/local/include -L/opt/local/lib -lgmp -o miller-rabin mrlib.o miller-rabin.c

gprime: mrlib.o gprime.c
	gcc -I/opt/local/include -L/opt/local/lib -lgmp -o gprime mrlib.o gprime.c

clean:
	rm -f mrlib.o
	rm -f aks miller-rabin gprime

