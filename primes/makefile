
all: aks miller-rabin

aks: aks.c
	gcc -I/opt/local/include -L/opt/local/lib -lmpfr -lgmp -o aks aks.c

miller-rabin: miller-rabin.c
	gcc -I/opt/local/include -L/opt/local/lib -lmpfr -lgmp -o miller-rabin miller-rabin.c

