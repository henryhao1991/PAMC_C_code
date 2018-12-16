CC=gcc

all: pamc

pamc: dSFMT.o PAMC_main.o
	$(CC)  -o PAMC_main PAMC_main.o dSFMT.o -lm

dSFMT.o: dSFMT.c
	$(CC) -O  -DDSFMT_MEXP=521 -o dSFMT.o  -c dSFMT.c

PAMC_main.o: PAMC_main.c ./Lorenz96/model_L96.h
	$(CC) -O -DDSFMT_MEXP=521 -o PAMC_main.o -c PAMC_main.c

clean:
	rm *.o
