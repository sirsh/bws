#SIMPLE MAKE with or without -w
all:
	gcc -lm -D__USE_MINGW_ANSI_STDIO .\bws.c .\lattice_core.c .\rng\mt19937ar_clean_bkp.c  -O3 -o bws.exe




