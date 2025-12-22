OPTIONS=-O3 --std=c++17

all: main

debug: OPTIONS=-O0 -g
debug: all

build/libsais.o: libsais/src/libsais.c
	gcc -O3 -c -I libsais/include/ ./libsais/src/libsais.c -o build/libsais.o

main: build/libsais.o src/main.cpp
	g++ ${OPTIONS} -I libsais/include/ src/main.cpp build/libsais.o -o build/smr

clean:
	rm -rf ./build/*
