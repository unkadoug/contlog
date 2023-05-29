CFLAGS = -g

all: contlog

contlog: main.o contlog.o
	cc -o $@ -g main.o contlog.o

contlog.o: contlog.c contlog.h

main.o: main.c contlog.h

clean:
	rm contlog main.o contlog.o
