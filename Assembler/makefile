CC=g++
CFLAGS=-std=c++11
# CFLAGS=-Wall -stdlib=libc++ -std=c++11 -O3


main: main.o BWT.o
	# Link
	$(CC) $(CFLAGS) -Ofast -o main main.o BWT.o

tests: tests.o BWT.o
	# Link
	$(CC) $(CFLAGS) -I "../boost_1_66_0/" -o tests BWT.o  tests.o

main.o: main.cpp
	# Complile main.cpp
	$(CC) $(CFLAGS) -Ofast -c main.cpp

tests.o: tests.cpp
	# Complile tests.cpp
	$(CC) $(CFLAGS) -I "../boost_1_66_0/" -c tests.cpp

BWT.o: BWT.cpp
	# Compile BWT Class
	$(CC) $(CFLAGS) -Ofast -c BWT.cpp

clean:
	rm -f main.o
	rm -f BWT.o
	rm -f ./main

	rm -f tests.o
	rm -f ./tests
