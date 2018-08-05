CC=g++
CFLAGS=-std=c++11 -pthread -Ofast
# CFLAGS=-Wall -stdlib=libc++ -std=c++11 -O3


main: main.o ProblemDescription.o ConstructSA.o LinearSA.o
	# Link
	$(CC) $(CFLAGS) -o main ProblemDescription.o LinearSA.o ConstructSA.o  main.o

tests: tests.o ProblemDescription.o ConstructSA.o LinearSA.o
	# Link
	$(CC) $(CFLAGS) -I "../boost_1_66_0/" -o tests ProblemDescription.o LinearSA.o ConstructSA.o  tests.o

main.o: main.cpp
	# Complile main.cpp
	$(CC) $(CFLAGS) -c main.cpp

tests.o: tests.cpp
	# Complile tests.cpp
	$(CC) $(CFLAGS) -I "../boost_1_66_0/" -c tests.cpp

ProblemDescription.o: ProblemDescription.cpp
	# Compile Problem Description Class
	$(CC) $(CFLAGS) -c ProblemDescription.cpp

ConstructSA.o: ConstructSA.cpp
	# Compile ConstructSA Class
	$(CC) $(CFLAGS) -c ConstructSA.cpp


LinearSA.o: LinearSA.cpp
	# Compile LinearSA Class
	$(CC) $(CFLAGS) -c LinearSA.cpp

CompressedBWT.o: CompressedBWT.cpp
	# Compile CompressedBWT Class
	$(CC) $(CFLAGS) -c CompressedBWT.cpp

clean:
	rm -f main.o
	rm -f ConstructSA.o
	rm -f LinearSA.o
	rm -f ProblemDescription.o
	rm -f CompressedBWT.o
	rm -f ./main

	rm -f tests.o
	rm -f ./tests
