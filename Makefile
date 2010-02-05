
CC = g++ -Wall
LIB = -lntl
HDRS = -I .
OBJS = polynomial_selection.o sieve.o linear_algebra.o square_root.o

default: $(OBJS)
	$(CC) $(HDRS) prime_gen.cpp -o primegen $(LIB)
	$(CC) $(HDRS) n_gen.cpp -o ngen $(LIB)
	$(CC) $(HDRS) $(OBJS) gnfs.cpp -lntl -o gnfs $(LIB)


%.o: %.cpp %.hpp
	$(CC) $(HDRS) -c $< -o $@

clean:
	rm -f *.o
	rm -f primegen
	rm -f gnfs
	rm -f ngen

