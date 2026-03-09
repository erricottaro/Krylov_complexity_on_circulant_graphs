CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = triangle.h

main: main.o triangle.o
	$(CC) main.o triangle.o -o main $(CFLAGS)

%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm *.o

