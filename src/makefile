CC = g++
CFLAGS = -O3 --std=c++11
HEADERS = triangle.h

main.x: main.o triangle.o
	$(CC) main.o triangle.o -o main.x $(CFLAGS)

%.o: %.cpp $(HEADERS)
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm *.o

