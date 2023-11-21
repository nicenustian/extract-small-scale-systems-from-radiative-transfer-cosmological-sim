CFLAGS = -Wall -fopenmp -g -C  
LIBFLAGS = -lm
OBJ = main.o sys_3d.o
EXEC = sys_3d.x

all: $(EXEC)

$(EXEC): $(OBJ)
	g++ $(CFLAGS) -o $(EXEC) $(OBJ) $(LIBFLAGS)

main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp

sys_3d.o: sys_3d.cpp
	g++ $(CFLAGS) -c sys_3d.cpp

clean:
	rm -f $(OBJ) $(EXEC)
