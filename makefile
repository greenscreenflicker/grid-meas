CC=cpp
CFLAGS=-lm -I.
DEPS = 
OBJ = main.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

PMU: $(OBJ)
	g++ -o $@ $^ $(CFLAGS)
