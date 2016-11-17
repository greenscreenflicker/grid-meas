CC=cpp
CFLAGS=-lm -I.
DEPS = 
OBJ = main.o  fir1.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

PMU: $(OBJ)
	g++ -o $@ $^ $(CFLAGS)
