CC=gcc

LDFLAGS=-lm
BINARY=swe

SOURCE_FILES=run_swe.c\
             shallow_water.c

OBJECTS=$(SOURCE_FILES:.c=.o)

all: $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $(BINARY)

$(OBJECTS):
	$(CC) $(CFLAGS) -o $@ -c $(@:.o=.c)

clean:
	rm $(OBJECTS) $(BINARY)
