OBJECTS = config.o midas2smol.o midas-format.o grif-format.o \
          reorder.o default_sort.o

CFLAGS  = -g -O2

midas2smol: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) -rdynamic -lz -ldl -lm -lpthread

.c.o:
	$(CC) -c $(CFLAGS) $<

config.o:       config.c config.h midas2smol.h
grif-format.o:  grif-format.c midas2smol.h grif-format.h midas-format.h
midas-format.o: midas-format.c midas2smol.h midas-format.h
default_sort.o: default_sort.c config.h grif-format.h
midas2smol.o:  midas2smol.c config.h grif-format.h midas-format.h
reorder.o:      reorder.c midas2smol.h midas-format.h

#SOURCES = config.c midas2smol.c midas-format.c grif-format.c \
            reorder.c default_sort.c

clean:
	rm -f *.o midas2smol

# if there is a file called "clean", above will fail without this ...
.PHONY: clean
