CC = gcc
CFLAGS = -W -Wall 
TARGET = hw3

$(TARGET) : hw3.o bessj0.o bessj1.o nrutil.o rtbis.o rtflsp.o rtnewt.o rtsafe.o rtsec.o zbrak.o
			 $(CC) $(CFLAGS) -o $@ $^ -lm

hw3.o : hw3.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

bessj0.o : bessj0.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

bessj1.o : bessj1.c
		$(CC) $(CFLAGS) -c -o $@ $^ -lm

nrutil.o : nrutil.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtbis.o : rtbis.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtflsp.o : rtflsp.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtnewt.o : rtnewt.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtsafe.o : rtsafe.c
		$(CC) $(CFLAGS) -c -o $@ $^

rtsec.o : rtsec.c
		$(CC) $(CFLAGS) -c -o $@ $^

zbrak.o : zbrak.c
		$(CC) $(CFLAGS) -c -o $@ $^

clean :
		rm *.o 
