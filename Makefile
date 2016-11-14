CC    = g++
CCFLAGS = -O3 -W -Wall -g

MSRCS   = Fambl.c
SRCS	= Read.c Family.c Classify.c Common.c Metrics.c
OBJS	= $(SRCS:.c=.o)

Fambl: Fambl.o Fambllib.a
	$(CC) $(LDFLAGS) $^ -lm -o $@

%.o: %.c
	$(CC) -c $(CCFLAGS) -o $@ $<

Fambllib.a: $(OBJS)
	ar ruv $@ $?
	ranlib $@

clean:
	-rm Fambl
	-rm *.o
	-rm Fambllib.a
