# Makefile -- StarcodeV, sequential build (plain C, gcc, NO CUDA).
#
# See CONVENTIONS_STARCODEV.md for the file layout. Main targets:
#
#   make seq        build the sequential binary ./starcodev_seq
#   make test       build, then run the built-in smoke test (5 lines)
#   make clean      remove object files and the binary
#
# This Makefile builds the sequential program only (no CUDA).

CC       ?= gcc
CFLAGS   ?= -O2 -std=gnu99 -Wall -Wextra
LDFLAGS  ?=
LDLIBS   := -lm

SRC_SEQ  := pairgen.c starcodev.c autocal.c exactcons.c output_writer.c main-starcodev.c
OBJ_SEQ  := $(SRC_SEQ:.c=.o)
BIN_SEQ  := starcodev_seq

.PHONY: all seq clean test

all: seq

seq: $(BIN_SEQ)

$(BIN_SEQ): $(OBJ_SEQ)
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Explicit header dependencies (the project is small enough that this is
# plenty without pulling in an auto-depend setup -- just enough to force a
# rebuild whenever a header changes).
pairgen.o: pairgen.c pairgen.h
starcodev.o: starcodev.c starcodev.h pairgen.h
autocal.o: autocal.c autocal.h starcodev.h pairgen.h
exactcons.o: exactcons.c exactcons.h starcodev.h
output_writer.o: output_writer.c output_writer.h starcodev.h
main-starcodev.o: main-starcodev.c starcodev.h pairgen.h autocal.h exactcons.h output_writer.h

test: seq
	@printf 'AAAACCCCGGGGTTTT\t3\nAAAACCCCGGGGTTTA\t2\nTTTTGGGGCCCCAAAA\t5\nTTTTGGGGCCCCAAAT\t1\nGATTACAGATTACAGA\t4\n' > /tmp/starcodev_smoke.txt
	./$(BIN_SEQ) -i /tmp/starcodev_smoke.txt -q
	@echo "smoke test done (you should see 3 clusters above: AAAA.../AAAA..., GATTACA..., TTTT.../TTTT...)"

clean:
	rm -f $(OBJ_SEQ) $(BIN_SEQ)
