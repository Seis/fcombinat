CFLAGS = -g -D_COMPILE_LIB -Iinclude -Iinclude/grammar
LDFLAGS = -lflint -lm

GRAMMAR_SRC = src/grammar/absyn.c src/grammar/node.c src/grammar/object.c src/grammar/utils.c src/grammar/parser.tab.c src/grammar/lex.yy.c
SOLVER_SRC = src/solver/math.c \
             $(wildcard src/solver/count/*.c) \
             $(wildcard src/solver/rank/*.c) \
             $(wildcard src/solver/unrank/*.c) \
             $(wildcard src/solver/draw/*.c) \
             $(wildcard src/solver/boltzmann/*.c)
MAIN_SRC = src/api.c src/cli.c

SRC = $(MAIN_SRC) $(GRAMMAR_SRC) $(SOLVER_SRC)
OBJ = $(SRC:.c=.o)

TARGET = fcombinat

.PHONY: all test test-all clean

include config.mk
include Makefile.figures

all: $(TARGET)

src/grammar/parser.tab.c include/grammar/parser.tab.h: src/grammar/parser.y
	@bison -d $< --defines=include/grammar/parser.tab.h -o src/grammar/parser.tab.c

src/grammar/lex.yy.c: src/grammar/lexer.l include/grammar/parser.tab.h
	@flex -o $@ $<

%.o: %.c include/grammar/parser.tab.h
	@$(CC) $(CFLAGS) -c -o $@ $<

$(TARGET): $(OBJ)
	@$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)



test:
	$(MAKE) -C exp/correctness FCOMBINAT_BIN="$(abspath $(FCOMBINAT_BIN))" test
	$(MAKE) -C exp/correctness FCOMBINAT_BIN="$(abspath $(FCOMBINAT_BIN))" test-bijection
	$(MAKE) -C exp/correctness FCOMBINAT_BIN="$(abspath $(FCOMBINAT_BIN))" test-unique

test-all:
	$(MAKE) -C exp/correctness FCOMBINAT_BIN="$(abspath $(FCOMBINAT_BIN))" all

clean:
	@rm -f src/*.o src/grammar/*.o src/solver/*.o src/solver/rank/*.o src/solver/unrank/*.o src/solver/count/*.o src/solver/draw/*.o src/solver/boltzmann/*.o src/grammar/parser.tab.c include/grammar/parser.tab.h src/grammar/lex.yy.c $(TARGET)