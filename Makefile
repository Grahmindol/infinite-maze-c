NAME = infinite_maze

CC = gcc
MINGW_CC = x86_64-w64-mingw32-gcc

HDR = headers
BIN = bin
LIB = $(BIN)/lib
TEST = $(BIN)/test

SRC = $(NAME).c
TEST_SRC = test/main.c

SO  = $(LIB)/lib$(NAME).so
DLL = $(LIB)/$(NAME).dll
TEST_EXE = $(TEST)/test

all: dirs so dll

# -------- Linux SO --------
so:
	$(CC) -O3 -shared -fPIC -I$(HDR) $(SRC) -o $(SO)

# -------- Windows DLL --------
dll:
	$(MINGW_CC) -O3 -shared -municode -I$(HDR) $(SRC) -o $(DLL)

# -------- TEST --------
test: so
	$(CC) -I$(HDR) $(TEST_SRC) -L$(LIB) -l$(NAME) -o $(TEST_EXE)
	LD_LIBRARY_PATH=$(LIB) ./$(TEST_EXE)

# -------- DIRS --------
dirs:
	mkdir -p $(LIB) $(TEST)

clean:
	rm -rf bin

.PHONY: all so dll test clean dirs
