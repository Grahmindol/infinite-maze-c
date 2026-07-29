CC ?= gcc

CFLAGS ?= -Wall -Wextra -std=c99 -O3
BUILD_DIR := build

LIB := $(BUILD_DIR)/libinfinite_maze.so
HEADER := infinite_maze.h

all: $(LIB)

$(LIB): $(HEADER)
	@mkdir -p $(BUILD_DIR)
	$(CC) -x c \
		-DINFINITE_MAZE_IMPLEMENTATION \
		-shared -fPIC \
		$< \
		-o $@

test:
	gcc -g3 -O0 -Wall -Wextra -fsanitize=address,undefined  -DINFINITE_MAZE_IMPLEMENTATION tests/main.c  -o build/tests && ./build/tests

wasm:
	echo '#include "infinite_maze.h"' | emcc -x c \
	-DINFINITE_MAZE_IMPLEMENTATION \
	-O3 \
	-I. \
	-s MODULARIZE=1 \
	-s EXPORT_NAME="InfiniteMazeModule" \
	-s EXPORTED_RUNTIME_METHODS="['ccall','cwrap']" \
	-s SINGLE_FILE=1 \
	- \
	-o build/infinite_maze.js

clean:
	rm -rf $(BUILD_DIR)

docs: Doxyfile README.md
	doxygen Doxyfile
	cp -r demo docs/html/

.PHONY: all clean docs wasm test