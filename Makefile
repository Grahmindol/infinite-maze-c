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
	g++ -g3 -O0 -Wall -Wextra -fsanitize=address,undefined tests/main.c  -o build/tests && ./build/tests

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
	-o build/js/infinite_maze.js

#py:
#	mkdir -p build/py
#	swig -c++ -python -outdir build/py -o build/py/wrap.cpp infinite_maze.i
#	g++ -DINFINITE_MAZE_IMPLEMENTATION \
#		build/py/wrap.cpp \
#		-I/usr/include/python3.14/ \
#		-I. \
#		-fPIC -shared \
#		-o build/py/_infinite_maze.so 

py:
	rm -r build/py
	mkdir -p build/py
	cp infinite_maze.i build/py/infinite_maze.i
	cp infinite_maze.h build/py/infinite_maze.h
	cp README.md build/py/README.md
	cp template/pyproject.toml build/py/pyproject.toml
	cp template/setup.py build/py/setup.py
	cd build/py && python -m build --sdist


clean:
	rm -rf $(BUILD_DIR)

docs: Doxyfile README.md
	doxygen Doxyfile
	cp -r demo docs/html/

.PHONY: all clean docs wasm test py