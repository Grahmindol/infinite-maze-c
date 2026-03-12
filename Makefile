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

wasm:
	emcc infinite_maze.c -O3 -o $(BIN)/maze.js -I$(HDR) \
	-s EXPORTED_FUNCTIONS="['_infinite_maze_new','_infinite_maze_free','_infinite_maze_is_walkable','_infinite_maze_get_cell','_infinite_path_walk']" \
	-s EXPORTED_RUNTIME_METHODS="['ccall','cwrap','addFunction','removeFunction']" \
	-s ALLOW_TABLE_GROWTH \
	-s ENVIRONMENT=web \
	-s SINGLE_FILE=1 \
	-s WASM=1 \
	-s MODULARIZE=1 \
	-s EXPORT_NAME="MazeModule"
	echo "class InfiniteMaze{#e=null;#i=null;#t=null;path=[];static async create(e=1234){const i=await MazeModule();return new InfiniteMaze(e,i)}constructor(e,i){this.#t=i,this.path=[],this.#e=i._infinite_maze_new(e),this.#i=i.addFunction(((e,i,t)=>{this.path.push({x:e,y:i})}),\"viii\")}isWalkable(e,i){if(!this.#e)throw new Error(\"Maze already freed\");return!!this.#t._infinite_maze_is_walkable(e,i,this.#e)}getPath(e,i,t,a){if(!this.#e)throw new Error(\"Maze already freed\");return this.path=[],this.#t._infinite_path_walk(e,i,t,a,this.#e,this.#i),this.path}free(){this.#e&&(this.#t._infinite_maze_free(this.#e),this.#e=null,this.#t.removeFunction(this.#i),this.#i=null)}}window.InfiniteMaze=InfiniteMaze;" >> $(BIN)/maze.js
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
