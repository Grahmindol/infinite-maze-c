# Infinite Maze 🧩

Single-header C library for infinite procedurally generated mazes  
Chunk-based, deterministic, recursive DFS.

---

## Features 🚀

This project implements an **infinite maze generator** ♾️.

The maze is **perfect**:  
for **any two points** in the maze, there exists **one and only one path** connecting them 🧵.

### Key properties
- 🌐 Infinite, procedurally generated maze  
- 🛤️ **Unique path** between any two positions  
- 📏 The connecting path is always **finite**, even though the maze itself is infinite  
- 🌲 Structurally equivalent to an infinite spanning tree  

This guarantees **no loops**, **no ambiguity**, and deterministic navigation throughout the maze.


## Proof of Concept ⚠️

This project is a **proof of concept (PoC)** 🧪.

It demonstrates that the core idea **works correctly**, but it is **not yet optimized** 🚧.  
Performance, memory usage, and edge cases have **not** been the main focus so far.

### What this means
- ✅ Correctness over performance  
- ❌ Not production-ready  
- 🔧 Further **optimizations and refactoring are required**

### Next steps
- I found a limitation in the current maze generation algorithm: every generated maze has its main dead end located at the center, because the generation process always starts from the center cell.

This reduces the local variety of generated mazes, as there are effectively fewer than 32 distinct maze patterns available around a given area.

However, this discovery also reveals an opportunity for a major optimization. Since the number of possible local patterns is limited, it will be possible to precompute and store the 32 possible patterns, then refactor the data structure to reuse them efficiently.

A future improvement is planned to make this behavior configurable, allowing users to choose between the current procedural generation approach and the optimized pattern-based generation mode.

Contributions and optimization suggestions are **welcome** 🙌


## Demo :
  Try the interactive WebAssembly version:

- [Browser Demo](https://grahmindol.github.io/infinite-maze-c/html/demo/index.html)

## Doxygen page :
- [Documentation !!](https://grahmindol.github.io/infinite-maze-c/html/infinite__maze_8h.html)

## 🧪 Example

```c
#include <stdio.h>
#include <stdlib.h>

#define MAZE_RADIUS 2
#define INFINITE_MAZE_IMPLEMENTATION
#include "infinite_maze.h"

enum biom_t {
  WALKABLE  = 0b001,
  DEAD_END  = 0b010,
  TINY_ROOM = 0b100,
};

/** ASCII print for testing */
void infinite_maze_print_ascii(void* root, int W, int H) {
  int minx = -W / 2, miny = -H / 2;
  for (int y = H - 1; y >= 0; y--) {
    for (int x = 0; x < W; x++) {
      char l = infinite_maze_get_cell(minx + x, miny + y, root);
      putchar(l & WALKABLE ? (l & DEAD_END ? '*' : ' ') : 219);
      putchar(l & WALKABLE ? (l & TINY_ROOM ? '/' : ' ') : 219);
    }
     putchar('\n');
  }
  putchar('\n');
}

int main(void) {
  void* root = infinite_maze_new(0xC0FFEE);

  infinite_maze_print_ascii(root, 64, 64);

  infinite_maze_free(root);
  
  return 0;
}
```
