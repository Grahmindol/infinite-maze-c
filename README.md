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

---

## Proof of Concept ⚠️

This project is a **proof of concept (PoC)** 🧪.

It demonstrates that the core idea **works correctly**, but it is **not yet optimized** 🚧.  
Performance, memory usage, and edge cases have **not** been the main focus so far.

### What this means
- ✅ Correctness over performance  
- ❌ Not production-ready  
- 🔧 Further **optimizations and refactoring are required**

### Next steps
- Improve algorithmic efficiency ⚡  
- Reduce memory footprint 🧠  
- Add benchmarks and profiling 📊  

Contributions and optimization suggestions are **welcome** 🙌

---

## Usage 🔧

### 📚 API

1. `maze_t* infinite_maze_new(uint64_t seed)` -- Creates a new infinite maze instance.
  - 🌱 `seed` is used to initialize the procedural generation
  - ♾️ The maze is generated lazily, cells are created on demand
  - 🔁 Same seed ⇒ same maze structure
Returns a pointer to the root maze structure.


2. `void infinite_maze_free(maze_t* root)` -- Frees **all memory** associated with the maze.
  - 🧹 Recursively deallocates all generated chunks and nodes
  - ⚠️ Invalidate all `node_t*` generated

3. `node_t* infinite_maze_get_cell(int x, int y, maze_t* root)` -- Returns the cell at coordinates `(x, y)`.
- 📍 Coordinates are **unbounded** (the maze is infinite in both negative and positive coordinate)
- 🧠 Cells are generated on demand if they do not exist yet
Returns `NULL` only if allocation fails.

---

### 🧪 Example

```c
#include <stdio.h>
#include <stdlib.h>

#define MAZE_RADIUS 2
#define INFINITE_MAZE_IMPLEMENTATION
#include "infinite_maze.h"

/** ASCII print for testing */
void infinite_maze_print_ascii(maze_t* root, int W, int H) {
  int minx = -W / 2, miny = -H / 2;
  for (int y = H - 1; y >= 0; y--) {
    for (int x = 0; x < W; x++) {
      node_t* c = infinite_maze_get_cell(minx + x, miny + y, root);
      putchar('+');
      putchar((c && !c->is_open[DIR_N]) ? '-' : ' ');
      putchar((c && !c->is_open[DIR_N]) ? '-' : ' ');
      putchar((c && !c->is_open[DIR_N]) ? '-' : ' ');
    }
    putchar('+'); putchar('\n');
    for (int x = 0; x < W; x++) {
      node_t* c = infinite_maze_get_cell(minx + x, miny + y, root);
      putchar((c && !c->is_open[DIR_W]) ? '|' : ' ');
      putchar(' '); putchar(' '); putchar(' ');
    }
    putchar('|'); putchar('\n');
  }
  for (int x = 0; x < W; x++) {
    putchar('+'); putchar('-'); putchar('-'); putchar('-');
  }
  putchar('+'); putchar('\n');
}

int main(void) {
  maze_t* root = infinite_maze_new(0xC0FFEE);

  infinite_maze_print_ascii(root, 32, 32);

  infinite_maze_free(root);
  return 0;
}
