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

## 🧩 Public API :

- ### infinite_maze_new
  
  ```c
  void* infinite_maze_new(int seed);
  ```
  
  Create a new infinite procedurally generated maze.
  
  - Allocates and initializes the root chunk.
  - Maze generation is **deterministic** with respect to `seed`.
  - Additional chunks are generated **lazily** on access.
  
  **Parameters**
  - `seed` — Initial random seed controlling maze generation.
  
  **Returns**
  - Opaque pointer to the maze instance.
  - `NULL` on allocation failure.
  
  **Complexity**
  - **O(1)**

- ### infinite_maze_free
  
  ```c
  void infinite_maze_free(void* m);
  ```
  
  Free an infinite procedurally generated maze.
  
  - Recursively frees all generated chunks.
  - Includes inner and outer mazes.
  
  **Parameters**
  - `m` — Maze instance returned by `infinite_maze_new()`.
  
  **Complexity**
  - **O(log(r))**, where *r* is the distance from the origin to the farthest generated cell.
  
  **Warning**
  - All pointers obtained from this maze become invalid.

- ### infinite_maze_is_walkable
  
  ```c
  bool infinite_maze_is_walkable(int wx, int wy, void* maze_p);
  ```
  
  Test whether a world-space cell is walkable.
  
  - Returns whether the cell is an open passage or a wall.
  - May trigger lazy chunk generation.
  
  **Parameters**
  - `wx` — World X coordinate.
  - `wy` — World Y coordinate.
  - `maze_p` — Maze instance.
  
  **Returns**
  - `true` if walkable
  - `false` otherwise
  
  **Complexity**
  - Average: **O(1)**
  - Worst case: **O(log(wx² + wy²))**


- ### infinite_maze_get_cell

  ```c
  uint8_t infinite_maze_get_cell(int wx, int wy, void* maze_p);
  ```

  Retrieve hierarchical dead-end information for a world cell.

  Each bit of the returned byte encodes maze information:

  - Bit 0 : Walkability (1 = walkable, 0 = wall)
  - Bit 1 : Dead-end at local chunk level
  - Bits 2–7 : Dead-end status propagated through parent maze levels

  **Parameters**
  - `wx` — World-space X coordinate.
  - `wy` — World-space Y coordinate.
  - `maze_p` — Pointer to the root maze instance.

  **Returns**
  - Encoded walkability and dead-end hierarchy information
  - `0` if `maze_p` is NULL

  **Complexity**
  - **O(log(wx² + wy²))**

  **Notes**
  - May trigger lazy generation of parent maze chunks.
  - See `infinite_maze_is_walkable`
---

### 🧪 Example

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
