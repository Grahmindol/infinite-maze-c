# Infinite Maze 🧩

Single-header C library for infinite procedurally generated mazes  
Chunk-based, deterministic, recursive DFS.

## Features 🚀
- Infinite maze (lazy generation)
- Single-header (stb-style)
- Deterministic via seed
- No dependencies (libc only)

## Usage 🔧

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
      putchar((c && !c->is_open[DIR_W]) ? '|' : ' ');  putchar(' '); putchar((c && c->is_fixed) ? ' ' : ' '); putchar(' ');
    }
    putchar('|'); putchar('\n');
  }
  for (int x = 0; x < W; x++) { putchar('+'); putchar('-'); putchar('-'); putchar('-'); }
  putchar('+'); putchar('\n');
}

int main(void) {
  maze_t* root = infinite_maze_new(0xC0FFEE);

  infinite_maze_print_ascii(root, 32, 32);

  infinite_maze_free(root);
  return 0;
}

```
