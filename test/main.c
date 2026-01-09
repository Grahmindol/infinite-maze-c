#include <stdio.h>
#include <stdlib.h>

#include "infinite_maze.h"

/** ASCII print for testing */
void infinite_maze_print_ascii(void* root, int W, int H) {
  int minx = -W / 2, miny = -H / 2;
  for (int y = H - 1; y >= 0; y--) {
    for (int x = 0; x < W; x++) {
      char l = infinite_maze_get_cell(minx + x, miny + y, root);
      putchar(l&1 ? (l&2 ? '*' : ' ') : 219);
      putchar(l&1 ? (l&4 ? '/' : ' ') : 219);
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