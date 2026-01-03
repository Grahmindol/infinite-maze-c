#ifndef INFINITE_MAZE_H
#define INFINITE_MAZE_H

/**
 * Infinite procedurally generated maze
 * Single-header library
 * 
 * Usage :
 * 
 * #define INFINITE_MAZE_IMPLEMENTATION // in ONE c file
 * #include "infinite_maze.h"
 * 
 */

#ifndef MAZE_RADIUS
#define MAZE_RADIUS 2
#endif

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* =======================
   ======= API ===========
   ======================= */

typedef enum {
    DIR_N = 0,
    DIR_W = 1,
    DIR_E = 2,
    DIR_S = 3,
    DIR_COUNT
} direction_t;

struct maze_t;

typedef struct node_t {
    direction_t parent;
    bool is_open[DIR_COUNT];
    bool is_fixed;
    struct maze_t* inner_maze;
    struct maze_t* outer_maze;
} node_t;

typedef struct maze_t {
    int seed;
    node_t* data;
    node_t* outer_node;
} maze_t;

maze_t* infinite_maze_new(int seed);
void infinite_maze_free(maze_t* m);
node_t* infinite_maze_get_cell(int wx, int wy, maze_t* root);

/* =======================
   === IMPLEMENTATION ====
   ======================= */

#ifdef INFINITE_MAZE_IMPLEMENTATION

static const int DIR_DX[DIR_COUNT] = {0, -1, 1, 0};
static const int DIR_DY[DIR_COUNT] = {1, 0, 0, -1};
static const direction_t DIR_OPP[DIR_COUNT] = {DIR_S, DIR_E, DIR_W, DIR_N};

/** Internal helpers */

static inline node_t* _node(maze_t* m, int x, int y) {
  if (!m) return NULL;
  int size = 2 * MAZE_RADIUS + 1;
  int ix = x + MAZE_RADIUS;
  int iy = y + MAZE_RADIUS;
  if (ix < 0 || ix >= size || iy < 0 || iy >= size) return NULL;
  return &m->data[iy * size + ix];
}

static inline void _carve(maze_t* m, int x, int y, direction_t d) {
  node_t* a = _node(m, x, y);
  node_t* b = _node(m, x + DIR_DX[d], y + DIR_DY[d]);
  if (!a || !b) return;
  a->is_open[d] = true;
  b->is_open[DIR_OPP[d]] = true;
  b->parent = DIR_OPP[d];
}

static void _explore(maze_t* m, int x, int y) {
  node_t* c = _node(m, x, y);
  if (!c || c->outer_maze != NULL) return;
  c->outer_maze = m;

  direction_t start = rand() % DIR_COUNT;
  for (int i = 0; i < DIR_COUNT; i++) {
    direction_t d = (start + i) % DIR_COUNT;
    int nx = x + DIR_DX[d];
    int ny = y + DIR_DY[d];
    node_t* n = _node(m, nx, ny);
    if (n && n->outer_maze == NULL) {
      _carve(m, x, y, d);
      _explore(m, nx, ny);
    }
  }
}

static void _fix_parent_path(maze_t* m, int lx, int ly){
  node_t* n = _node(m, lx, ly);
  n->is_fixed = true;
  if (n->parent < 0 || n->parent >= DIR_COUNT) return;

  lx += DIR_DX[n->parent];
  ly += DIR_DY[n->parent];
  _fix_parent_path(m, lx,ly);
}

static void _update_chunk_aperture(maze_t* m) {
  if (!m || !m->outer_node) return;

  for (int d = 0; d < DIR_COUNT; d++) {
    int dx = MAZE_RADIUS * DIR_DX[d];
    int dy = MAZE_RADIUS * DIR_DY[d];
    node_t* n = _node(m, dx, dy);
    if(n->is_open[d] = m->outer_node->is_open[d]) 
    _fix_parent_path(m,dx,dy);
    
    if (n->inner_maze) _update_chunk_aperture(n->inner_maze);
  }
}

/** Get cell in world coordinates with chunk generation */
static node_t* _raw_get_cell(int wx, int wy, maze_t* root) {
  if (abs(wx) <= MAZE_RADIUS && abs(wy) <= MAZE_RADIUS)
    return _node(root, wx, wy);

  const int size = 2 * MAZE_RADIUS + 1;
  int sgn_wx = (wx > 0) - (wx < 0);
  int cx = (wx + sgn_wx * MAZE_RADIUS) / size;
  int sgn_wy = (wy > 0) - (wy < 0);
  int cy = (wy + sgn_wy * MAZE_RADIUS) / size;

  int lx = wx - cx * size;
  int ly = wy - cy * size;

  if (!root->outer_node) {
    maze_t* outer_maze = infinite_maze_new(root->seed + 1);
    root->outer_node = _node(outer_maze, 0, 0);
    root->outer_node->inner_maze = root;
    _update_chunk_aperture(root);
  }

  maze_t* outer_maze = root->outer_node->outer_maze;
  node_t* outer_node = _raw_get_cell(cx, cy, outer_maze);

  if (!outer_node->inner_maze) {
    outer_node->inner_maze = infinite_maze_new(
        outer_maze->seed ^ (uint64_t)(cx) * 0x9E3779B185EBCA87ULL ^
        (uint64_t)(cy) * 0xC2B2AE3D27D4EB4FULL);
    outer_node->inner_maze->outer_node = outer_node;
    _update_chunk_aperture(outer_node->inner_maze);
  }

  return _node(outer_node->inner_maze, lx, ly);
}

static node_t* _generate_final_cell(int wx, int wy, maze_t* root) {
  const int size = 2 * MAZE_RADIUS + 1;
  int sgn_wx = (wx > 0) - (wx < 0);
  int cx = (wx + sgn_wx * MAZE_RADIUS) / size;
  int sgn_wy = (wy > 0) - (wy < 0);
  int cy = (wy + sgn_wy * MAZE_RADIUS) / size;

  int lx = wx - cx * size;
  int ly = wy - cy * size;

  node_t* base_node = _raw_get_cell(wx, wy, root);
  if (abs(lx) != MAZE_RADIUS && abs(ly) != MAZE_RADIUS || base_node->is_fixed)
    return base_node;

  srand(base_node->outer_maze->seed ^ (lx) * 0x9E3779B185EBCA87ULL ^
        (ly) * 0xC2B2AE3D27D4EB4FULL);
  if (lx == -MAZE_RADIUS && ly != 0 && rand() & 0b10) {
    base_node->is_open[DIR_W] = true;
    _raw_get_cell(wx - 1, wy, root)->is_open[DIR_E] = true;

    base_node->is_open[base_node->parent] = false;
    _raw_get_cell(wx + DIR_DX[base_node->parent],
                           wy + DIR_DY[base_node->parent], root)
        ->is_open[DIR_OPP[base_node->parent]] = false;

    base_node->parent = DIR_W;
  } else if (ly == MAZE_RADIUS && lx != 0 && rand() & 1) {
    base_node->is_open[DIR_N] = true;
    _raw_get_cell(wx, wy + 1, root)->is_open[DIR_S] = true;

    base_node->is_open[base_node->parent] = false;
    _raw_get_cell(wx + DIR_DX[base_node->parent],
                           wy + DIR_DY[base_node->parent], root)
        ->is_open[DIR_OPP[base_node->parent]] = false;

    base_node->parent = DIR_N;
  }

  base_node->is_fixed = true;
  return base_node;
}

maze_t* infinite_maze_new(int seed) {
  int size = 2 * MAZE_RADIUS + 1;
  maze_t* m = calloc(1, sizeof(*m));
  m->seed = seed;
  m->data = calloc(size * size, sizeof(node_t));
  srand(m->seed);
  _node(m, 0, 0)->parent = -1;
  _explore(m, 0, 0);
  return m;
}

void infinite_maze_free(maze_t* m) {
  if (!m) return;

  if (m->outer_node) {
    infinite_maze_free(m->outer_node->outer_maze);
    return;
  }

  for (int y = -MAZE_RADIUS; y <= MAZE_RADIUS; y++) {
    for (int x = -MAZE_RADIUS; x <= MAZE_RADIUS; x++) {
      node_t* n = _node(m, x, y);
      if (n->inner_maze) {
        n->inner_maze->outer_node = NULL;
        infinite_maze_free(n->inner_maze);
      }
    }
  }

  free(m->data);
  free(m);
}

node_t* infinite_maze_get_cell(int wx, int wy, maze_t* root){
  _generate_final_cell(wx+1, wy, root);
  _generate_final_cell(wx-1, wy, root);
  _generate_final_cell(wx, wy+1, root);
  _generate_final_cell(wx, wy-1, root);

  return _generate_final_cell(wx, wy, root);
}

#endif /* INFINITE_MAZE_IMPLEMENTATION */
#endif /* INFINITE_MAZE_H */
