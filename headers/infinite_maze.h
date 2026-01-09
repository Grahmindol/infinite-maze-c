#ifndef INFINITE_MAZE_H
#define INFINITE_MAZE_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/*
 * Infinite Maze — Single-header C library
 *
 * Procedurally generates an infinite perfect maze (tree: one unique path
 * between any two cells), lazily and deterministically from a seed.
 *
 * Features:
 *  - Infinite perfect maze in all directions
 *  - Deterministic generation (same seed → same maze)
 *
 * ------------------------------------------------------------
 * Usage:
 *
 *   // In EXACTLY ONE .c file:
 *   #define INFINITE_MAZE_IMPLEMENTATION
 *   #define MAZE_RADIUS 1 // the higher this is, the harder the maze is
 *   #include "infinite_maze.h"
 *
 *   // Everywhere else:
 *   #include "infinite_maze.h"
 *
 * ------------------------------------------------------------
 * Configuration:
 *
 *   MAZE_RADIUS : Radius (in cells) of a chunk.
 *                 Chunk size = (2 * MAZE_RADIUS + 1)²
 *                 Default = 1
 *
 * ------------------------------------------------------------
 * Thread safety: Not thread-safe (shared internal state)
 *
 * ------------------------------------------------------------
 * License:
 * Copyright (c) 2026 Grahmindol
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * Repository: https://github.com/Grahmindol/infinite-maze-c
 */

#ifndef MAZE_RADIUS
#define MAZE_RADIUS 1
#endif

#ifdef _WIN32
#define API __declspec(dllexport)
#else
#define API
#endif
/* =======================
   ======= API ===========
   ======================= */

/**
 * @brief Create a new infinite procedurally generated maze.
 *
 * Allocates and initializes the root chunk of an infinite maze.
 * The maze is generated deterministically from the given seed.
 * Additional chunks are generated lazily when accessed.
 *
 * @param seed Initial random seed controlling maze generation.
 * @return Pointer to the maze instance (opaque handle), or NULL on failure.
 *
 * @note Complexity: O(1)
 *
 */
API void* infinite_maze_new(int seed);

/**
 * @brief Free an infinite maze and all generated chunks.
 *
 * Recursively frees all chunks that have been generated so far,
 * including inner and outer mazes.
 *
 * @param m Maze instance returned by infinite_maze_new().
 *
 * @note Complexity: O(log(r))
 *
 * where r is distance to orgine to the farthest point generated
 *
 * @warning All pointers obtained from this maze become invalid.
 */
API void infinite_maze_free(void* m);

/**
 * @brief Test whether a world-space cell is walkable.
 *
 * Determines if the given world coordinate corresponds to an open
 * passage or a wall. Chunk generation may occur lazily if the
 * required area has not been generated yet.
 *
 * @param wx World X coordinate.
 * @param wy World Y coordinate.
 * @param maze_p Maze instance.
 *
 * @return true if the cell is walkable, false otherwise.
 *
 * @note Complexity:
 *
 * - average: O(1)
 *
 * - worst: O(ln(wx² + wy²)) if a new chunk is generated
 *
 */
API bool infinite_maze_is_walkable(int wx, int wy, void* maze_p);

/**
 * @brief Retrieve hierarchical dead-end information for a world cell.
 *
 * Encodes walkability and hierarchical dead-end levels of the cell located
 * at world coordinates (`wx`, `wy`).
 * Each bit in the returned byte represents a dead-end status at a given
 * hierarchy depth of the infinite maze.
 *
 * ### Bit layout of the returned byte
 *
 * - Bit 0 : Cell walkability (1 = walkable, 0 = wall)
 *
 * - Bit 1 : Dead-end status at local chunk level
 *
 * - Bits 2–7 : Dead-end status propagated through parent maze levels
 *
 * A bit set to `1` indicates that the corresponding cell (or its ancestor
 * at that hierarchy level) is a dead-end.
 *
 * @param wx World-space X coordinate.
 * @param wy World-space Y coordinate.
 * @param maze_p Pointer to the root maze instance.
 *
 * @return uint8_t
 * - Encoded walkability and dead-end hierarchy information.
 * - Returns `0` if `maze_p` is NULL.
 *
 * @note Time complexity: O(log(wx² + wy²))
 *
 * @warning This function may trigger lazy generation of parent maze chunks.
 *
 * @see infinite_maze_is_walkable
 */
API uint8_t infinite_maze_get_cell(int wx, int wy, void* maze_p);

/* =======================
   === IMPLEMENTATION ====
   ======================= */

#ifdef INFINITE_MAZE_IMPLEMENTATION

typedef enum {
  DIR_N = 0,
  DIR_W = 1,
  DIR_E = 2,
  DIR_S = 3,
  DIR_COUNT
} direction_t;

struct maze_t;

typedef struct node_t {
  bool is_open[DIR_COUNT];
  bool is_fixed;
  direction_t parent_direction;
  struct maze_t* inner_maze;
  struct maze_t* outer_maze;
} node_t;

typedef struct maze_t {
  int seed;
  node_t* data;
  node_t* outer_node;
} maze_t;

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
  a->is_open[d] = 1;
  b->is_open[DIR_OPP[d]] = 1;
  b->parent_direction = DIR_OPP[d];
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

static inline void _world_to_chunk(const int wx, const int wy, int* cx, int* cy,
                                   int* lx, int* ly) {
  const int size = 2 * MAZE_RADIUS + 1;
  const int sx = (wx > 0) - (wx < 0);
  const int sy = (wy > 0) - (wy < 0);

  *cx = (wx + sx * MAZE_RADIUS) / size;
  *cy = (wy + sy * MAZE_RADIUS) / size;
  *lx = wx - (*cx) * size;
  *ly = wy - (*cy) * size;
}

static void _fix_parent_path(maze_t* m, int lx, int ly) {
  node_t* n = _node(m, lx, ly);
  n->is_fixed = 1;
  if (n->parent_direction < 0 || n->parent_direction >= DIR_COUNT) return;

  lx += DIR_DX[n->parent_direction];
  ly += DIR_DY[n->parent_direction];
  _fix_parent_path(m, lx, ly);
}

static void _update_chunk_aperture(maze_t* m) {
  if (!m || !m->outer_node) return;

  for (int d = 0; d < DIR_COUNT; d++) {
    int dx = MAZE_RADIUS * DIR_DX[d];
    int dy = MAZE_RADIUS * DIR_DY[d];
    node_t* n = _node(m, dx, dy);
    if ((n->is_open[d] = m->outer_node->is_open[d]))
      _fix_parent_path(m, dx, dy);

    if (n->inner_maze) _update_chunk_aperture(n->inner_maze);
  }
}

/** Get cell in world coordinates with chunk generation */
static node_t* _get_raw_cell(int wx, int wy, maze_t* root) {
  if (abs(wx) <= MAZE_RADIUS && abs(wy) <= MAZE_RADIUS)
    return _node(root, wx, wy);

  int cx, cy, lx, ly;
  _world_to_chunk(wx, wy, &cx, &cy, &lx, &ly);

  if (!root->outer_node) {
    maze_t* outer_maze = (maze_t*)infinite_maze_new(root->seed + 1);
    root->outer_node = _node(outer_maze, 0, 0);
    root->outer_node->inner_maze = root;
    _update_chunk_aperture(root);
  }

  maze_t* outer_maze = root->outer_node->outer_maze;
  node_t* outer_node = _get_raw_cell(cx, cy, outer_maze);

  if (!outer_node->inner_maze) {
    outer_node->inner_maze = (maze_t*)infinite_maze_new(
        outer_maze->seed ^ (long long)(cx) * 0x9E3779B185EBCA87ULL ^
        (long long)(cy) * 0xC2B2AE3D27D4EB4FULL);
    outer_node->inner_maze->outer_node = outer_node;
    _update_chunk_aperture(outer_node->inner_maze);
  }

  return _node(outer_node->inner_maze, lx, ly);
}

static inline node_t* _get_refined_cell(int wx, int wy, maze_t* root) {
  node_t* base_node = _get_raw_cell(wx, wy, root);
  if (base_node->is_fixed) return base_node;

  int cx, cy, lx, ly;
  _world_to_chunk(wx, wy, &cx, &cy, &lx, &ly);

  if ((abs(lx) != MAZE_RADIUS && abs(ly) != MAZE_RADIUS)) return base_node;

  srand(base_node->outer_maze->seed ^ (lx) * 0x9E3779B185EBCA87ULL ^
        (ly) * 0xC2B2AE3D27D4EB4FULL);

  if (lx == -MAZE_RADIUS && ly != 0 && rand() & 0b10) {
    base_node->is_open[DIR_W] = 1;
    _get_raw_cell(wx - 1, wy, root)->is_open[DIR_E] = 1;

    base_node->is_open[base_node->parent_direction] = 0;
    _get_raw_cell(wx + DIR_DX[base_node->parent_direction],
                  wy + DIR_DY[base_node->parent_direction], root)
        ->is_open[DIR_OPP[base_node->parent_direction]] = 0;

    base_node->parent_direction = DIR_W;
  } else if (ly == MAZE_RADIUS && lx != 0 && rand() & 1) {
    base_node->is_open[DIR_N] = 1;
    _get_raw_cell(wx, wy + 1, root)->is_open[DIR_S] = 1;

    base_node->is_open[base_node->parent_direction] = 0;
    _get_raw_cell(wx + DIR_DX[base_node->parent_direction],
                  wy + DIR_DY[base_node->parent_direction], root)
        ->is_open[DIR_OPP[base_node->parent_direction]] = 0;

    base_node->parent_direction = DIR_N;
  }

  base_node->is_fixed = 1;
  return base_node;
}

API void* infinite_maze_new(int seed) {
  int size = 2 * MAZE_RADIUS + 1;
  maze_t* m = calloc(1, sizeof(*m));
  m->seed = seed;
  m->data = calloc(size * size, sizeof(node_t));
  srand(m->seed);
  _node(m, 0, 0)->parent_direction = -1;
  _explore(m, 0, 0);
  return (void*)m;
}

API void infinite_maze_free(void* maze_p) {
  if (!maze_p) return;
  maze_t* m = (maze_t*)maze_p;

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

API bool infinite_maze_is_walkable(int wx, int wy, void* maze_p) {
  if (!maze_p) return 0;
  if ((wx ^ wy) & 1) return wx & 1;

  maze_t* root = (maze_t*)maze_p;

  int x = wx >> 1;
  int y = wy >> 1;
  _get_refined_cell(x + 1, y, root);
  _get_refined_cell(x - 1, y, root);
  _get_refined_cell(x, y + 1, root);
  _get_refined_cell(x, y - 1, root);
  return _get_refined_cell(x, y, root)->is_open[!(wx & 1)];
}

static inline bool _is_dead_end(node_t* n) {
  return (n->is_open[0] + n->is_open[1] + n->is_open[2] + n->is_open[3]) == 1;
}

API uint8_t infinite_maze_get_cell(int wx, int wy, void* maze_p) {
  if (!maze_p) return 0;
  maze_t* root = maze_p;

  const int size = 2 * MAZE_RADIUS + 1;
  int x = wx >> 1, y = wy >> 1;

  uint8_t res = infinite_maze_is_walkable(wx, wy, maze_p);
  node_t* n = _get_raw_cell(x, y, root);

  if ((wx ^ wy) & 1 && wx & 1) {
    res |= _is_dead_end(n) << 1;
  }

  for (int i = 2; i < 8; i++) {
    if (!root->outer_node) _get_raw_cell(x, y + (2 * MAZE_RADIUS + 1), root);
    root = root->outer_node->outer_maze;

    int sx = (x > 0) - (x < 0), sy = (y > 0) - (y < 0);
    x = (x + sx * MAZE_RADIUS) / size;
    y = (y + sy * MAZE_RADIUS) / size;
    node_t* parent = _get_raw_cell(x, y, root);
    res |= _is_dead_end(parent) << i;
  }
  return res;
}

#endif /* INFINITE_MAZE_IMPLEMENTATION */
#endif /* INFINITE_MAZE_H */