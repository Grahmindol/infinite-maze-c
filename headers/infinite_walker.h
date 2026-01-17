#ifndef INFINITE_WALKER_H
#define INFINITE_WALKER_H

/*
 * Infinite Walker — Single-header C pathfinding extension
 *
 * Hierarchicalpathfinding for the Infinite Maze library.
 *
 * Computes deterministic shortest paths between two world-space
 * coordinates inside an infinite perfect maze, using a multi-level
 * strategy
 *
 * Path reconstruction is streamed via a user callback.
 *
 * ------------------------------------------------------------
 * Usage:
 *
 *   // In EXACTLY ONE .c file:
 *   #define INFINITE_MAZE_IMPLEMENTATION
 *   #include "infinite_maze.h"
 *   #include "infinite_walker.h" // after infinite_maze !!!!
 *
 *   // Everywhere else:
 *   #include "infinite_maze.h"
 *   #include "infinite_walker.h" // after infinite_maze !!!!
 *
 * ------------------------------------------------------------
 * Thread safety:
 *
 *  - Thread-safe as long as each thread operates on its own
 *    maze instance (`maze_p`).
 *
 * ------------------------------------------------------------
 * License:
 *   Copyright (c) 2026 Grahmindol
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 *   furnished to do so, subject to the following conditions:
 *
 *   The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 *
 * Repository: https://github.com/Grahmindol/infinite-maze-c
 */

#ifndef INFINITE_MAZE_H
#error \
    "Invalid include order: <infinite_path.h> depends on <infinite_maze.h>. Please include <infinite_maze.h> first."
#endif

/**
 * @brief Stream the unique shortest path between two world coordinates.
 *
 * Computes and streams the deterministic shortest path between two
 * world-space coordinates inside the given infinite maze.
 *
 * The maze being perfect (acyclic), the path is guaranteed to be unique.
 *
 * @param fwx Starting world X coordinate (cell or corridor).
 * @param fwy Starting world Y coordinate (cell or corridor).
 * @param twx Target world X coordinate (cell or corridor).
 * @param twy Target world Y coordinate (cell or corridor).
 * @param maze_p Pointer to a maze instance created by infinite_maze_new().
 * @param walker Callback invoked for each world coordinate along the path.
 * @param user_data User-defined pointer passed unchanged to the callback.
 *
 * @note Complexity: O(L)
 *
 * where L is the number of world cells along the path.
 *
 * @note The path is streamed in traversal order from start to target.
 *
 * @warning This function performs no internal synchronization.
 *          Concurrent calls must operate on distinct maze instances.
 */
API void infinite_path_walk(int fwx, int fwy, int twx, int twy, void* maze_p,
                            void (*walker)(int x, int y, void* user_data),
                            void* user_data);

/* =======================
   === IMPLEMENTATION ====
   ======================= */

#ifdef INFINITE_MAZE_IMPLEMENTATION

#include <stdbool.h>
#include <stdlib.h>

typedef struct point_t {
  int x, y;
} point_t;

typedef struct {
  point_t* data;
  int size;
  int capacity;
} path_t;

static inline path_t _path_create(int initial_cap) {
  path_t p;
  p.data = malloc(sizeof(point_t) * initial_cap);
  p.size = 0;
  p.capacity = initial_cap;
  return p;
}

static inline void _path_free(path_t* p) {
  free(p->data);
  p->data = NULL;
  p->size = p->capacity = 0;
}

static inline void _path_push(path_t* p, point_t pt) {
  if (p->size == p->capacity) {
    p->capacity *= 2;
    p->data = realloc(p->data, sizeof(point_t) * p->capacity);
  }
  p->data[p->size++] = pt;
}

static inline void _path_pop(path_t* p) {
  if (p->size == 0) return;
  p->size--;
}

static inline void _path_reverse(path_t* p) {  // very costly to avoid
  for (int i = 0, j = p->size - 1; i < j; ++i, --j) {
    point_t tmp = p->data[i];
    p->data[i] = p->data[j];
    p->data[j] = tmp;
  }
}

static inline void _path_concat(path_t* dst, const path_t* src) {
  for (int i = 0; i < src->size; ++i) _path_push(dst, src->data[i]);
}

static inline path_t _get_path_to_root(maze_t* m, int wx, int wy) {
  int cap = (2 * MAZE_RADIUS + 1);
  cap = cap * cap;

  path_t path = _path_create(cap);

  int cx = wx, cy = wy;
  while (1) {
    _path_push(&path, (point_t){cx, cy});

    node_t* cell = _get_raw_cell(cx, cy, m);
    if (cell->parent_direction < 0 || cell->parent_direction >= DIR_COUNT)
      break;

    cx += DIR_DX[cell->parent_direction];
    cy += DIR_DY[cell->parent_direction];
  }
  return path;
}

static inline path_t _concat_and_merge(const path_t* AB, const path_t* BC) {
  path_t AC = _path_create(AB->size + BC->size);

  int i = AB->size - 1;
  int j = 0;

  while (i >= 0 && j < BC->size && AB->data[i].x == BC->data[j].x &&
         AB->data[i].y == BC->data[j].y) {
    i--;
    j++;
  }

  for (int k = 0; k <= i; ++k) {
    _path_push(&AC, AB->data[k]);
  }

  if (j > 0) {
    _path_push(&AC, BC->data[j - 1]);
  }

  for (int k = j; k < BC->size; ++k) {
    _path_push(&AC, BC->data[k]);
  }
  return AC;
}

static inline path_t _get_local_path(maze_t* m, point_t a, point_t b) {
  path_t pa = _get_path_to_root(m, a.x, a.y);
  path_t pb = _get_path_to_root(m, b.x, b.y);

  _path_reverse(&pb);

  path_t res = _concat_and_merge(&pa, &pb);

  _path_free(&pa);
  _path_free(&pb);

  return res;
}

static path_t _hierarchical_path(int fwx, int fwy, int twx, int twy,
                                 maze_t* maze) {
  path_t pa = _get_path_to_root(maze, fwx, fwy);
  path_t pb = _get_path_to_root(maze, twx, twy);

  point_t ra = pa.data[pa.size - 1];
  point_t rb = pb.data[pb.size - 1];

  // même racine → chemin local direct
  if (ra.x == rb.x && ra.y == rb.y) {
    _path_reverse(&pb);
    path_t res = _concat_and_merge(&pa, &pb);
    _path_free(&pa);
    _path_free(&pb);
    return res;
  }

  // passage niveau supérieur
  int acrx, acry, alx, aly;
  int bcrx, bcry, blx, bly;

  _world_to_chunk(ra.x, ra.y, &acrx, &acry, &alx, &aly);
  _world_to_chunk(rb.x, rb.y, &bcrx, &bcry, &blx, &bly);

  path_t high =
      _hierarchical_path(acrx, acry, bcrx, bcry, maze->outer_node->outer_maze);

  path_t middle = _path_create((2 * MAZE_RADIUS + 1) * (2 * MAZE_RADIUS + 1));

  int dim = 2 * MAZE_RADIUS + 1;

  for (int i = 0; i < high.size; ++i) {
    point_t h = high.data[i];

    point_t in = {h.x * dim, h.y * dim};
    point_t out = {h.x * dim, h.y * dim};

    if (i > 0) {
      point_t hp = high.data[i - 1];
      int dx = (hp.x > h.x) - (hp.x < h.x);
      int dy = (hp.y > h.y) - (hp.y < h.y);

      in.x += dx * MAZE_RADIUS;
      in.y += dy * MAZE_RADIUS;
    } else {
      in.x = fwx;
      in.y = fwy;
    }

    if (i + 1 < high.size) {
      point_t hn = high.data[i + 1];
      int dx = (hn.x > h.x) - (hn.x < h.x);
      int dy = (hn.y > h.y) - (hn.y < h.y);

      out.x += dx * MAZE_RADIUS;
      out.y += dy * MAZE_RADIUS;
    } else {
      out.x = twx;
      out.y = twy;
    }

    path_t seg = _get_local_path(maze, in, out);
    _path_concat(&middle, &seg);
    _path_free(&seg);
  }

  _path_free(&pa);
  _path_free(&pb);
  _path_free(&high);
  return middle;
}

API void infinite_path_walk(int fwx, int fwy, int twx, int twy, void* maze_p,
                            void (*walker)(int x, int y, void* user_data),
                            void* user_data) {
  if (!maze_p || !walker) return;

  path_t path = _hierarchical_path(fwx >> 1, fwy >> 1, twx >> 1, twy >> 1,
                                   (maze_t*)maze_p);

  if (path.size == 0) {
    _path_free(&path);
    return;
  }

  int i = 0;

  if (!((fwx ^ fwy) & 1)) {
    walker(fwx, fwy, user_data);

    int x1 = (path.data[0].x << 1) | 1;
    int y1 = (path.data[0].y << 1);

    if (path.size > 1) {
      int x2 = (path.data[1].x << 1) | 1;
      int y2 = (path.data[1].y << 1);
      if ((x1 + x2) / 2 == fwx && (y1 + y2) / 2 == fwy) i = 1;
    }
  }

  // parcours linéaire
  for (; i < path.size; ++i) {
    int x1 = (path.data[i].x << 1) | 1;
    int y1 = (path.data[i].y << 1);
    walker(x1, y1, user_data);

    if (i + 1 < path.size) {
      int x2 = (path.data[i + 1].x << 1) | 1;
      int y2 = (path.data[i + 1].y << 1);

      int mid_x = (x1 + x2) / 2;
      int mid_y = (y1 + y2) / 2;

      walker(mid_x, mid_y, user_data);

      if (mid_x == twx && mid_y == twy) break;
    } else if (x1 != twx || y1 != twy) {
      walker(twx, twy, user_data);
    }
  }

  _path_free(&path);
}

#endif
#endif