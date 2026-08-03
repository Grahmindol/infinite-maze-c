#ifndef INFINITE_MAZE_H
#define INFINITE_MAZE_H

#ifdef __cplusplus
#include <cstdint>
#include <utility>
#include <vector>
#else
#include <stdbool.h>
#include <stdint.h>
#endif
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
 * Thread safety:
 *
 *  - Thread-safe as long as each thread operates on its own
 *    maze instance (`maze_p`).
 *
 *  - Multiple maze instances can safely be created with the
 *    same seed in different threads and will generate
 *    identical mazes deterministically.
 *
 *  - Accessing the same maze instance from multiple threads
 *    concurrently is NOT safe.
 *
 * ------------------------------------------------------------
 * License:
 *  Copyright (c) 2026 Grahmindol
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this software and associated documentation files (the "Software"), to
 *  deal in the Software without restriction, including without limitation the
 *  rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included in
 *  all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 *  IN THE SOFTWARE.
 *
 * Repository: https://github.com/Grahmindol/infinite-maze-c
 */

#ifndef MAZE_RADIUS
#define MAZE_RADIUS 1
#endif

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

#ifdef _WIN32
#define API __declspec(dllexport)
#elif defined(__EMSCRIPTEN__)
#define API EMSCRIPTEN_KEEPALIVE
#elif defined(__GNUC__)
#define API __attribute__((visibility("default")))
#else
#define API
#endif

#ifdef __cplusplus
extern "C" {
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
API bool infinite_maze_is_walkable(void* maze_p, int wx, int wy);

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
API uint8_t infinite_maze_get_cell(void* maze_p, int wx, int wy);

/**
 * @brief Retrieve hierarchical dead-end information for a world cell.
 *
 * Encodes walkability and hierarchical fixedness of the cell located
 * at world coordinates (`wx`, `wy`).
 *
 * A cell is considered *fixed* if it lies on the path connecting regions
 * at the corresponding hierarchy level. Consequently, the more bits that
 * are set, the more globally important the cell is: highly fixed cells are
 * shared by increasingly large-scale paths through the maze.
 *
 * A bit set to `1` indicates that the corresponding cell (or its ancestor
 * at that hierarchy level) is fixed.
 *
 * @param wx World-space X coordinate.
 * @param wy World-space Y coordinate.
 * @param maze_p Pointer to the root maze instance.
 *
 * @return uint8_t
 * - Encoded walkability and hierarchy information.
 * - Returns `0` if `maze_p` is NULL.
 *
 * @note Time complexity: O(log(wx² + wy²))
 *
 * @warning This function may trigger lazy generation of parent maze chunks.
 *
 * @see infinite_maze_is_walkable
 * @see infinite_maze_get_cell
 */
API uint8_t infinite_maze_get_fixedness(void* maze_p, int wx, int wy);

/**
 * @brief Compute the hierarchical center of a region.
 *
 * Computes the center coordinates of the region containing
 * the given world cell at a specified hierarchical level.
 * Higher levels correspond to larger structural groupings in the maze.
 *
 * @param wx World X coordinate of the cell.
 * @param wy World Y coordinate of the cell.
 * @param level Dead-end hierarchy level (0 = cell itself).
 * @param cx Output pointer receiving the center X coordinate.
 * @param cy Output pointer receiving the center Y coordinate.
 * @param maze_p Pointer to the root maze instance.
 *
 * @note Complexity: O(level)
 *
 * @note For level = 0, the returned center is (wx, wy).
 * @note For level = 1, the returned center is (wx | 1, wy &~1).
 */
API void infinite_maze_get_region_center(void* maze_p, int wx, int wy, int level, int* cx, int* cy);

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
 * @note Complexity: O(L) = O(exp(d))
 *
 * where L is the number of world cells along the path.
 *   and d is the distance of the farest point to the origine.
 *
 * @note The path is streamed in traversal order from start to target.
 *
 * @warning This function performs no internal synchronization.
 *          Concurrent calls must operate on distinct maze instances.
 */
API void infinite_maze_walk_from_to(void* maze_p, int fwx, int fwy, int twx, int twy,
                                    void (*walker)(int x, int y, void* user_data), void* user_data);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
class InfiniteMaze {
 public:
  // Create a new infinite procedurally generated maze.
  // @related infinite_maze_new
  explicit InfiniteMaze(uint64_t seed) { maze_ = infinite_maze_new(seed); }

  // Free an infinite maze and all generated chunks.
  // @related infinite_maze_free
  ~InfiniteMaze() { infinite_maze_free(maze_); }

  // Copy constructor
  InfiniteMaze(const InfiniteMaze& other) : InfiniteMaze(other.seed()) {}
  InfiniteMaze& operator=(const InfiniteMaze& other) {
    if (this != &other) {
      infinite_maze_free(maze_);
      maze_ = infinite_maze_new(other.seed());
    }
    return *this;
  }
  bool operator==(const InfiniteMaze& other) const noexcept { return seed() == other.seed(); }
  bool operator!=(const InfiniteMaze& other) const noexcept { return !(*this == other); }

  [[nodiscard]]
  uint64_t seed() const noexcept {
    return *static_cast<const uint64_t*>(maze_);
  }

  // Test whether a world-space cell is walkable.
  // @related infinite_maze_is_walkable
  [[nodiscard]]
  bool is_walkable(int wx, int wy) const {
    return infinite_maze_is_walkable(maze_, wx, wy);
  }

  // Retrieve hierarchical dead-end information for a world cell.
  // @related infinite_maze_get_cell
  [[nodiscard]]
  uint8_t get_cell(int wx, int wy) const {
    return infinite_maze_get_cell(maze_, wx, wy);
  }

  // @related infinite_maze_get_fixedness
  [[nodiscard]]
  uint8_t get_fixedness(int wx, int wy) const {
    return infinite_maze_get_fixedness(maze_, wx, wy);
  }

  // Compute the hierarchical center of a region.
  // @related infinite_maze_get_region_center
  [[nodiscard]]
  std::pair<int, int> get_region_center(int wx, int wy, int level) const {
    int cx, cy;
    infinite_maze_get_region_center(maze_, wx, wy, level, &cx, &cy);
    return {cx, cy};
  }

  // Stream the unique shortest path between two world coordinates.
  // @related infinite_maze_walk_from_to
  std::vector<std::pair<int, int>> get_path(int from_x, int from_y, int to_x, int to_y) {
    std::vector<std::pair<int, int>> path;
    infinite_maze_walk_from_to(
        maze_, from_x, from_y, to_x, to_y,
        [](int x, int y, void* user) {
          auto* path = static_cast<std::vector<std::pair<int, int>>*>(user);
          path->emplace_back(x, y);
        },
        &path);
    return path;
  }

  class _column_proxy {
   public:
    _column_proxy(const InfiniteMaze& maze, int x) : maze_(maze), x_(x) {}
    uint8_t operator[](int y) const { return maze_.get_cell(x_, y); }
   private:
    const InfiniteMaze& maze_;
    int x_;
  };
  _column_proxy operator[](int x) const { return _column_proxy(*this, x); }

 private:
  void* maze_ = nullptr;

};

#endif
/* =======================
   === IMPLEMENTATION ====
   ======================= */

#ifdef INFINITE_MAZE_IMPLEMENTATION

#include <stdlib.h>

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
  uint64_t seed;  // should be first !!!
  node_t* data;
  node_t* outer_node;
} maze_t;

static const int DIR_DX[DIR_COUNT] = {0, -1, 1, 0};
static const int DIR_DY[DIR_COUNT] = {1, 0, 0, -1};
static const direction_t DIR_OPP[DIR_COUNT] = {DIR_S, DIR_E, DIR_W, DIR_N};

/** Internal helpers */

static inline uint64_t rng(uint64_t* s) {
  uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

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

  rng(&m->seed);
  direction_t start = (direction_t)(m->seed % DIR_COUNT);
  for (int i = 0; i < (int)DIR_COUNT; i++) {
    direction_t d = (direction_t)((start + i) % DIR_COUNT);
    int nx = x + DIR_DX[d];
    int ny = y + DIR_DY[d];
    node_t* n = _node(m, nx, ny);
    if (n && n->outer_maze == NULL) {
      _carve(m, x, y, d);
      _explore(m, nx, ny);
    }
  }
}

static inline void _world_to_chunk(const int wx, const int wy, int* cx, int* cy, int* lx, int* ly) {
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
    if ((n->is_open[d] = m->outer_node->is_open[d])) _fix_parent_path(m, dx, dy);

    if (n->inner_maze) _update_chunk_aperture(n->inner_maze);
  }
}

/** Get cell in world coordinates with chunk generation */
static node_t* _get_raw_cell(int wx, int wy, maze_t* root) {
  if (abs(wx) <= MAZE_RADIUS && abs(wy) <= MAZE_RADIUS) return _node(root, wx, wy);

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
    outer_node->inner_maze =
        (maze_t*)infinite_maze_new(outer_maze->seed ^ (long long)(cx) * 0x9E3779B185EBCA87ULL ^
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

  uint64_t s =
      (base_node->outer_maze->seed ^ (lx) * 0x9E3779B185EBCA87ULL ^ (ly) * 0xC2B2AE3D27D4EB4FULL);

  if (lx == -MAZE_RADIUS && ly != 0 && rng(&s) & 0b11) {
    base_node->is_open[DIR_W] = 1;
    _get_raw_cell(wx - 1, wy, root)->is_open[DIR_E] = 1;

    base_node->is_open[base_node->parent_direction] = 0;
    _get_raw_cell(wx + DIR_DX[base_node->parent_direction],
                  wy + DIR_DY[base_node->parent_direction], root)
        ->is_open[DIR_OPP[base_node->parent_direction]] = 0;

    base_node->parent_direction = DIR_W;
  } else if (ly == MAZE_RADIUS && lx != 0 && rng(&s) & 0b11) {
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
  maze_t* m = (maze_t*)calloc(1, sizeof(*m));
  m->seed = seed;
  m->data = (node_t*)calloc(size * size, sizeof(node_t));
  _node(m, 0, 0)->parent_direction = DIR_COUNT;
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

API bool infinite_maze_is_walkable(void* maze_p, int wx, int wy) {
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

API void infinite_maze_get_region_center(void* maze_p, int wx, int wy, int level, int* cx,
                                         int* cy) {
  *cx = wx;
  *cy = wy;

  if (level <= 0) return;

  const int size = 2 * MAZE_RADIUS + 1;

  *cx >>= 1;
  *cy >>= 1;

  if (level <= 1) {
    *cx = (*cx << 1) + 1;
    *cy = (*cy << 1);
    return;
  }

  while (1) {
    node_t* cell = _get_refined_cell(*cx, *cy, (maze_t*)maze_p);
    if (cell->parent_direction < 0 || cell->parent_direction >= DIR_COUNT) break;

    *cx += DIR_DX[cell->parent_direction];
    *cy += DIR_DY[cell->parent_direction];
  }

  for (int i = 1; i < level; ++i) {
    int sx = (*cx > 0) - (*cx < 0);
    int sy = (*cy > 0) - (*cy < 0);

    *cx = (*cx + sx * MAZE_RADIUS) / size;
    *cy = (*cy + sy * MAZE_RADIUS) / size;
  }

  for (int i = 1; i < level; ++i) {
    *cx *= size;
    *cy *= size;
  }

  *cx = ((uint32_t)*cx << 1) + 1;
  *cy = ((uint32_t)*cy << 1);
}

typedef uint8_t (*node_predicate_t)(const node_t*);

static uint8_t _get_cell_impl(void* maze_p, int wx, int wy, node_predicate_t predicate) {
  if (!maze_p) return 0;
  maze_t* root = (maze_t*)maze_p;
  const int size = 2 * MAZE_RADIUS + 1;

  node_t* n = _get_raw_cell(wx >> 1, wy >> 1, root);
  uint8_t res = infinite_maze_is_walkable(maze_p, wx, wy);
  if (wy & 1 || ~wx & 1) {
    int xa, ya, xb, yb;
    if (wy & 1) {
      uint8_t neib = _get_cell_impl(maze_p, wx, wy - 1, predicate);
      for (int i = 2; i < 8; i++) {
        infinite_maze_get_region_center(maze_p, wx, wy + 1, i, &xa, &ya);
        infinite_maze_get_region_center(maze_p, wx, wy - 1, i, &xb, &yb);

        if (xa == xb && ya == yb) {
          res |= neib & (1 << i);
        } else {
          res &= ~(1 << i);
        }
      }
    }
    if (~wx & 1) {
      uint8_t neib = _get_cell_impl(maze_p, wx + 1, wy, predicate);
      for (int i = 2; i < 8; i++) {
        infinite_maze_get_region_center(maze_p, wx + 1, wy, i, &xa, &ya);
        infinite_maze_get_region_center(maze_p, wx - 1, wy, i, &xb, &yb);

        if (xa == xb && ya == yb) {
          res |= neib & (1 << i);
        } else {
          res &= ~(1 << i);
        }
      }
    }
    if (~wx & 1 && wy & 1) {
      uint8_t neib = _get_cell_impl(maze_p, wx + 1, wy - 1, predicate);
      for (int i = 2; i < 8; i++) {
        infinite_maze_get_region_center(maze_p, wx + 1, wy - 1, i, &xa, &ya);
        infinite_maze_get_region_center(maze_p, wx - 1, wy + 1, i, &xb, &yb);

        if (xa == xb && ya == yb) {
          infinite_maze_get_region_center(maze_p, wx - 1, wy - 1, i, &xa, &ya);
          infinite_maze_get_region_center(maze_p, wx + 1, wy + 1, i, &xb, &yb);
          if (xa == xb && ya == yb) {
            res |= neib & (1 << i);
          } else {
            res &= ~(1 << i);
          }
        } else {
          res &= ~(1 << i);
        }
      }
    }
    return res;
  }

  res |= predicate(n) << 1;
  infinite_maze_get_region_center(maze_p, wx, wy, 2, &wx, &wy);
  int x = wx >> 1, y = wy >> 1;

  for (int i = 2; i < 8; i++) {
    if (!root->outer_node) _get_raw_cell(x, y + (2 * MAZE_RADIUS + 1), root);
    root = root->outer_node->outer_maze;

    const int sx = (x > 0) - (x < 0), sy = (y > 0) - (y < 0);
    const int cx = (x + sx * MAZE_RADIUS) / size;
    const int cy = (y + sy * MAZE_RADIUS) / size;

    x = cx;
    y = cy;

    node_t* parent = _get_raw_cell(cx, cy, root);
    res |= predicate(parent) << i;
  }
  return res;
}

static unsigned char _is_dead_end(const node_t* n) {
  return (n->is_open[0] + n->is_open[1] + n->is_open[2] + n->is_open[3]) == 1;
}

API uint8_t infinite_maze_get_cell(void* maze_p, int wx, int wy) {
  return _get_cell_impl(maze_p, wx, wy, _is_dead_end);
}

static unsigned char _is_fixed(const node_t* n) { return !!n->is_fixed; }

API uint8_t infinite_maze_get_fixedness(void* maze_p, int wx, int wy) {
  return _get_cell_impl(maze_p, wx, wy, _is_fixed);
}

// Path finding module !

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
  p.data = (point_t*)malloc(sizeof(point_t) * initial_cap);
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
    p->data = (point_t*)realloc(p->data, sizeof(point_t) * p->capacity);
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
    if (cell->parent_direction < 0 || cell->parent_direction >= DIR_COUNT) break;

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

static path_t _hierarchical_path(int fwx, int fwy, int twx, int twy, maze_t* maze) {
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

  path_t high = _hierarchical_path(acrx, acry, bcrx, bcry, maze->outer_node->outer_maze);

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

API void infinite_maze_walk_from_to(void* maze_p, int fwx, int fwy, int twx, int twy,
                                    void (*walker)(int x, int y, void* user_data),
                                    void* user_data) {
  if (!maze_p || !walker) return;

  path_t path = _hierarchical_path(fwx >> 1, fwy >> 1, twx >> 1, twy >> 1, (maze_t*)maze_p);

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

#endif /* INFINITE_MAZE_IMPLEMENTATION */
#endif /* INFINITE_MAZE_H */