#ifndef TGRID2D_H
#define TGRID2D_H

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

namespace utils {

template <typename T> class Grid2D {
  std::vector<T> grid;
  size_t Nx, Nz;

public:
  Grid2D() : Nx(0), Nz(0) {}

  Grid2D(size_t Nx, size_t Nz) { setGridSize(Nx, Nz); }

  Grid2D( const Grid2D &other_grid) {
    size_t Nx = other_grid.getNx();
    size_t Nz = other_grid.getNz();
    setGridSize(Nx, Nz);
  }

  void setGridSize(size_t Nx, size_t Nz) {
    this->Nx = Nx;
    this->Nz = Nz;
    grid.resize(Nx * Nz);
  }

  /** Accessor / Mutator */
  T &get(size_t ix, size_t iz) { return grid[ix * Nz + iz]; }
  T &get(const size_t &i) { return grid[i]; }

  /* Min / Max */
  T max() const { return *max_element(grid.begin(), grid.end()); }
  T min() const { return *min_element(grid.begin(), grid.end()); }

  /** Accessor */
  const T &get(size_t ix, size_t iz) const {
    assert(ix >= 0 && ix < Nx);
    assert(iz >= 0 && iz < Nz);
    return grid[ix * Nz + iz];
  }

  const T &get(const size_t &i) const { return grid[i]; }

  T getValue(size_t ix, size_t iz) { return grid[ix * Nz + iz]; }

  /** Return a reference to the grid values */
  std::vector<T> &getGrid() { return grid; }

  void clearGrid() { grid.clear(); }

  size_t getNx() const { return Nx; }
  void setNx(size_t nx) { Nx = nx; }
  size_t getNz() const { return Nz; }
  void setNz(size_t nz) { Nz = nz; }
  size_t getSize() const { return Nx * Nz; }
};

} // namespace utils

#endif /* TGRID2D_H_ */