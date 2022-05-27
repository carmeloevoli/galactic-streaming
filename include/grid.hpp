#ifndef TGRID2D_H
#define TGRID2D_H

#include <algorithm>
#include <functional>
#include <vector>
#include <iostream>
#include <cassert>

namespace utils {

template <typename T>
class Grid {
private:
  std::vector<T> m_data;
  size_t m_Nx, m_Nz;

 public:
  Grid() : m_Nx(0), m_Nz(0) {}

  Grid(size_t Nx, size_t Nz) { set_grid_size(Nx, Nz); }

  Grid(size_t Nx, size_t Nz, T value) {
    set_grid_size(Nx, Nz);
    std::fill(m_data.begin(), m_data.end(), value);
  }

  Grid(const Grid &t_grid) {
    size_t Nx = t_grid.get_Nx();
    size_t Nz = t_grid.get_Nz();
    set_grid_size(Nx, Nz);
    m_data = t_grid.m_data;
  }

  void set_grid_size(size_t Nx, size_t Nz) {
    m_Nx = Nx;
    m_Nz = Nz;
    m_data.resize(Nx * Nz);
  }

  /** Accessor / Mutator */
  T &get(size_t ix, size_t iz) { return m_data[ix * m_Nz + iz]; }
  T &get(const size_t &i) { return m_data[i]; }
  const T &get(size_t ix, size_t iz) const { return m_data[ix * m_Nz + iz]; }
  const T &get(const size_t &i) const { return m_data[i]; }

  /* Min / Max */
  T max() const { return *max_element(m_data.begin(), m_data.end()); }
  T min() const { return *min_element(m_data.begin(), m_data.end()); }

  T get_value(size_t ix, size_t iz) const { return m_data[ix * m_Nz + iz]; }
  void set_value(size_t ix, size_t iz, T value) { m_data[ix * m_Nz + iz] = value; }

  /** Return a reference to the grid values */
  std::vector<T> &get_data() const { return m_data; }

  void clear_grid() { m_data.clear(); }

  size_t get_Nx() const { return m_Nx; }
  void set_Nx(size_t nx) { m_Nx = nx; }
  size_t get_Nz() const { return m_Nz; }
  void set_Nz(size_t nz) { m_Nz = nz; }

  size_t get_size() const { return m_Nx * m_Nz; }

  /** Calculates the total size of the grid in bytes */
  size_t get_size_of() const { return sizeof(m_data) + (sizeof(m_data[0]) * m_data.size()); }

  /** Overload operators **/
  Grid<T> operator=(const Grid<T> &t_grid) {
	  const int t_Nx = t_grid.get_Nx();
	  const int t_Nz = t_grid.get_Nz();
	  set_grid_size(t_Nx, t_Nz);
	  m_data = t_grid.m_data;
	  return *this;
  }

  Grid<T> operator+=(const Grid<T> &grid) {
	  assert(this->get_size() == grid.get_size());
//	  auto data = grid.get_data();
	  transform(m_data.begin(), m_data.end(), grid.m_data.begin(), m_data.begin(), std::plus<T>());
//	  for(int i=0; i<this->get_size(); ++i) {
//		  this->m_data[i] += grid.m_data[i];
//	  }
	  return *this;
  }

  Grid<T> operator-=(Grid<T> &grid) {
    auto data = grid.get_data();
    transform(m_data.begin(), m_data.end(), data.begin(), m_data.begin(), std::minus<T>());
    return *this;
  }

  Grid<T> operator*=(T value) {
    transform(m_data.begin(), m_data.end(), m_data.begin(), [value](T &c) { return c * value; });
    return *this;
  }

  Grid<T> operator/=(T value) {
    transform(m_data.begin(), m_data.end(), m_data.begin(), [value](T &c) { return c / value; });
    return *this;
  }

  Grid<T> operator+(const Grid<T> &t_grid) const {
	  assert(this->get_size() == t_grid.get_size());
	  Grid<T> tmp = (*this);
	  tmp += t_grid;
	  return tmp;
  }

  Grid<T> operator*(T value) {
	  Grid<T> tmp = (*this);
	  tmp *= value;
	  return tmp;
  }

};

}  // namespace utils

#endif /* TGRID2D_H_ */
