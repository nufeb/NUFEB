#ifndef LMP_BOX_H
#define LMP_BOX_H

#include <array>

namespace LAMMPS_NS {
template <typename T, std::size_t N>
struct Box {
  Box() {
    for (int i = 0; i < N; i++) {
      lower[i] = T(0);
      upper[i] = T(0);
    }
  }
  Box(const Box &other) = default;
  Box(const std::array<T, N> &l, const std::array<T, N> &u) : lower(l), upper(u) {}
  template <typename ForwardIterator>
  Box(ForwardIterator l, ForwardIterator u) { 
    for(int i = 0; i < N; i++) {
      lower[i] = *l++;
      upper[i] = *u++;
    }
  }

  std::array<T, N> lower;
  std::array<T, N> upper;	
};

template <typename T, std::size_t N>
inline Box<T, N> intersect(const Box<T, N> &b0, const Box<T, N> &b1) {
  Box<T, N> result;
  for (int i = 0; i < N; i++) {
    result.lower[i] = b0.lower[i] > b1.lower[i] ? b0.lower[i] : b1.lower[i];
    result.upper[i] = b0.upper[i] < b1.upper[i] ? b0.upper[i] : b1.upper[i];
  }
  return result;
}

template <typename T, std::size_t N>
inline Box<T, N> extend(const Box<T, N> &b, T n = T(1))
{
  Box<T, N> result;
  for (int i = 0; i < N; i++) {
    result.lower[i] = b.lower[i] - n;
    result.upper[i] = b.upper[i] + n;
  }
  return result;
}

template <typename T, std::size_t N>
inline bool is_empty(const Box<T, N> &b)
{
  for (int i = 0; i < N; i++) {
    if (b.lower[i] >= b.upper[i])
      return true;
  }
  return false;
}

template <typename T, std::size_t N>
inline Box<T, N> translate(const Box<T, N> &b, const std::array<T, N> &d) {
  Box<T, N> result(b);
  for (int i = 0; i < N; i++) {
    result.lower[i] += d[i];
    result.upper[i] += d[i];
  }
  return result;
}

template <typename T, std::size_t N>
inline std::array<T, N> size(const Box<T, N> &b) {
  std::array<T, N> result;
  for (int i = 0; i < N; i++) {
    result[i] = b.upper[i] - b.lower[i];
  }
  return result;
}

template <typename T, std::size_t N>
inline bool is_inside(const Box<T, N> &b, const std::array<T, N> &point) {
  bool result = true;
  for (int i = 0; i < N; i++) {
    result &= point[i] >= b.lower[i] && point[i] < b.upper[i];
  }
  return result;
}
} // namespace LAMMPS_NS
#endif // LMP_BOX_H

