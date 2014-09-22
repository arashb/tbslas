// Copyright (C) 2014 by Arash Bakhtiari

// *************************************************************************
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef COMMON_H
#define COMMON_H

#include <vector>
#include <cassert>
#include <cmath>

inline bool
is_little_endian() {
  int n = 1;
  return (*reinterpret_cast<char*>(&n) == 1);
}

template<typename real_t>
inline int
find_grid_index_1d(const std::vector<real_t>& grid,
                   const real_t query) {
  int num_el = grid.size();
  assert(num_el >= 2);
  // assume grid is sorted
  real_t init    = grid[0];
  real_t spacing = grid[1] - grid[0];
  assert((query - init) >= 0);
  int index = (int)std::floor((query - init) / spacing);
  assert(index < num_el);
  return index;
}

#endif /* COMMON_H */
