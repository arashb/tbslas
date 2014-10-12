// *************************************************************************
// Copyright (C) 2014 by Arash Bakhtiari

// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_SEMILAG_UTILS_H_
#define SRC_SEMILAG_UTILS_H_

#include <cstdlib>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cmath>

namespace slas {

template <typename real_t, int dim>
void
get_reg_grid_points(const size_t N, real_t* out);

template<typename real_t, int sdim>
void
get_gaussian_field(const real_t* pnts_pos,
                   int num_points,
                   real_t* out);

template<typename real_t, int sdim>
void
get_vorticity_field(const real_t* points_pos,
                    int num_points,
                    real_t* out);

template <typename real_t, int dim>
std::vector<real_t>
generate_reg_grid_points(const size_t N);

template<typename real_t>
std::vector<real_t>
generate_vorticity_field(const std::vector<real_t>& pnts_pos,
                         const real_t time  = 0,
                         const real_t omega = 1);

template<typename real_t>
std::vector<real_t>
generate_gaussian_field(const std::vector<real_t>& pnts_pos,
                        const real_t xcenter = 0.6,
                        const real_t ycenter = 0.6,
                        const real_t theta   = 0.0,
                        const real_t sigma_x = 0.06,
                        const real_t sigma_y = 0.06);

template<typename real_t>
inline int
find_grid_index_1d(const std::vector<real_t>& grid,
                   const real_t query);

inline bool
is_little_endian() {
  int n = 1;
  return (*reinterpret_cast<char*>(&n) == 1);
}

}  // namespace slas

#include "utils.inc"

#endif  // SRC_SEMILAG_UTILS_H_
