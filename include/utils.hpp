// Copyright (C) 2014 by Arash Bakhtiari

// *************************************************************************
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdlib>
#include <vector>
#include <sstream>
#include <iomanip>

namespace tbslas {

// unit length
template <typename real_t, int dim>
std::vector<real_t>
gen_reg_grid_points(size_t N);

template<typename real_t>
std::vector<real_t>
vorticity_field(std::vector<real_t>& pnts_pos,
                real_t omega = 1);

template<typename real_t>
std::vector<real_t>
gaussian_field(std::vector<real_t>& pnts_pos,
               real_t xcenter = 0.6,
               real_t ycenter = 0.6,
               real_t theta   = 0.0,
               real_t sigma_x = 0.06,
               real_t sigma_y = 0.06);

}  // namespace tbslas

#include "utils.ipp"

#endif /* UTILS_HPP */
