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

#ifndef SEMILAG_HPP
#define SEMILAG_HPP

#include <cstdlib>
#include <vector>

namespace semilag {

  // unit length
  template <typename real_t>
  std::vector<real_t> gen_reg_grid_points(size_t N);

  // template<typename real_t>
  // std::vector<real_t> traj_rk2(std::vector<real_t> ipos,
  //                              fvel,
  //                              real_t ti,
  //                              real_t tf,
  //                              int n);

} // namespace: semilag

#include "semilag.ipp"

#endif /* SEMILAG_HPP */
