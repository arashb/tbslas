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

namespace semilag {

  // unit length
  template <typename real_t, int dim>
  std::vector<real_t> gen_reg_grid_points(size_t N);

}  // namespace semilag

#include "utils.ipp"

#endif /* UTILS_HPP */
