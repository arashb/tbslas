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


/**
 * @file   semilag.cpp
 * @author Arash Bakhtiari <bakhtiar@in.tum.de>
 * @date   Wed Aug 27 19:53:22 2014
 * 
 * @brief  Contains the implementation of semilag library
 * 
 * 
 */

#include <semilag.hpp>
#include <vector>

namespace semilag {

  template <typename real_t>
  std::vector<real_t> reg_grid_points(size_t N) {
    std::vector<real_t> rgp_pos;    // regular grid points positions
    size_t total = N*N*N;
    for (size_t i = 0; i < total; i++) {
      rgp_pos.push_back(((real_t)((i / 1) %N))/N);  // x cooridnate
      rgp_pos.push_back(((real_t)((i / N) %N))/N);  // y coordinate
      rgp_pos.push_back(((real_t)((i /(N*N)) %N))/N);  // z coordinate
    }
    return rgp_pos;
  }

  // explicit instantiation
  template std::vector<double> reg_grid_points<double>(size_t);

}  // namespace semilag
