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

#ifndef INCLUDE_CUBIC_INTERPOLANT_H_
#define INCLUDE_CUBIC_INTERPOLANT_H_

#include <vector>

namespace tbslas {

template<typename real_t>
class CubicInterpolant {
 public:
  CubicInterpolant();
  virtual ~CubicInterpolant();

 public:
  void
  interp(std::vector<real_t>& grid_points,
         int sdim,
         std::vector<real_t>& grid_values,
         int vdim,
         std::vector<real_t>& query_points,
         std::vector<real_t>& query_values);
};  // class CubicInterpolant

}  // namespace tbslas

#include "cubic_interpolant.inc"

#endif  // INCLUDE_CUBIC_INTERPOLANT_H_
