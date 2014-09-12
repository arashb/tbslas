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

#ifndef INCLUDE_CUBIC_INTERP_POLICY_H_
#define INCLUDE_CUBIC_INTERP_POLICY_H_

#include <vector>

namespace tbslas {

template<typename real_t>
class CubicInterpPolicy {
 public:
  CubicInterpPolicy();
  virtual ~CubicInterpPolicy();

 public:
  void
  interp(const std::vector<real_t>& grid_points,
         const int sdim,
         const std::vector<real_t>& grid_values,
         const int vdim,
         const std::vector<real_t>& query_points,
         std::vector<real_t>& query_values) const;
};  // class CubicInterpPolicy

}  // namespace tbslas

#include "cubic_interp_policy.inc"

#endif  // INCLUDE_CUBIC_INTERP_POLICY_H_
