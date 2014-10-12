// *************************************************************************
// Copyright (C) 2014 by Arash Bakhtiari

// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_SEMILAG_CUBIC_INTERP_POLICY_H_
#define SRC_SEMILAG_CUBIC_INTERP_POLICY_H_

#include <vector>

namespace slas {

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

  void
  interp(const std::vector<real_t>& grid_points,
         const int sdim,
         const std::vector< std::vector<real_t> >& grid_values,
         const int vdim,
         const std::vector<real_t>& time,
         const std::vector<real_t>& query_points,
         const real_t query_time,
         std::vector<real_t>& query_values) const;

 private:

  inline static real_t
  cubic_interp_1d(real_t p[4], real_t x);

  inline static real_t
  cubic_interp_2d(real_t p[4][4], real_t x, real_t y);

};  // class CubicInterpPolicy

}  // namespace slas

#include "cubic_interp_policy.inc"

#endif  // SRC_SEMILAG_CUBIC_INTERP_POLICY_H_
