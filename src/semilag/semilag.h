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

#ifndef SRC_SEMILAG_SEMILAG_H_
#define SRC_SEMILAG_SEMILAG_H_

#include <vector>

namespace tbslas {

  template<typename real_t,
    typename VFieldFunctor,
    typename CFieldFunctor>
    void
    SolveSemilagRK2(VFieldFunctor vel_evaluator,
		    CFieldFunctor con_evaluator,
		    const std::vector<real_t>& points_pos,
		    const int sdim,
		    const int timestep,
		    const real_t dt,
		    const int num_rk_step,
		    std::vector<real_t>& points_vals);

  template<typename real_t,
    typename VFieldFunctor,
    typename VExtrapFunctor,
    typename CFieldFunctor>
    void
    SolveSemilagRK2(VFieldFunctor vel_evaluator,
		    VExtrapFunctor extrap_evaluator,
		    CFieldFunctor con_evaluator,
		    const std::vector<real_t>& points_pos,
		    const int sdim,
		    const int timestep,
		    const real_t dt,
		    const int num_rk_step,
		    std::vector<real_t>& points_vals);

}  // namespace tbslas

#include "semilag.inc"

#endif  // SRC_SEMILAG_SEMILAG_H_
