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

#include "semilag/common.h"

namespace tbslas {

template<typename real_t, typename FieldFunctor>
void
semilag_rk2(const FieldFunctor vel_evaluator,
            const FieldFunctor con_evaluator,
            const std::vector<real_t>& points_pos,
            const int timestep,
            const real_t dt,
            const int num_rk_step,
            std::vector<real_t>& points_vals);

template<typename real_t,
         typename InterpPolicy,
         typename VField,
         typename CField>
void
semilag_rk2(const VField& vel_field,
            const InterpPolicy& interp_policy,
            const int timestep,
            const real_t dt,
            const int num_rk_step,
            CField& conc_field);



}  // namespace tbslas

#include "semilag.inc"

#endif  // SRC_SEMILAG_SEMILAG_H_
