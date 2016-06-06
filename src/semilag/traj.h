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

#ifndef SRC_SEMILAG_TRAJ_H_
#define SRC_SEMILAG_TRAJ_H_

#include <vector>


namespace tbslas {

  //**********************************************************************
  // INTEFACE FUNCTIONS USING FUNCTORS
  //**********************************************************************
  template<typename real_t, typename FieldFunctor>
    void
    IntegrateRK2(FieldFunctor field_fn,
		 const std::vector<real_t>& xinit,
		 const real_t tinit,
		 const real_t dt,
		 std::vector<real_t>& xsol);

  template<typename real_t, typename Functor>
    void
    ComputeTrajRK2(Functor field_fn,
		   const std::vector<real_t>& xinit,
		   const real_t tinit,
		   const real_t tfinal,
		   const int num_rk_step,
		   std::vector<real_t>& xsol);

  template<typename real_t, typename FieldFunctor, typename ExtrapFunctor>
    void
    IntegrateRK2(FieldFunctor field_fn,
		 ExtrapFunctor extrap_fn,
		 const std::vector<real_t>& xinit,
		 const real_t tinit,
		 const real_t dt,
		 std::vector<real_t>& xsol);

  template<typename real_t, typename FieldFunctor, typename ExtrapFunctor>
    void
    ComputeTrajRK2(FieldFunctor field_fn,
		   ExtrapFunctor extrap_fn,
		   const std::vector<real_t>& xinit,
		   const real_t tinit,
		   const real_t tfinal,
		   const int num_rk_step,
		   std::vector<real_t>& xsol);

}  // namespace tbslas

#include "traj.inc"

#endif  // SRC_SEMILAG_TRAJ_H_
