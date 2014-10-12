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

#include <vector>

#include "semilag/vec_field.h"
//#include "semilag/cubic_interp_policy.h"

namespace slas {

template<typename real_t, typename FieldFunctor>
void
integrate_rk2(const FieldFunctor field_fn,
              const std::vector<real_t>& xinit,
              const real_t tinit,
              const real_t dt,
              std::vector<real_t>& xsol) {
  int num_points = xinit.size()/3;
  std::vector<real_t> xtmp(xinit.size());
  std::vector<real_t> vtmp(xinit.size());

  // FIXME: currently only time independent
  field_fn(xinit.data(), num_points, vtmp.data());
  for (int i = 0; i < xinit.size(); i++)
    xtmp[i] = xinit[i] + 0.5*dt*vtmp[i];

  // FIXME: currently only time independent
  field_fn(xtmp.data(), num_points, vtmp.data());
  for (int i = 0; i < xinit.size(); i++)
    xsol[i] = xinit[i] + dt*vtmp[i];
}

template<typename real_t, typename FieldFunctor>
void
traj_rk2(const FieldFunctor field_fn,
         const std::vector<real_t>& xinit,
         const real_t tinit,
         const real_t tfinal,
         const int num_rk_step,
         std::vector<real_t>& xsol) {
  // FIXME: currently only time independent
  real_t tau = (tfinal - tinit)/num_rk_step;

  for (int i = 0; i < xinit.size(); i++)
    xsol[i] = xinit[i];

  real_t current_time = tinit;
  // n-step RK2
  for (int taustep = 0; taustep < num_rk_step; taustep++) {
    // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
    integrate_rk2(field_fn, xsol, current_time, tau, xsol);
    current_time = current_time + tau;
  }
}

template<typename real_t, typename InterpPolicy, int sdim, int vdim>
void
integrate_rk2(const VecField<real_t, sdim, vdim>& field,
              const std::vector<real_t>& xinit,
              const real_t init_time,
              const real_t dt,
              const InterpPolicy& interp_policy,
              std::vector<real_t>& xsol) {
  std::vector<real_t> xtmp(xinit.size());
  std::vector<real_t> vtmp(xinit.size());

  field.interp(xinit, init_time, interp_policy, vtmp);
  for (int i = 0; i < xinit.size(); i++)
    xtmp[i] = xinit[i] + 0.5*dt*vtmp[i];

  field.interp(xtmp, init_time+0.5*dt, interp_policy, vtmp);
  for (int i = 0; i < xinit.size(); i++)
    xsol[i] = xinit[i] + dt*vtmp[i];
}

template<typename real_t, typename InterpPolicy, int sdim, int vdim>
void
traj_rk2(const VecField<real_t, sdim, vdim>& field,
         const std::vector<real_t>& xinit,
         const real_t tinit,
         const real_t tfinal,
         const int num_rk_step,
         const InterpPolicy& interp_policy,
         std::vector<real_t>& xsol) {
  real_t tau = (tfinal - tinit)/num_rk_step;

  for (int i = 0; i < xinit.size(); i++)
    xsol[i] = xinit[i];

  real_t current_time = tinit;
  // n-step RK2
  for (int taustep = 0; taustep < num_rk_step; taustep++) {
    // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
    integrate_rk2(field, xsol, current_time, tau, interp_policy, xsol);
    current_time = current_time + tau;
  }
}

}  // namespace slas
