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

#include <VecField.hpp>
#include <CubicInterpolant.hpp>

#include <vector>

namespace tbslas {

template<typename real_t, int sdim, int vdim>
void
integrate_rk2(VecField<real_t, sdim, vdim>& field,
              std::vector<real_t>& xinit,
              real_t dt,
              std::vector<real_t>& xsol) {
  // FIXME: currently time-independent vector field is implemented
  CubicInterpolant<real_t> intrplnt;
  std::vector<real_t> xtmp(xinit.size());
  std::vector<real_t> vtmp(xinit.size());

  field.interp(xinit, intrplnt, vtmp);
  for (int i = 0; i < xinit.size(); i++)
    xtmp[i] = xinit[i] + 0.5*dt*vtmp[i];

  field.interp(xtmp, intrplnt, vtmp);
  for (int i = 0; i < xinit.size(); i++)
    xsol[i] = xinit[i] + dt*vtmp[i];
}

template<typename real_t, int sdim, int vdim>
void
traj_rk2(VecField<real_t, sdim, vdim>& field,
         std::vector<real_t>& xinit,
         real_t tinit,
         real_t tfinal,
         const int n,
         std::vector<real_t>& xsol) {
  // FIXME: currently time-independent vector field is implemented
  real_t tau = (tfinal - tinit)/n;

  for (int i = 0; i < xinit.size(); i++)
      xsol[i] = xinit[i];

  // n-step RK2
  for (int taustep = 0; taustep < n; taustep++) {
    // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
    integrate_rk2(field, xsol, tau, xsol);
  }
}

}  // namespace tbslas
