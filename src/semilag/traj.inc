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

#include <profile.hpp>

namespace tbslas {

  template<typename real_t, typename FieldFunctor>
  void
  IntegrateRK2(FieldFunctor field_fn,
	       const std::vector<real_t>& xinit,
	       const real_t tinit,
	       const real_t dt,
	       std::vector<real_t>& xsol) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    // pvfmm::Profile::Tic("IntegrateRK2", &sim_config->comm, false,5);
    int num_points = xinit.size()/3;
    std::vector<real_t> xtmp(xinit.size());
    std::vector<real_t> vtmp(xinit.size());

    // FIXME: currently only time independent
    // field_fn(xinit.data(), num_points, vtmp.data());
    // field_fn(xinit.data(), num_points, tinit+2*dt, vtmp.data());
    field_fn(xinit.data(), num_points, tinit, vtmp.data());
#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xtmp[i] = xinit[i] + 0.5*dt*vtmp[i];

    // FIXME: currently only time independent
    // field_fn(xtmp.data(), num_points, vtmp.data());
    field_fn(xtmp.data(), num_points, tinit+0.5*dt, vtmp.data());
#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xsol[i] = xinit[i] + dt*vtmp[i];

    pvfmm::Profile::Add_FLOP(2*3*xinit.size());
    // pvfmm::Profile::Toc();
  }

  template<typename real_t, typename FieldFunctor>
  void
  ComputeTrajRK2(FieldFunctor field_fn,
		 const std::vector<real_t>& xinit,
		 const real_t tinit,
		 const real_t tfinal,
		 const int num_rk_step,
		 std::vector<real_t>& xsol) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Tic("TrajRK2", &sim_config->comm, false,5);
    // FIXME: currently only time independent
    real_t tau = (tfinal - tinit)/num_rk_step;

#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xsol[i] = xinit[i];

    real_t current_time = tinit;
    // n-step RK2
    for (int taustep = 0; taustep < num_rk_step; taustep++) {
      // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
      IntegrateRK2(field_fn, xsol, current_time, tau, xsol);
      current_time = current_time + tau;
    }
    pvfmm::Profile::Toc();
  }

  template<typename real_t, typename FieldFunctor, typename ExtrapFunctor>
  void
  IntegrateRK2(FieldFunctor field_fn,
	       ExtrapFunctor extrap_fn,
	       const std::vector<real_t>& xinit,
	       const real_t tinit,
	       const real_t dt,
	       std::vector<real_t>& xsol) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    // pvfmm::Profile::Tic("IntegrateRK2", &sim_config->comm, false,5);
    int num_points = xinit.size()/3;
    std::vector<real_t> xtmp(xinit.size());
    std::vector<real_t> vtmp(xinit.size());

    field_fn(xinit.data(), num_points, vtmp.data());
#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xtmp[i] = xinit[i] + 0.5*dt*vtmp[i];

    // FIXME: currently only time independent
    extrap_fn(xtmp.data(), num_points, vtmp.data());
#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xsol[i] = xinit[i] + dt*vtmp[i];

    pvfmm::Profile::Add_FLOP(2*3*xinit.size());
    // pvfmm::Profile::Toc();
  }

  template<typename real_t,
	   typename FieldFunctor,
	   typename ExtrapFunctor>
  void
  ComputeTrajRK2(FieldFunctor field_fn,
		 ExtrapFunctor extrap_fn,
		 const std::vector<real_t>& xinit,
		 const real_t tinit,
		 const real_t tfinal,
		 const int num_rk_step,
		 std::vector<real_t>& xsol) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Tic("TrajRK2", &sim_config->comm, false,5);
    // FIXME: currently only time independent
    real_t tau = (tfinal - tinit)/num_rk_step;

#pragma omp parallel for
    for (int i = 0; i < xinit.size(); i++)
      xsol[i] = xinit[i];

    real_t current_time = tinit;
    // n-step RK2
    for (int taustep = 0; taustep < num_rk_step; taustep++) {
      // std::cout << xsol[0] << " " << xsol[1] << " " << xsol[2] << std::endl;
      IntegrateRK2(field_fn, extrap_fn, xsol, current_time, tau, xsol);
      current_time = current_time + tau;
    }
    pvfmm::Profile::Toc();
  }

}  // namespace tbslas
