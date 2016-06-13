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

#include <mpi.h>
#include <vector>
#include <cstdlib>
#include <iostream>

#include <profile.hpp>

#include "utils/common.h"
#include "semilag/traj.h"

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
		  std::vector<real_t>& points_vals) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    // pvfmm::Profile::Tic("SemilagRK2", &sim_config->comm, false,5);
    const size_t num_points = points_pos.size()/sdim;
    real_t tinit  = timestep*dt;
    real_t tfinal = tinit - dt;  // backward in time

    std::vector<real_t> departure_points_pos(points_pos.size());

    // computing the departure points
    ComputeTrajRK2(vel_evaluator,
		   points_pos,
		   tinit,
		   tfinal,
		   num_rk_step,
		   departure_points_pos);

    con_evaluator(departure_points_pos.data(), num_points,  points_vals.data());
    // pvfmm::Profile::Toc();
  }

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
		  std::vector<real_t>& points_vals) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    // pvfmm::Profile::Tic("SemilagRK2", &sim_config->comm, false,5);
    const size_t num_points = points_pos.size()/sdim;
    real_t tinit  = timestep*dt;
    real_t tfinal = tinit - dt;  // backward in time

    std::vector<real_t> departure_points_pos(points_pos.size());

    // computing the departure points
    ComputeTrajRK2(vel_evaluator,
		   extrap_evaluator,
		   points_pos,
		   tinit,
		   tfinal,
		   num_rk_step,
		   departure_points_pos);

    con_evaluator(departure_points_pos.data(), num_points,  points_vals.data());
    // pvfmm::Profile::Toc();
  }

}  // namespace tbslas
