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
#include <cstdlib>
#include <iostream>

namespace tbslas {

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
            CField& conc_field) {
  real_t tinit  = timestep*dt;
  real_t tfinal = tinit - dt;  // backward in time
  const std::vector<real_t> points_pos = conc_field.get_points();
  const size_t num_points = conc_field.get_num_points();
  std::vector<real_t> departure_pos(points_pos.size());

  std::cout << "-> semi-Lagrangian" << std::endl;
  // computing the departure points
  std::cout << "--> computing the backward trajectory" << std::endl;
  traj_rk2(vel_field, points_pos, tinit, tfinal, num_rk_step, interp_policy,
           departure_pos);
  std::cout << "--> interplating at the departure points" << std::endl;
  // interpolating the concectration values at the departure points
  std::vector<real_t> departure_points_values(num_points);
  conc_field.interp(departure_pos,
                    interp_policy,
                    departure_points_values,
                    timestep-1);             // using previeous time step values
                                             // to interpolate
  conc_field.push_back_values(departure_points_values, dt*timestep);
}

}
