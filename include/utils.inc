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

#include <iostream>
#include <vector>
#include <cmath>

namespace tbslas {

template <typename real_t>
std::vector<real_t>
gen_reg_grid_points(size_t N) {
  size_t tN = N*N*N;
  real_t dx = 1.0/(N-1);
  // regular grid points positions
  std::vector<real_t> reg_grid_points_pos(tN*3);
  for (size_t i = 0; i < tN; i++) {
    // x cooridnate
    reg_grid_points_pos[i]      = (((real_t)((i / 1)    %N))*dx);
    // y coordinate
    reg_grid_points_pos[i+tN]   = (((real_t)((i / N)    %N))*dx);
    // z coordinate
    reg_grid_points_pos[i+2*tN] = (((real_t)((i /(N*N)) %N))*dx);
  }
  return reg_grid_points_pos;
}

template<typename real_t>
std::vector<real_t>
vorticity_field(std::vector<real_t>& points_pos, real_t omega) {
  size_t tN = points_pos.size()/3;
  std::vector<real_t> points_values(3*tN);
  for (size_t i = 0; i < tN; i++) {
    points_values[i]      = omega*(0.5         - points_pos[i+tN]) ;
    points_values[i+tN]   = omega*(points_pos[i] - 0.5           ) ;
    points_values[i+2*tN] = 0                                    ;
  }
  return points_values;
}

template<typename real_t>
std::vector<real_t>
gaussian_field(std::vector<real_t>& points_pos,
               real_t xc      ,
               real_t yc      ,
               real_t theta   ,
               real_t sigma_x ,
               real_t sigma_y) {
  real_t A           = 1.0;
  real_t cos_theta   = cos(theta);
  real_t cos_2theta  = cos(2*theta);
  real_t cos_theta_2 = cos_theta * cos_theta;
  real_t sin_theta   = sin(theta);
  real_t sin_2theta  = sin(2*theta);
  real_t sin_theta_2 = sin_theta * sin_theta;
  real_t sigma_x_2   = sigma_x*sigma_x;
  real_t sigma_y_2   = sigma_y*sigma_y;

  real_t a = cos_theta_2*0.5/sigma_x_2 + sin_theta_2*0.5/sigma_y_2;
  real_t b = -sin_2theta*0.5/sigma_x_2 + sin_2theta*0.5/sigma_y_2;
  real_t c = sin_theta_2*0.5/sigma_x_2 + cos_theta_2*0.5/sigma_y_2;

  size_t tN = points_pos.size()/3;
  std::vector<real_t> points_values(tN);
  for (size_t i = 0; i < tN; i++) {
    points_values[i]  = A*exp(-(a * (points_pos[i]-xc)    * (points_pos[i]-xc)    +
                                b * (points_pos[i]-xc)    * (points_pos[i+tN]-yc) +
                                c * (points_pos[i+tN]-yc) * (points_pos[i+tN]-yc)
                                )
                              );
  }
  return points_values;
}

}  // namespace tbslas
