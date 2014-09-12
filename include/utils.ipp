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
  std::vector<real_t> rgp_pos(tN*3);    // regular grid points positions
  for (size_t i = 0; i < tN; i++) {
    rgp_pos[i]      = (((real_t)((i / 1)    %N))*dx);  // x cooridnate
    rgp_pos[i+tN]   = (((real_t)((i / N)    %N))*dx);  // y coordinate
    rgp_pos[i+2*tN] = (((real_t)((i /(N*N)) %N))*dx);  // z coordinate
  }
  return rgp_pos;
}

template<typename real_t>
std::vector<real_t>
vorticity_field(std::vector<real_t>& pnts_pos, real_t omega) {
  size_t tN = pnts_pos.size()/3;
  std::vector<real_t> pnts_vals(3*tN);
  for (size_t i = 0; i < tN; i++) {
    pnts_vals[i]      = omega*(0.5         - pnts_pos[i+tN]) ;
    pnts_vals[i+tN]   = omega*(pnts_pos[i] - 0.5           ) ;
    pnts_vals[i+2*tN] = 0                                    ;
  }
  return pnts_vals;
}

template<typename real_t>
std::vector<real_t>
gaussian_field(std::vector<real_t>& pnts_pos,
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

  size_t tN = pnts_pos.size()/3;
  std::vector<real_t> pnts_vals(tN);
  for (size_t i = 0; i < tN; i++) {
    pnts_vals[i]  = A*exp(-(a * (pnts_pos[i]-xc)    * (pnts_pos[i]-xc)    +
                            b * (pnts_pos[i]-xc)    * (pnts_pos[i+tN]-yc) +
                            c * (pnts_pos[i+tN]-yc) * (pnts_pos[i+tN]-yc)
                            )
                          );
  }
  return pnts_vals;
}

}  // namespace tbslas
