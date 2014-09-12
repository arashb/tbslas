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

#include <TriCubicInterpolator.h>
#include <vector>
#include <cmath>

namespace tbslas {

template<typename real_t>
CubicInterpolant<real_t>::CubicInterpolant() {
}

template<typename real_t>
CubicInterpolant<real_t>::~CubicInterpolant() {
}

template<typename real_t>
void
CubicInterpolant<real_t>::interp(std::vector<real_t>& grid_points,
                                 int sdim,
                                 std::vector<real_t>& grid_values,
                                 int vdim,
                                 std::vector<real_t>& query_points,
                                 std::vector<real_t>& query_values) {
  int num_grid_points    = grid_points.size()/sdim;  // number of grid points
  // num_grid_points in one dimension
  float dnum_grid_points = pow(num_grid_points, 1./sdim);
  real_t spacing           = 1.0/(dnum_grid_points-1); // spacing

  // FIXME: move the construction of Interpolator
  // to the constructor of this class.
  // Here only initialize it with values.
  int num_query_points = query_points.size()/sdim; // number of query points
  for (int vdim_cnt = 0; vdim_cnt < vdim; vdim_cnt++) {
    real_t* dgrid_values_array = &grid_values[vdim_cnt*num_grid_points];
    likely::TriCubicInterpolator tc_intrpltr(dgrid_values_array, spacing,
                                             dnum_grid_points);
    for(int i = 0; i < num_query_points; i++)
      query_values[i+vdim_cnt*num_query_points] =
          tc_intrpltr(query_points[i],
                      query_points[i+num_query_points],
                      query_points[i+2*num_query_points]);
  }
}

}  // namespace tbslas
