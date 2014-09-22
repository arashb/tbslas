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
#include <cassert>
#include "utils.h"

namespace tbslas {

template<typename real_t>
CubicInterpPolicy<real_t>::CubicInterpPolicy() {
}

template<typename real_t>
CubicInterpPolicy<real_t>::~CubicInterpPolicy() {
}

template<typename real_t>
void
CubicInterpPolicy<real_t>::interp(const std::vector<real_t>& grid_points,
                                  const int sdim,
                                  const std::vector<real_t>& grid_values,
                                  const int vdim,
                                  const std::vector<real_t>& query_points,
                                  std::vector<real_t>& query_values) const {
  int num_grid_points  = grid_points.size()/sdim;  // number of grid points
  int num_query_points = query_points.size()/sdim; // number of query points
  switch(sdim) {
    case 3: {
      // num_grid_points in one dimension
      float dnum_grid_points = pow(num_grid_points, 1./sdim);
      real_t spacing         = 1.0/(dnum_grid_points-1); // spacing

      // FIXME: construction of Interpolator here is not efficient
      for (int vdim_cnt = 0; vdim_cnt < vdim; vdim_cnt++) {
        const real_t* dgrid_values_array = &grid_values[vdim_cnt*num_grid_points];
        likely::TriCubicInterpolator tc_intrpltr(dgrid_values_array, spacing,
                                                 dnum_grid_points);
        for(int i = 0; i < num_query_points; i++)
          query_values[i+vdim_cnt*num_query_points] =
              tc_intrpltr(query_points[i],
                          query_points[i+num_query_points],
                          query_points[i+2*num_query_points]);
      }
      break;
    }
    case 2: {
      // TODO:
    }
    case 1:{
      // TODO:
    }
  }
}

template<typename real_t>
void
CubicInterpPolicy<real_t>::interp(const std::vector<real_t>& grid_points,
                                  const int sdim,
                                  const std::vector< std::vector<real_t> >& grid_values,
                                  const int vdim,
                                  const std::vector<real_t>& time,
                                  const std::vector<real_t>& query_points,
                                  const real_t query_time,
                                  std::vector<real_t>& query_values) const {
  size_t num_time_steps = time.size();
  size_t num_query_points = query_points.size()/sdim;
  // in case of 1d cubic interpolation, at least four elements are required
  assert (num_time_steps >= 4);

  int stencil[4] = {-1,0,1,2};
  int time_index = find_grid_index_1d(time, query_time);
  // change stencil at the begining and the end of the interval
  if(time_index == 0)
    stencil[0] = 0;  // repeat the frist point
  else if (time_index == num_time_steps-2)
    stencil[3] = 1;  // repeat the last ponit
  else if (time_index == num_time_steps-1) {
    stencil[0] = -1;
    stencil[1] = 0;
    stencil[2] = 0;
    stencil[3] = 0;
  }

// #ifndef NDEBUG
//   std::cout << stencil[0] << stencil[1] << stencil[2] << stencil[3] << std::endl;
// #endif

  // FIXME: optimize here by allocation of one big vector to avoid
  //        copy operation by push_back
  std::vector< std::vector<real_t> > qs_collec;
  std::vector<real_t> qs_values(num_query_points*vdim);
  for(int i = 0 ; i < 4; i++) {
    this->interp(grid_points,
                 sdim,
                 grid_values[time_index+stencil[i]],
                 vdim,
                 query_points,
                 qs_values);
    qs_collec.push_back(qs_values);
  }

  real_t p[4];
  // assume equidistant time spacing
  // query time should be in [0,1] interval -> normalize tq
  real_t tq = (query_time - time[time_index])/(time[time_index+1] - time[time_index]);
  int shift, point_value_index;
  for(int ivdim = 0; ivdim < vdim; ivdim++) {
    shift = ivdim*num_query_points;
    for( int iqp = 0; iqp < num_query_points; iqp++) {
      point_value_index = iqp + shift;
      // FIXME: inefficient memory access pattern
      p[0] = qs_collec[0][point_value_index];
      p[1] = qs_collec[1][point_value_index];
      p[2] = qs_collec[2][point_value_index];
      p[3] = qs_collec[3][point_value_index];
      query_values[point_value_index] = CubicInterpPolicy::cubic_interp_1d(p, tq);
    }
  }
}

template<typename real_t>
inline real_t
CubicInterpPolicy<real_t>::cubic_interp_1d(real_t p[4], real_t x) {
  return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3]
                                          + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

template<typename real_t>
inline real_t
CubicInterpPolicy<real_t>::cubic_interp_2d(real_t p[4][4],
                                           real_t x, real_t y) {
  double arr[4];
  arr[0] = cubic_interp_1d(p[0], y);
  arr[1] = cubic_interp_1d(p[1], y);
  arr[2] = cubic_interp_1d(p[2], y);
  arr[3] = cubic_interp_1d(p[3], y);
  return cubic_interp_1d(arr, x);
}


}  // namespace tbslas
