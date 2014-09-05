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
std::vector<real_t>
CubicInterpolant<real_t>::interp(std::vector<real_t>& grid_pnts, int sdim,
                                 std::vector<real_t>& grid_vals, int vdim,
                                 std::vector<real_t>& query_pnts) {
  // FIXME: at the moment only support 3d
  // if (sdim != 3) return;
  // FIXME: at moment only suppert one dimensional values
  // if (vdim != 1) return;

  int pnt_cnt = grid_pnts.size()/sdim;
  // FIXME: is correct only for the regular grid
  float dpnt_cnt = pow(pnt_cnt, 1./sdim);
  real_t spacing = 1.0/(dpnt_cnt-1);

  // FIXME: move this the construction of Interpolator
  // to the constructor of this class.
  // Here only initialize it with values.
  likely::TriCubicInterpolator tc_intrpltr(grid_vals, spacing, dpnt_cnt);

  int qpnt_cnt = query_pnts.size()/sdim;
  std::vector<real_t> query_vals;
  for(int i=0; i < qpnt_cnt; i += sdim)
    query_vals.push_back(tc_intrpltr(query_pnts[i],
                                     query_pnts[i+1],
                                     query_pnts[i+2]));
  return query_vals;
}

}  // namespace tbslas
