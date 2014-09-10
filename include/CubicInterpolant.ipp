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
CubicInterpolant<real_t>::interp(std::vector<real_t>& grd_pnts,
                                 int sdim,
                                 std::vector<real_t>& grd_vls,
                                 int vdim,
                                 std::vector<real_t>& qry_pnts) {
  std::vector<real_t> qry_vls;          // interpolated values of query points
  // FIXME: at the moment only support 3d
  if (sdim != 3) return qry_vls;
  int ngrd_pnts = grd_pnts.size()/sdim;  // number of grid points
  float dngrd_pnts = pow(ngrd_pnts, 1./sdim);
  real_t spcng = 1.0/(dngrd_pnts-1);  // spacing

  // FIXME: move the construction of Interpolator
  // to the constructor of this class.
  // Here only initialize it with values.
  for (int vdim_cnt = 0; vdim_cnt < vdim; vdim_cnt++) {
    // see the header file to understand the memory layout of the
    // vector field points and values of correspoding points.
    real_t* dgrd_vls_array = &grd_vls[vdim_cnt*ngrd_pnts];
    likely::TriCubicInterpolator tc_intrpltr(dgrd_vls_array, spcng, dngrd_pnts);
    int nqry_pnts = qry_pnts.size()/sdim; // number of query points
    for(int i=0; i < nqry_pnts; i += sdim)
      qry_vls.push_back(tc_intrpltr(qry_pnts[i],      // x coordinate
                                    qry_pnts[i+1],    // y coordinate
                                    qry_pnts[i+2]));  // z coordinate
  }

  return qry_vls;
}

}  // namespace tbslas
