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
CubicInterpolant<real_t>::interp(std::vector<real_t>& grd_pnts,
                                 std::vector<real_t>& grd_vls,
                                 std::vector<real_t>& qry_pnts,
                                 std::vector<real_t>& qry_vls) {
  int ngp = grd_pnts.size()/3;  // number of grid points
  float dngp = pow(ngp, 1./3);
  real_t spcng = 1.0/(dngp-1);  // spacing

  // FIXME: move the construction of Interpolator
  // to the constructor of this class.
  // Here only initialize it with values.
  int nqp = qry_pnts.size()/3; // number of query points
  for (int vdim_cnt = 0; vdim_cnt < 3; vdim_cnt++) {
    real_t* dgrd_vls_array = &grd_vls[vdim_cnt*ngp];
    likely::TriCubicInterpolator tc_intrpltr(dgrd_vls_array, spcng, dngp);
    for(int i=0; i < nqp; i++)
      qry_vls[i+vdim_cnt*nqp] = tc_intrpltr(qry_pnts[i], qry_pnts[i+nqp],
                                            qry_pnts[i+2*nqp]);
  }
}

}  // namespace tbslas
