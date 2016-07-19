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

#ifndef SRC_KERNEL_MOD_LAPLACE_H_
#define SRC_KERNEL_MOD_LAPLACE_H_

#include <utility>
#include <cmath>
#include <string>

#include <math_utils.hpp>

extern double TBSLAS_ALPHA;

char kernel_name[256];

namespace tbslas {

// alpha*u - \Delta u = alpha*f
template <class Real_t>
void modified_laplace_poten(Real_t* r_src,
                            int src_cnt,
                            Real_t* v_src,
                            int dof,
                            Real_t* r_trg,
                            int trg_cnt,
                            Real_t* k_out,
                            pvfmm::mem::MemoryManager* mem_mgr) {
  Real_t lambda = sqrt(TBSLAS_ALPHA);
  const Real_t OOFP = 1.0/(4.0*pvfmm::const_pi<Real_t>());
  for (int t=0;t<trg_cnt;t++) {
    for (int i=0;i<dof;i++) {
      Real_t p=0;
      for(int s=0;s<src_cnt;s++) {
        Real_t dX_reg = r_trg[3*t  ]-r_src[3*s  ];
        Real_t dY_reg = r_trg[3*t+1]-r_src[3*s+1];
        Real_t dZ_reg = r_trg[3*t+2]-r_src[3*s+2];
        Real_t invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
        Real_t R = 0;
        if (invR!=0) {
          R = sqrt(invR);
          invR = 1.0/R;
        }
        p += v_src[s*dof+i]*invR*exp(-lambda*R)*TBSLAS_ALPHA;
      }
      k_out[t*dof+i] += p*OOFP;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// CONSTRUCT THE KERNEL NAME
//////////////////////////////////////////////////////////////////////
template<typename Real_t>
const char* GetModfiedLaplaceKernelName(Real_t alpha) {
  snprintf(kernel_name, sizeof(kernel_name),
           "modified_laplace_alpha_%5.2f",alpha);
  return kernel_name;
}

//////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLE FOR KERNEL
//////////////////////////////////////////////////////////////////////
const pvfmm::Kernel<double> modified_laplace_potn_d =
    pvfmm::BuildKernel<double, modified_laplace_poten>
    (GetModfiedLaplaceKernelName(TBSLAS_ALPHA), 3, std::pair<int,int>(1,1));

//////////////////////////////////////////////////////////////////////
// KERNEL STRUCTURE
//////////////////////////////////////////////////////////////////////
template<class Real_t>
struct ModifiedLaplaceKernel{
  inline static const pvfmm::Kernel<Real_t>& potn_ker();
};

template<> const pvfmm::Kernel<double>& ModifiedLaplaceKernel<double>::potn_ker() {
  return modified_laplace_potn_d;
}

}  // namespace tbslas

#endif  // SRC_KERNEL_MOD_LAPLACE_H_
