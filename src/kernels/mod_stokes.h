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

#ifndef SRC_KERNEL_MOD_STOKES_H_
#define SRC_KERNEL_MOD_STOKES_H_

#include <utility>
#include <cmath>
#include <string>

#include <math_utils.hpp>

extern double TBSLAS_MOD_STOKES_ALPHA;
extern double TBSLAS_MOD_STOKES_DIFF_COEFF;

char mod_stokes_kernel_name[256];

namespace tbslas {

  /**
   * \brief Green's function for the Stokes's equation. Kernel tensor
   * dimension = 3x3.
   */
  template <class T>
    void modified_stokes_vel(T* r_src, int src_cnt, T* v_src_, int dof, T* r_trg, int trg_cnt, T* v_trg, pvfmm::mem::MemoryManager* mem_mgr){
/* #ifndef __MIC__ */
/*     pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof)); */
/* #endif */

    const T mu=TBSLAS_MOD_STOKES_DIFF_COEFF;
    const T invMu = 1.0/mu;
    const T lambda2=TBSLAS_MOD_STOKES_ALPHA/TBSLAS_MOD_STOKES_DIFF_COEFF;
    const T lambda =sqrt(TBSLAS_MOD_STOKES_ALPHA/TBSLAS_MOD_STOKES_DIFF_COEFF);
    const T COEEF = 1.0/(4.0*pvfmm::const_pi<T>()*lambda2);
    /* const T OOEPMU = 1.0/(8.0*pvfmm::const_pi<T>()*mu); */

    for(int t=0;t<trg_cnt;t++){
      for(int i=0;i<dof;i++){
	T p[3]={0,0,0};
	for(int s=0;s<src_cnt;s++){
	  T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
		   r_trg[3*t+1]-r_src[3*s+1],
		   r_trg[3*t+2]-r_src[3*s+2]};
	  T R2 = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
	  if (R2!=0){
	    T invR2 = 1.0/R2;
	    T R     = sqrt(R2);
	    T invR  = 1.0/R;
	    T invR3 = invR2*invR;

	    T explambda = exp(-lambda*R);
	    T fr = COEEF*( invR2*(explambda-1)   + invR*lambda*explambda );
	    T fr2= COEEF*( 2*invR3*(1-explambda) - 2*lambda*invR2*explambda - lambda2*invR*explambda );
	    T G = -fr2 - invR*fr;
	    T H = invR2*fr2 - invR3*fr;

	    T v_src[3]={v_src_[(s*dof+i)*3  ],
			v_src_[(s*dof+i)*3+1],
			v_src_[(s*dof+i)*3+2]};
	    T inner_prod=(v_src[0]*dR[0] +
			  v_src[1]*dR[1] +
			  v_src[2]*dR[2]);
	    p[0] += (G*v_src[0] + H*dR[0]*inner_prod);
	    p[1] += (G*v_src[1] + H*dR[1]*inner_prod);
	    p[2] += (G*v_src[2] + H*dR[2]*inner_prod);
	  }
	}
	v_trg[(t*dof+i)*3+0] += p[0]*invMu;
	v_trg[(t*dof+i)*3+1] += p[1]*invMu;
	v_trg[(t*dof+i)*3+2] += p[2]*invMu;
      }
    }
  }

  template<typename Real_t>
    const char* GetModfiedStokesKernelName(Real_t alpha, Real_t diff) {
    snprintf(mod_stokes_kernel_name, sizeof(mod_stokes_kernel_name),
	     "modified_stokes_vel_alpha_%05.5f_%05.9f",alpha, diff);
    return mod_stokes_kernel_name;
  }

  const pvfmm::Kernel<double> modified_stokes_vel_d =
    pvfmm::BuildKernel<double, modified_stokes_vel>(GetModfiedStokesKernelName(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3)); /* TODO: change the dimension */

  template<class Real_t>
    struct ModifiedStokesKernel{
      inline static const pvfmm::Kernel<Real_t>& velocity();
    };

  template<> const pvfmm::Kernel<double>& ModifiedStokesKernel<double>::velocity(){
    return modified_stokes_vel_d;
  }

}  // namespace tbslas

#endif  // SRC_KERNEL_MOD_STOKES_H_
