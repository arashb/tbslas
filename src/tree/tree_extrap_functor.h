// *************************************************************************
// Copyright (C) 2016 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_TREE_EXTRAP_FUNCTOR_H_
#define SRC_TREE_EXTRAP_FUNCTOR_H_

#include <vector>

#include <pvfmm_common.hpp>
#include <cheb_node.hpp>
#include <profile.hpp>

#include "utils/common.h"
#include "utils/cubic.h"

namespace tbslas {

  template<typename Real_t,
    class Tree_t>
    class FieldExtrapFunctor {

  public:
      explicit FieldExtrapFunctor(Tree_t* tp,
				  Tree_t* tc):
      tp_(tp),
	tc_(tc) {
	  typedef typename Tree_t::Node_t Node_t;
	  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

	  //////////////////////////////////////////////////////////////////////
	  // GET THE TREES PARAMETERS
	  //////////////////////////////////////////////////////////////////////
	  Node_t* n_curr = tp_->PostorderFirst();
	  while (n_curr != NULL) {
	    if(!n_curr->IsGhost() && n_curr->IsLeaf())
	      break;
	    n_curr = tp_->PostorderNxt(n_curr);
	  }
	  data_dof_ = n_curr->DataDOF();
	}

      virtual ~FieldExtrapFunctor() {
      }

      void operator () (const Real_t* query_points_pos,
			int num_points,
			Real_t* out) {
	////////////////////////////////////////////////////////////////////////
	// EXTRAPOLATE IN TIME FOR ALL QUERY POINTS
	////////////////////////////////////////////////////////////////////////
	// ===================================
	// CONSTRUCT THE EVALUATORS
	// ===================================
	tbslas::NodeFieldFunctor<Real_t,Tree_t> tp_evaluator(tp_);
	tbslas::NodeFieldFunctor<Real_t,Tree_t> tc_evaluator(tc_);
	// ===================================
	// EVALUATE AT T^N
	// ===================================
	std::vector<Real_t> tnc_pnts_val;
	tnc_pnts_val.resize(num_points*data_dof_);
	tc_evaluator(query_points_pos, num_points, tnc_pnts_val.data());
	// ===================================
	// EVALUATE AT T^(N-1)
	// ===================================
	std::vector<Real_t> tnp_pnts_val;
	tnp_pnts_val.resize(num_points*data_dof_);
	tp_evaluator(query_points_pos, num_points, tnp_pnts_val.data());
	// ===================================
	// COMBINE AND STORE THE VALUES
	// ===================================
	Real_t ccoeff = 3.0/2;
	Real_t pcoeff = 0.5;

#pragma omp parallel for
	for (int i = 0; i < tnc_pnts_val.size(); i++) {
	  out[i] = ccoeff*tnc_pnts_val[i] - pcoeff*tnp_pnts_val[i] ;
	}
      }

      void update(Tree_t* new_tree, Real_t time) {
	// POP FRONT (DEALLOCATE AND REMOVE) THE FIRST TREE IN THE QUEUE
	delete tp;
	tp_ = tc_;
	tc_ = new_tree;
      }

  
  private:
      Tree_t* tc_;
      Tree_t* tp_;

      int data_dof_;
    };

}      // namespace tbslas
#endif  // SRC_TREE_EXTRAP_FUNCTOR_H_
