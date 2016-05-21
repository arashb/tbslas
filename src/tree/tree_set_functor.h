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

#ifndef SRC_TREE_FIELD_SET_FUNCTOR_H_
#define SRC_TREE_FIELD_SET_FUNCTOR_H_

#include <vector>

#include <pvfmm_common.hpp>
#include <cheb_node.hpp>
#include <profile.hpp>

#include "utils/common.h"
#include "utils/cubic.h"

namespace tbslas {

template<typename Real_t,
         class Tree_t>
class FieldSetFunctor {

 public:
  explicit FieldSetFunctor(std::vector<Tree_t*> field_set_elems,
                           std::vector<Real_t>  field_set_times):
      field_set_elems_(field_set_elems),
      field_set_times_(field_set_times) {
    typedef typename Tree_t::Node_t Node_t;
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

    //////////////////////////////////////////////////////////////////////
    // GET THE TREES PARAMETERS
    //////////////////////////////////////////////////////////////////////
    Node_t* n_curr = field_set_elems[0]->PostorderFirst();
    while (n_curr != NULL) {
      if(!n_curr->IsGhost() && n_curr->IsLeaf())
        break;
      n_curr = field_set_elems[0]->PostorderNxt(n_curr);
    }
    data_dof_ = n_curr->DataDOF();
  }

  virtual ~FieldSetFunctor() {
  }

  void operator () (const Real_t* query_points_pos,
                    int num_points,
                    Real_t query_time,
                    Real_t* out) {
    ////////////////////////////////////////////////////////////////////////
    // INTERPOLATE IN SPACE FOR ALL QUERY POINTS IN CORRESPONDING TIME STEPS
    ////////////////////////////////////////////////////////////////////////
    std::vector<std::vector<Real_t>*> field_set_points_val;
    for (int i = 0 ; i < field_set_elems_.size(); i++) {
      std::vector<Real_t>* pPoints_val = new std::vector<Real_t>(num_points*data_dof_);
      tbslas::NodeFieldFunctor<Real_t,Tree_t> tree_func(field_set_elems_[i]);
      tree_func(query_points_pos, num_points, (*pPoints_val).data());
      field_set_points_val.push_back(pPoints_val);
    }

    //////////////////////////////////////////////////////////////////////
    // INTERPOLATE IN TIME (CUBIC)
    //////////////////////////////////////////////////////////////////////
    for (int i = 0; i < field_set_points_val[0]->size(); i++) {
      Real_t grid_vals[4] = {(*field_set_points_val[0])[i],
                             (*field_set_points_val[1])[i],
                             (*field_set_points_val[2])[i],
                             (*field_set_points_val[3])[i]};
      out[i] = tbslas::CubicInterpPolicy<Real_t>::InterpCubic1D
          (query_time, field_set_times_.data(), grid_vals);
    }

    //////////////////////////////////////////////////////////////////////
    // DEALLOCATTE THE TMP VALUES
    //////////////////////////////////////////////////////////////////////
    for (int i = 0 ; i < field_set_elems_.size(); i++)
      delete field_set_points_val[i];
  }

  void update(Tree_t* new_tree, Real_t time) {
    // POP FRONT (DEALLOCATE AND REMOVE) THE FIRST TREE IN THE QUEUE
    delete field_set_elems_[0];
    field_set_elems_.erase(field_set_elems_.begin());
    field_set_times_.erase(field_set_times_.begin());

    // PUSH BACK THE NEW TREE
    field_set_times_.push_back(time);
    field_set_elems_.push_back(new_tree);
  }


 private:
  std::vector<Real_t>  field_set_times_;
  std::vector<Tree_t*> field_set_elems_;

  int data_dof_;
};

}      // namespace tbslas
#endif  // SRC_TREE_FIELD_SET_FUNCTOR_H_
