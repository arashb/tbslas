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

#ifndef SRC_TREE_ADVECT_TREE_SEMILAG_H_
#define SRC_TREE_ADVECT_TREE_SEMILAG_H_

#include <vector>

#include "semilag/semilag.h"
#include "tree/node_field_functor.h"
#include "utils/common.h"
#include "utils/profile.h"

namespace tbslas {

template <typename real_t,
          class NodeType,
          class TreeType>
void advect_tree_semilag(TreeType& tvel_curr,
                         TreeType& tree_curr,
                         TreeType& tree_next,
                         int timestep,
                         real_t dt,
                         int num_rk_step = 1) {
  Profile<double>::Tic("advect_tree_semilag",false,5);
  NodeType* n_curr = tree_curr.PostorderFirst();
  NodeType* n_next = tree_next.PostorderFirst();
  int data_dof = n_curr->DataDOF();
  int cheb_deg = n_curr->ChebDeg();
  int sdim     = tree_curr.Dim();

  // compute chebychev points positions on the fly
  std::vector<real_t> cheb_pos = pvfmm::cheb_nodes<real_t>(cheb_deg, sdim);
  int num_points               = cheb_pos.size()/sdim;

  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
    n_curr = tree_curr.PostorderNxt(n_curr);
  }
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf()) break;
    n_next = tree_next.PostorderNxt(n_next);
  }

  while (n_curr != NULL && n_next != NULL) {
    if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
      real_t length      = static_cast<real_t>(std::pow(0.5, n_curr->Depth()));
      real_t* node_coord = n_curr->Coord();

      // TODO: figure out a way to optimize this part.
      std::vector<real_t> points_pos(cheb_pos.size());
      // scale the cheb points
      for (int i = 0; i < num_points; i++) {
        points_pos[i*sdim+0] = node_coord[0] + length * cheb_pos[i*sdim+0];
        points_pos[i*sdim+1] = node_coord[1] + length * cheb_pos[i*sdim+1];
        points_pos[i*sdim+2] = node_coord[2] + length * cheb_pos[i*sdim+2];
      }

      std::vector<real_t> points_val(num_points*data_dof);
      tbslas::semilag_rk2(tbslas::NodeFieldFunctor<real_t, TreeType>(&tvel_curr),
                          tbslas::NodeFieldFunctor<real_t, TreeType>(&tree_curr),
                          points_pos,
                          sdim,
                          timestep,
                          dt,
                          num_rk_step,
                          points_val);

      pvfmm::cheb_approx<real_t, real_t>(points_val.data(),
                                         cheb_deg,
                                         data_dof,
                                         &(n_next->ChebData()[0])
                                         );
    }
    n_curr = tree_curr.PostorderNxt(n_curr);
    n_next = tree_next.PostorderNxt(n_next);
  }
  Profile<double>::Toc();
}

}  // namespace tbslas

#endif  // SRC_TREE_ADVECT_TREE_SEMILAG_H_
