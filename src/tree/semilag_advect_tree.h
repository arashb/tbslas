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

#ifndef SRC_TREE_SEMILAG_ADVECT_TREE_H_
#define SRC_TREE_SEMILAG_ADVECT_TREE_H_

#include <vector>

#include "semilag/semilag.h"
#include "tree/tree_common.h"

namespace tbslas {

template <typename real_t>
void semilag_advect_tree(Tree_t<real_t>& tvel_curr,
                         Tree_t<real_t>& tree_curr,
                         Tree_t<real_t>& tree_next) {
  int sdim        = 3;
  int timestep    = 1;
  real_t dt       = 0.5;
  int num_rk_step = 1;

  Node_t<real_t>* n_curr = tree_curr.PostorderFirst();
  Node_t<real_t>* n_next = tree_next.PostorderFirst();
  int data_dof   = n_curr->DataDOF();
  int cheb_deg   = n_curr->ChebDeg();

  // compute chebychev points positions on the fly
  std::vector<real_t> cheb_pos = pvfmm::cheb_nodes<real_t>(cheb_deg, sdim);
  int num_points               = cheb_pos.size()/COORD_DIM;

  while (n_curr != NULL) {
    if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
      real_t length      = static_cast<real_t>(std::pow(0.5, n_curr->Depth()));
      real_t* node_coord = n_curr->Coord();

      printf("NODE: [%f, %f, %f]\n",
             node_coord[0],
             node_coord[1],
             node_coord[2]);
      // TODO: figure out a way to optimize this part.
      std::vector<real_t> points_pos(cheb_pos.size());
      // scale the cheb points
      for (int i = 0; i < num_points; i++) {
        points_pos[i*COORD_DIM+0] = node_coord[0] + length * cheb_pos[i*COORD_DIM+0];
        points_pos[i*COORD_DIM+1] = node_coord[1] + length * cheb_pos[i*COORD_DIM+1];
        points_pos[i*COORD_DIM+2] = node_coord[2] + length * cheb_pos[i*COORD_DIM+2];
      }

      std::vector<real_t> points_val(num_points);
      tbslas::NodeFieldFunctor<real_t> con_field(tree_curr.RootNode());
      tbslas::NodeFieldFunctor<real_t> vel_field(tvel_curr.RootNode());

      slas::semilag_rk2(vel_field,
                        con_field,
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
}

}  // namespace tbslas

#endif  // SRC_TREE_SEMILAG_ADVECT_TREE_H_
