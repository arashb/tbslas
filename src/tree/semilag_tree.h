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
#include "tree/utils_tree.h"

#include "utils/common.h"
#include "utils/profile.h"

namespace tbslas {

template <class TreeType>
void SolveSemilagTree(TreeType& tvel_curr,
                      TreeType& tree_curr,
                      TreeType& tree_next,
                      const int timestep,
                      const typename TreeType::Real_t dt,
                      int num_rk_step = 1) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  Profile<double>::Tic("SolveSemilagTree",false,5);
  ////////////////////////////////////////////////////////////////////////
  // (1) collect the starting positions of the backward traj computation
  ////////////////////////////////////////////////////////////////////////
  // get the tree parameters from current tree
  NodeType* n_curr = tree_curr.PostorderFirst();
  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf())
      break;
    n_curr = tree_curr.PostorderNxt(n_curr);
  }
  int data_dof = n_curr->DataDOF();
  int cheb_deg = n_curr->ChebDeg();
  int sdim     = tree_curr.Dim();

  // compute chebychev points positions on the fly
  std::vector<RealType> cheb_pos = pvfmm::cheb_nodes<RealType>(cheb_deg, sdim);
  int num_points_per_node        = cheb_pos.size()/sdim;

#ifndef NDEBUG
  CountNumLeafNodes(tree_next);
#endif

  // compute total number of tree leaf nodes
  NodeType* n_next = tree_next.PostorderFirst();
  int num_leaf_nodes = 0;
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      num_leaf_nodes++;
    n_next = tree_next.PostorderNxt(n_next);
  }

  std::vector<RealType> points_pos_all_nodes;
  points_pos_all_nodes.resize(cheb_pos.size()*num_leaf_nodes);

  n_next = tree_next.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      break;
    n_next = tree_next.PostorderNxt(n_next);
  }

  int tree_next_node_counter = 0;
  while (n_next != NULL) {
    if (n_next->IsLeaf() && !n_next->IsGhost()) {
      RealType length      = static_cast<RealType>(std::pow(0.5, n_next->Depth()));
      RealType* node_coord = n_next->Coord();
      // scale the cheb points
      size_t shift = tree_next_node_counter*cheb_pos.size();
      for (int i = 0; i < num_points_per_node; i++) {
        points_pos_all_nodes[shift + i*sdim+0] = node_coord[0] + length * cheb_pos[i*sdim+0];
        points_pos_all_nodes[shift + i*sdim+1] = node_coord[1] + length * cheb_pos[i*sdim+1];
        points_pos_all_nodes[shift + i*sdim+2] = node_coord[2] + length * cheb_pos[i*sdim+2];
      }
      tree_next_node_counter++;
    }
    n_next = tree_next.PostorderNxt(n_next);
  }
  ////////////////////////////////////////////////////////////////////////
  // (2) solve semi-Lagrangian
  ////////////////////////////////////////////////////////////////////////
  int num_points_local_nodes = points_pos_all_nodes.size()/sdim;
  std::vector<RealType> points_val_local_nodes(num_points_local_nodes*data_dof);
  tbslas::SolveSemilagRK2(tbslas::NodeFieldFunctor<RealType,TreeType>(&tvel_curr),
                          tbslas::NodeFieldFunctor<RealType,TreeType>(&tree_curr),
                          points_pos_all_nodes,
                          sdim,
                          timestep,
                          dt,
                          num_rk_step,
                          points_val_local_nodes);

  ////////////////////////////////////////////////////////////////////////
  // (3) set the computed values
  ////////////////////////////////////////////////////////////////////////
  n_next = tree_next.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf()) break;
    n_next = tree_next.PostorderNxt(n_next);
  }

  tree_next_node_counter = 0;
  while (n_next != NULL) {
    if (n_next->IsLeaf() && !n_next->IsGhost()) {
      pvfmm::cheb_approx<RealType, RealType>(
          &points_val_local_nodes[tree_next_node_counter*num_points_per_node*data_dof],
          cheb_deg,
          data_dof,
          &(n_next->ChebData()[0])
                                             );
      tree_next_node_counter++;
    }
    n_next = tree_next.PostorderNxt(n_next);
  }
  Profile<double>::Toc();
}

}  // namespace tbslas

#endif  // SRC_TREE_SOLVE_SEMILAG_TREE_H_
