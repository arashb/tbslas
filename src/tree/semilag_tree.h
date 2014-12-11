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

#ifndef SRC_TREE_SEMILAG_TREE_H_
#define SRC_TREE_SEMILAG_TREE_H_

#include <vector>
#include <string>

#include "semilag/semilag.h"

#include "tree/node_field_functor.h"
#include "tree/utils_tree.h"

#include "utils/common.h"
#include "utils/profile.h"

namespace tbslas {

template<class RealType>
struct SimParam {

  SimParam():
      vtk_order(14),
      vtk_filename_format("%s/values_%d_") {
  }
  // time stepping
  int total_num_timestep;
  RealType dt;

  // semi-lagrangian solver
  int num_rk_step;

  // OUTPUT
  int vtk_order;
  std::string vtk_filename_format;
};

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

template <class TreeType>
void
RunSemilagSimulation(TreeType* vel_tree,
                     TreeType* con_tree_curr,
                     struct tbslas::SimParam<typename TreeType::Real_t>* sim_param,
                     bool adaptive = true,
                     bool save = true) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  int myrank;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myrank);

  // simulation parameters
  int tn = sim_param->total_num_timestep;
  RealType dt = sim_param->dt;
  int num_rk_step = sim_param->num_rk_step;

  //////////////////////////////////////////////////////////////////////
  // NEXT STEP TREE
  //////////////////////////////////////////////////////////////////////
  // clone tree
  TreeType* con_tree_next = new TreeType(comm);
  tbslas::CloneTree<TreeType>(*con_tree_curr, *con_tree_next, 1);

  // set the input_fn to NULL -> needed for adaptive refinement
  std::vector<NodeType*>  ncurr_list = con_tree_curr->GetNodeList();
  for(int i = 0; i < ncurr_list.size(); i++) {
    ncurr_list[i]->input_fn = NULL;
  }

  std::vector<NodeType*>  nnext_list = con_tree_next->GetNodeList();
  for(int i = 0; i < nnext_list.size(); i++) {
    nnext_list[i]->input_fn = NULL;
  }

  //////////////////////////////////////////////////////////////////////
  // TIME STEPPING
  //////////////////////////////////////////////////////////////////////
  int tstep = 0;
  TreeType* tconc_curr = con_tree_curr;
  TreeType* tconc_next = con_tree_next;
  // save current time step data
  char out_name_buffer[300];
  if (save) {
    snprintf(out_name_buffer, sizeof(out_name_buffer),
             sim_param->vtk_filename_format.c_str(),
             tbslas::get_result_dir().c_str(), tstep);
    tconc_curr->Write2File(out_name_buffer, sim_param->vtk_order);
  }

  for (tstep = 1; tstep < tn+1; tstep++) {
    if(!myrank) {
      printf("============================\n");
      printf("dt: %f tstep: %d\n", dt, tstep);
      printf("============================\n");
    }
    tbslas::SolveSemilagTree<TreeType>(*vel_tree,
                                       *tconc_curr,
                                       *tconc_next,
                                       tstep,
                                       dt,
                                       num_rk_step);

    if (adaptive) {
      // refine the tree according to the computed values
      tbslas::Profile<double>::Tic("RefineTree",false,5);
      tconc_next->RefineTree();
      tbslas::Profile<double>::Toc();

      // prepare the next step tree
      tbslas::SyncTreeRefinement(*tconc_next, *tconc_curr);
    }

    // save current time step data
    if (save) {
      snprintf(out_name_buffer, sizeof(out_name_buffer),
               sim_param->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(), tstep);
      tconc_next->Write2File(out_name_buffer, sim_param->vtk_order);
    }

    tbslas::swap_pointers(&tconc_curr, &tconc_next);
  }
  // clean up memory
  delete con_tree_next;
}

}  // namespace tbslas

#endif  // SRC_TREE_SEMILAG_TREE_H_