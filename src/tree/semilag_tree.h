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

#include <profile.hpp>

#include "semilag/semilag.h"
#include "tree/node_field_functor.h"
#include "tree/utils_tree.h"
#include "utils/common.h"

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
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("SolveSemilagTree", &sim_config->comm, false,5);
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
  int num_points_per_node = (cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);

  std::vector<RealType> points_pos_all_nodes;
  tbslas::CollectChebTreeGridPoints(tree_next, points_pos_all_nodes);

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
  NodeType* n_next = tree_next.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf()) break;
    n_next = tree_next.PostorderNxt(n_next);
  }

  int tree_next_node_counter = 0;
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
  pvfmm::Profile::Toc();
}

template <class TreeType>
void SolveSemilagInSitu(TreeType& tvel_curr,
                        TreeType& tree_curr,
                        const int timestep,
                        const typename TreeType::Real_t dt,
                        int num_rk_step = 1) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("SolveSemilagTree", &sim_config->comm, false,5);
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
  int num_points_per_node = (cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);

  std::vector<RealType> points_pos_all_nodes;
  tbslas::CollectChebTreeGridPoints(tree_curr, points_pos_all_nodes);

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
  NodeType* n_next = tree_curr.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf()) break;
    n_next = tree_curr.PostorderNxt(n_next);
  }

  int tree_next_node_counter = 0;
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
    n_next = tree_curr.PostorderNxt(n_next);
  }
  pvfmm::Profile::Toc();
}

template <class TreeType>
void
RunSemilagSimulation(TreeType* vel_tree,
                     TreeType* con_tree_curr,
                     TreeType* con_tree_next,
                     TreeType** result,
                     int tn,
                     typename TreeType::Real_t dt,
                     int num_rk_step = 1,
                     bool adaptive = true,
                     bool save = true) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  int myrank;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myrank);

  // simulation parameters
  SimConfig* sim_param = tbslas::SimConfigSingleton::Instance();

  //////////////////////////////////////////////////////////////////////
  // NEXT STEP TREE
  //////////////////////////////////////////////////////////////////////
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
    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             sim_param->vtk_filename_format.c_str(),
             tbslas::get_result_dir().c_str(),
             sim_param->vtk_filename_prefix.c_str(),
             sim_param->vtk_filename_variable.c_str(),
             tstep);
    tconc_curr->Write2File(out_name_buffer, sim_param->vtk_order);
  }

  pvfmm::Profile::Tic("RunSemilag", &sim_param->comm, false,5);
  for (tstep = 1; tstep < tn+1; tstep++) {
    if(!myrank) {
      printf("============================\n");
      printf("dt: %f tstep: %d time: %f\n", dt, tstep, dt*tstep);
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
      pvfmm::Profile::Tic("RefineTree", &sim_param->comm, false,5);
      tconc_next->RefineTree();
      pvfmm::Profile::Toc();

      // prepare the next step tree
      tbslas::SyncTreeRefinement(*tconc_next, *tconc_curr);
    }

    tbslas::swap_pointers(&tconc_curr, &tconc_next);

    // save current time step data
    *result = tconc_curr;
    if (save) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_param->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_param->vtk_filename_prefix.c_str(),
               sim_param->vtk_filename_variable.c_str(),
               tstep);
      (*result)->Write2File(out_name_buffer, sim_param->vtk_order);
    }
  }  // end of for
  pvfmm::Profile::Toc();
}

template <class TreeType>
void
RunSemilagSimulationInSitu(TreeType* vel_tree,
                           TreeType* con_tree_curr,
                           int tn,
                           typename TreeType::Real_t dt,
                           int num_rk_step = 1,
                           bool adaptive = true,
                           bool save     = true) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  int myrank;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myrank);

    // simulation parameters
  SimConfig* sim_param = tbslas::SimConfigSingleton::Instance();

  // set the input_fn to NULL -> needed for adaptive refinement
  std::vector<NodeType*>  ncurr_list = con_tree_curr->GetNodeList();
  for(int i = 0; i < ncurr_list.size(); i++) {
    ncurr_list[i]->input_fn = NULL;
  }

  //////////////////////////////////////////////////////////////////////
  // TIME STEPPING
  //////////////////////////////////////////////////////////////////////
  int tstep = 0;
  TreeType* tconc_curr = con_tree_curr;
  // save current time step data
  char out_name_buffer[300];
  if (save) {
    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             sim_param->vtk_filename_format.c_str(),
             tbslas::get_result_dir().c_str(),
             sim_param->vtk_filename_prefix.c_str(),
             sim_param->vtk_filename_variable.c_str(),
             tstep);
    tconc_curr->Write2File(out_name_buffer, sim_param->vtk_order);
  }

  pvfmm::Profile::Tic("RunSemilag", &sim_param->comm, false,5);
  for (tstep = 1; tstep < tn+1; tstep++) {
    if(!myrank) {
      printf("============================\n");
      printf("dt: %f tstep: %d time: %f\n", dt, tstep, dt*tstep);
      printf("============================\n");
    }
    tbslas::SolveSemilagInSitu<TreeType>(*vel_tree,
                                         *tconc_curr,
                                         tstep,
                                         dt,
                                         num_rk_step);

    if (adaptive) {
      // refine the tree according to the computed values
      pvfmm::Profile::Tic("RefineTree", &sim_param->comm, false, 5);
      tconc_curr->RefineTree();
      pvfmm::Profile::Toc();
    }

    // save current time step data
    if (save) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_param->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_param->vtk_filename_prefix.c_str(),
               sim_param->vtk_filename_variable.c_str(),
               tstep);
      tconc_curr->Write2File(out_name_buffer, sim_param->vtk_order);
    }
  }  // end of for
  pvfmm::Profile::Toc();
}

}  // namespace tbslas

#endif  // SRC_TREE_SEMILAG_TREE_H_
