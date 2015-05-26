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

#ifndef SRC_TREE_UTILS_TREE_H_
#define SRC_TREE_UTILS_TREE_H_

#include <mpi.h>
#include <vector>
#include <string>

#include <mortonid.hpp>
#include <tree.hpp>
#include <parUtils.h>
#include <profile.hpp>

#include "utils/common.h"
#include "tree/node_field_functor.h"

namespace tbslas {

template<typename TreeType>
void CloneTree(TreeType& tree_in,
               TreeType& tree_out,
               int data_dof
               ) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("CloneTree", &sim_config->comm, false, 5);
  typename NodeType::NodeData tree_data;
  tree_data.dim       = tree_in.Dim();
  tree_data.max_depth = MAX_DEPTH;
  tree_data.cheb_deg  = tree_in.RootNode()->ChebDeg();

  //Set input function pointer
  tree_data.input_fn = tbslas::dummy_fn<RealType>;
  tree_data.data_dof = data_dof;
  tree_data.tol      = tree_in.RootNode()->MaxErr();

  //Set source coordinates.
  std::vector<RealType> pt_coord;
  NodeType* node =
      static_cast< NodeType* >(tree_in.PreorderFirst());

  while (node!=NULL) {
    if (node->IsLeaf() && !node->IsGhost()) {
      RealType* c=node->Coord();
      RealType s=pow(0.5,node->Depth()+1);
      pt_coord.push_back(c[0]+s);
      pt_coord.push_back(c[1]+s);
      pt_coord.push_back(c[2]+s);
    }
    node = static_cast<NodeType*>(tree_in.PreorderNxt(node));
  }
  //Set source coordinates.
  tree_data.max_pts = 1; // Points per octant.
  tree_data.pt_coord=pt_coord;

  //Create Tree and initialize with input data.
  pvfmm::Profile::Tic("Initialize", &sim_config->comm, false,5);
  tree_out.Initialize(&tree_data);
  pvfmm::Profile::Toc();

  //2:1 Balancing
  pvfmm::Profile::Tic("Balance21", &sim_config->comm, false,5);
  tree_out.Balance21(pvfmm::FreeSpace);
  pvfmm::Profile::Toc();

  //Redistribute nodes.
  pvfmm::Profile::Tic("RedistNodes", &sim_config->comm, false,5);
  tree_out.RedistNodes();
  pvfmm::Profile::Toc();

  pvfmm::Profile::Toc();
}

template <typename TreeType,
          typename InputFunction>
void ConstructTree(const size_t N,
                   const size_t M,
                   const int cheb_deg,
                   const int depth,
                   const bool adap,
                   const typename TreeType::Real_t tol,
                   const MPI_Comm& comm,
                   const InputFunction input_fn,
                   const int data_dof,
                   TreeType& tree) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("ConstructTree", &sim_config->comm, false,5);
  //Various parameters.
  typename NodeType::NodeData tree_data;
  tree_data.dim       = COORD_DIM;
  tree_data.max_depth = depth;
  tree_data.cheb_deg  = cheb_deg;

  //Set input function pointer
  tree_data.input_fn  = input_fn;
  tree_data.data_dof  = data_dof;
  tree_data.tol       = tol;

  //Set source coordinates.
  std::vector<RealType> pt_coord;
  pt_coord           = tbslas::point_distrib<RealType>(UnifGrid,N,comm);
  tree_data.max_pts  = M; // Points per octant.
  tree_data.pt_coord = pt_coord;

  //initialize with input data.
  pvfmm::Profile::Tic("Initialize", &sim_config->comm, false, 5);
  tree.Initialize(&tree_data);
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
  tree.RefineTree();
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("Balance21", &sim_config->comm, false, 5);
  tree.Balance21(pvfmm::FreeSpace);
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("RedistNodes", &sim_config->comm, false, 5);
  tree.RedistNodes();
  pvfmm::Profile::Toc();

  pvfmm::Profile::Toc();
}

template<typename TreeType,
         typename InputFunction>
void
InitTree(TreeType& tree,
         InputFunction input_fn,
         int data_dof) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("InitTree", &sim_config->comm, false,5);
  NodeType* n_curr = tree.PostorderFirst();
  int cheb_deg     = n_curr->ChebDeg();
  int sdim         = tree.Dim();

  // compute chebychev points positions on the fly
  std::vector<RealType> cheb_pos = pvfmm::cheb_nodes<RealType>(cheb_deg, sdim);
  int num_points               = cheb_pos.size()/sdim;

  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
    n_curr = tree.PostorderNxt(n_curr);
  }

  while (n_curr != NULL) {
    if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
      RealType length      = static_cast<RealType>(std::pow(0.5, n_curr->Depth()));
      RealType* node_coord = n_curr->Coord();

      // TODO: figure out a way to optimize this part.
      std::vector<RealType> points_pos(cheb_pos.size());
      // scale the cheb points
      for (int i = 0; i < num_points; i++) {
        points_pos[i*sdim+0] = node_coord[0] + length * cheb_pos[i*sdim+0];
        points_pos[i*sdim+1] = node_coord[1] + length * cheb_pos[i*sdim+1];
        points_pos[i*sdim+2] = node_coord[2] + length * cheb_pos[i*sdim+2];
      }

      std::vector<RealType> points_val(num_points*data_dof);
      input_fn(points_pos.data(), num_points, points_val.data());

      pvfmm::cheb_approx<RealType, RealType>(points_val.data(),
                                             cheb_deg,
                                             data_dof,
                                             &(n_curr->ChebData()[0])
                                             );
    }
    n_curr = tree.PostorderNxt(n_curr);
  }
  pvfmm::Profile::Toc();
}

// in case of the multi-dimensional values
// maximum value of all dimensions together will be returend
// and not the maximum norm value.
template<typename TreeType>
void
GetMaxTreeValues(TreeType& tree,
                 typename TreeType::Real_t& max_value,
                 int& max_depth) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  NodeType* n_curr = tree.PostorderFirst();
  int cheb_deg     = n_curr->ChebDeg();
  int sdim         = tree.Dim();
  int data_dof     = n_curr->DataDOF();

  double node_max_value = 0;
  int nod_max_depth = 0;

  double lcl_max_values[] = {0,0};
  double glb_max_values[] = {0,0};

  // compute chebychev points positions on the fly
  std::vector<RealType> pos_x = pvfmm::cheb_nodes<RealType>(cheb_deg, 1);
  std::vector<RealType> pos_y = pos_x;
  std::vector<RealType> pos_z = pos_x;

  int size_1d = pos_x.size();
  int num_points = static_cast<int>(std::pow(size_1d,3));

  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
    n_curr = tree.PostorderNxt(n_curr);
  }

  while (n_curr != NULL) {
    if (n_curr->IsLeaf() && !n_curr->IsGhost()) {
      RealType length      = static_cast<RealType>(std::pow(0.5, n_curr->Depth()));
      RealType* node_coord = n_curr->Coord();
      std::vector<RealType> lcl_pos_x(size_1d);
      std::vector<RealType> lcl_pos_y(size_1d);
      std::vector<RealType> lcl_pos_z(size_1d);

      // scale the cheb points
      for (int i = 0; i < size_1d; i++) {
        lcl_pos_x[i] = node_coord[0] + length * pos_x[i];
        lcl_pos_y[i] = node_coord[1] + length * pos_y[i];
        lcl_pos_z[i] = node_coord[2] + length * pos_z[i];
      }

      // node maximum value
      std::vector<RealType> points_val(num_points*data_dof);
      n_curr->ReadVal(lcl_pos_x,
                      lcl_pos_y,
                      lcl_pos_z,
                      points_val.data());

      node_max_value = *std::max_element(points_val.begin(), points_val.end());
      if (node_max_value > lcl_max_values[0])
        lcl_max_values[0] = node_max_value;

      // node maximum depth
      nod_max_depth = n_curr->Depth();
      if (nod_max_depth > lcl_max_values[1])
        lcl_max_values[1] = nod_max_depth;
    }
    n_curr = tree.PostorderNxt(n_curr);
  }
  MPI_Allreduce(&lcl_max_values, &glb_max_values, 2,
                MPI_DOUBLE, MPI_MAX, *(tree.Comm()));

  max_value = glb_max_values[0];
  max_depth = static_cast<int>(glb_max_values[1]);
}

template<typename TreeType>
void
GetTreeMortonIdMins(TreeType& tree,
                    std::vector<pvfmm::MortonId>& mins) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;

  NodeType* n=tree.PreorderFirst();
  while (n != NULL) {
    if(!n->IsGhost() && n->IsLeaf()) break;
    n=tree.PreorderNxt(n);
  }
  ASSERT_WITH_MSG(n!=NULL,"No non-ghost nodes found on this process.");

  pvfmm::MortonId my_min;
  my_min = n->GetMortonId();

  int np;
  MPI_Comm_size(*tree.Comm(),&np);
  mins.resize(np);

  MPI_Allgather(&my_min , 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                &mins[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), *tree.Comm());

}

template<typename TreeType>
int
CountNumLeafNodes(TreeType& tree) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;

  int np, myrank;
  MPI_Comm_size(*tree.Comm(), &np);
  MPI_Comm_rank(*tree.Comm(), &myrank);

  int num_leaf_nodes = 0;
  // compute total number of tree leaf nodes
  NodeType* n_next = tree.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      num_leaf_nodes++;
    n_next = tree.PostorderNxt(n_next);
  }
  int total_num_leaf_nodes = 0;
  MPI_Reduce(&num_leaf_nodes, &total_num_leaf_nodes, 1, MPI_INT,
             MPI_SUM, 0, *tree.Comm());

  return total_num_leaf_nodes;
}

template<typename NodeType>
void
MergeNodeRefinement(NodeType* node_in,
                   NodeType* node_out) {
  if(node_in->IsGhost()) return;

  if(!node_in->IsLeaf()) {
    node_out->Subdivide();
    int n_child = 1UL<< node_in->Dim();
    for (int k = 0; k < n_child; k++) {
      MergeNodeRefinement(dynamic_cast<NodeType*>(node_in->Child(k)),
                          dynamic_cast<NodeType*>(node_out->Child(k)));
    }
  }
}

template<typename TreeType>
void
MergeTreeRefinement(TreeType& tree_in,
                    TreeType& tree_out) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("MergeTreeRefinement", &sim_config->comm, false,5);
  int np, myrank;
  MPI_Comm_size(*tree_in.Comm(), &np);
  MPI_Comm_rank(*tree_in.Comm(), &myrank);

  std::vector<pvfmm::MortonId> mins_in;
  GetTreeMortonIdMins(tree_in, mins_in);

  std::vector<pvfmm::MortonId> mins_out;
  GetTreeMortonIdMins(tree_out, mins_out);

  size_t range[2]={0,np};
  range[0]=std::lower_bound(&mins_in[0], &mins_in[np], mins_out[myrank  ])-&mins_in[0];
  if(myrank<np-1) range[1]=std::lower_bound(&mins_in[0], &mins_in[np], mins_out[myrank+1])-&mins_in[0];
  for(size_t i=range[0];i<range[1];i++) tree_out.FindNode(mins_in[i],true,NULL);

  tree_out.RedistNodes(&mins_in[myrank]);

  NodeType* node_in  = tree_in.PreorderFirst();
  NodeType* node_out = tree_out.PreorderFirst();
  MergeNodeRefinement<NodeType>(node_in, node_out);

  pvfmm::Profile::Toc();
}

template<typename NodeType>
void
SyncNodeRefinement(NodeType* node_in,
                   NodeType* node_out) {
  if(node_in->IsGhost()) return;

  if(node_in->IsLeaf()) {
    node_out->Truncate();
  } else {
    node_out->Subdivide();
    int n_child = 1UL<< node_in->Dim();
    for (int k = 0; k < n_child; k++) {
      SyncNodeRefinement(dynamic_cast<NodeType*>(node_in->Child(k)),
                         dynamic_cast<NodeType*>(node_out->Child(k)));
    }
  }
}

template<typename TreeType>
void
SyncTreeRefinement(TreeType& tree_in,
                   TreeType& tree_out) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("SyncTreeRefinement", &sim_config->comm, false,5);
  int np, myrank;
  MPI_Comm_size(*tree_in.Comm(), &np);
  MPI_Comm_rank(*tree_in.Comm(), &myrank);

  std::vector<pvfmm::MortonId> mins_in;
  GetTreeMortonIdMins(tree_in, mins_in);

  std::vector<pvfmm::MortonId> mins_out;
  GetTreeMortonIdMins(tree_out, mins_out);

  size_t range[2]={0,np};
  range[0]=std::lower_bound(&mins_in[0], &mins_in[np], mins_out[myrank  ])-&mins_in[0];
  if(myrank<np-1) range[1]=std::lower_bound(&mins_in[0], &mins_in[np], mins_out[myrank+1])-&mins_in[0];
  for(size_t i=range[0];i<range[1];i++) tree_out.FindNode(mins_in[i],true,NULL);

  tree_out.RedistNodes(&mins_in[myrank]);

  NodeType* node_in  = tree_in.PreorderFirst();
  NodeType* node_out = tree_out.PreorderFirst();
  SyncNodeRefinement<NodeType>(node_in, node_out);

  pvfmm::Profile::Toc();
}

template<typename TreeType>
void CollectChebTreeGridPoints(TreeType& tree,
                           std::vector<typename TreeType::Real_t>& grid_points) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("CollecGridPoints", &sim_config->comm, false,5);
  int cheb_deg = tree.RootNode()->ChebDeg();
  int sdim     = tree.Dim();

  // compute chebychev points positions on the fly
  std::vector<RealType> cheb_pos = pvfmm::cheb_nodes<RealType>(cheb_deg, sdim);
  int num_points_per_node        = cheb_pos.size()/sdim;

  // compute total number of tree leaf nodes
  NodeType* n_next = tree.PostorderFirst();
  int num_leaf_nodes = 0;
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      num_leaf_nodes++;
    n_next = tree.PostorderNxt(n_next);
  }

  // std::vector<RealType> grid_points;
  grid_points.resize(cheb_pos.size()*num_leaf_nodes);

  n_next = tree.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      break;
    n_next = tree.PostorderNxt(n_next);
  }

  int tree_node_counter = 0;
  while (n_next != NULL) {
    if (n_next->IsLeaf() && !n_next->IsGhost()) {
      RealType length      = static_cast<RealType>(std::pow(0.5, n_next->Depth()));
      RealType* node_coord = n_next->Coord();
      // scale the cheb points
      size_t shift = tree_node_counter*cheb_pos.size();
      for (int i = 0; i < num_points_per_node; i++) {
        grid_points[shift + i*sdim+0] = node_coord[0] + length * cheb_pos[i*sdim+0];
        grid_points[shift + i*sdim+1] = node_coord[1] + length * cheb_pos[i*sdim+1];
        grid_points[shift + i*sdim+2] = node_coord[2] + length * cheb_pos[i*sdim+2];
      }
      tree_node_counter++;
    }
    n_next = tree.PostorderNxt(n_next);
  }
  pvfmm::Profile::Toc();
}

template<typename TreeType>
void CollectTreeGridPoints(TreeType& tree,
                           std::vector<typename TreeType::Real_t> node_pos,
                           std::vector<typename TreeType::Real_t>& grid_points) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("CollectGridPoints", &sim_config->comm, false,5);
  int cheb_deg = tree.RootNode()->ChebDeg();
  int sdim     = tree.Dim();

  // compute chebychev points positions on the fly
  // std::vector<RealType> node_pos = pvfmm::cheb_nodes<RealType>(cheb_deg, sdim);
  int num_points_per_node        = node_pos.size()/sdim;

  // compute total number of tree leaf nodes
  NodeType* n_next = tree.PostorderFirst();
  int num_leaf_nodes = 0;
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      num_leaf_nodes++;
    n_next = tree.PostorderNxt(n_next);
  }

  // std::vector<RealType> grid_points;
  grid_points.resize(node_pos.size()*num_leaf_nodes);

  n_next = tree.PostorderFirst();
  while (n_next != NULL) {
    if(!n_next->IsGhost() && n_next->IsLeaf())
      break;
    n_next = tree.PostorderNxt(n_next);
  }

  int tree_node_counter = 0;
  while (n_next != NULL) {
    if (n_next->IsLeaf() && !n_next->IsGhost()) {
      RealType length      = static_cast<RealType>(std::pow(0.5, n_next->Depth()));
      RealType* node_coord = n_next->Coord();
      // scale the cheb points
      size_t shift = tree_node_counter*node_pos.size();
      for (int i = 0; i < num_points_per_node; i++) {
        grid_points[shift + i*sdim+0] = node_coord[0] + length * node_pos[i*sdim+0];
        grid_points[shift + i*sdim+1] = node_coord[1] + length * node_pos[i*sdim+1];
        grid_points[shift + i*sdim+2] = node_coord[2] + length * node_pos[i*sdim+2];
      }
      tree_node_counter++;
    }
    n_next = tree.PostorderNxt(n_next);
  }
  pvfmm::Profile::Toc();
}
template<typename TreeType,
         typename AnalyticalFunctor>
void
ComputeTreeError(TreeType& tree,
                 AnalyticalFunctor anal_func,
                 std::vector<typename TreeType::Real_t> grid_points,
                 typename TreeType::Real_t& l_inf_error,
                 typename TreeType::Real_t& l_two_error) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Tic("ComputeTreeError",&sim_config->comm, false,5);
  int data_dof = 1;
  int sdim     = tree.Dim();
  int num_points = grid_points.size()/sdim;

  // analytical values
  std::vector<RealType> points_val_analytical(num_points*data_dof);
  anal_func(grid_points.data(), num_points, points_val_analytical.data());

  // tree values
  std::vector<RealType> points_val(num_points*data_dof);
  auto tree_func = tbslas::NodeFieldFunctor<RealType, TreeType>(&tree);
  tree_func(grid_points.data(), num_points, points_val.data());

  // local absolute error vector
  RealType local_max = 0;
  RealType local_sum_2 = 0;
  std::vector<RealType> abs_error(points_val.size());
  for (int i = 0; i < points_val.size(); i++) {
    abs_error[i] = std::abs(points_val[i] - points_val_analytical[i]);
    local_sum_2 += abs_error[i]*abs_error[i];
    if (abs_error[i] > local_max)
      local_max = abs_error[i];
  }
  size_t lcl_num_points = grid_points.size()/sdim;
  size_t glb_num_points = 0;
  MPI_Allreduce(&lcl_num_points, &glb_num_points, 1,
                MPI_INT, MPI_SUM, *(tree.Comm()));
  // compute the error vector's norm
  RealType global_max;
  MPI_Allreduce(&local_max, &global_max, 1,
                MPI_DOUBLE, MPI_MAX, *(tree.Comm()));
  RealType global_sum_2;
  MPI_Allreduce(&local_sum_2, &global_sum_2, 1,
                MPI_DOUBLE, MPI_SUM, *(tree.Comm()));
  l_inf_error = global_max;
  l_two_error = sqrt(global_sum_2/glb_num_points);
  pvfmm::Profile::Toc();
}

}  // namespace tbslas
#endif  // SRC_TREE_UTILS_TREE_H_
