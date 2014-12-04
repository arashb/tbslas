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

#ifndef SRC_TREE_TREE_UTILS_H_
#define SRC_TREE_TREE_UTILS_H_

#include <mpi.h>
#include <vector>
#include <string>

#include "utils/profile.h"
#include "utils/common.h"

#include <mortonid.hpp>
#include <tree.hpp>
#include <parUtils.h>

namespace tbslas {

template<typename TreeType>
void clone_tree(TreeType& tree_in,
                TreeType& tree_out,
                int data_dof) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  Profile<double>::Tic("clone_tree",false,5);
  typename NodeType::NodeData tree_data;
  tree_data.dim       = tree_in.Dim();
  tree_data.max_depth = MAX_DEPTH;
  tree_data.cheb_deg  = tree_in.RootNode()->ChebDeg();

  //Set input function pointer
  tree_data.input_fn = tbslas::dummy_fn<RealType>;
  tree_data.data_dof = data_dof;
  tree_data.tol      = 1;

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
  Profile<double>::Tic("Initialize",false,5);
  tree_out.Initialize(&tree_data);
  Profile<double>::Toc();

  //2:1 Balancing
  Profile<double>::Tic("Balance21",false,5);
  tree_out.Balance21(pvfmm::FreeSpace);
  Profile<double>::Toc();

  //Redistribute nodes.
  Profile<double>::Tic("RedistNodes",false,5);
  tree_out.RedistNodes();
  Profile<double>::Toc();

  Profile<double>::Toc();
}

template <typename TreeType,
          typename InputFunction>
void construct_tree(const size_t N,
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

  Profile<double>::Tic("construct_tree",false,5);
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
  Profile<double>::Tic("Initialize",false,5);
  tree.Initialize(&tree_data);
  Profile<double>::Toc();

  Profile<double>::Tic("RefineTree",false,5);
  tree.RefineTree();
  Profile<double>::Toc();

  Profile<double>::Tic("Balance21",false,5);
  tree.Balance21(pvfmm::FreeSpace);
  Profile<double>::Toc();

  Profile<double>::Tic("RedistNodes",false,5);
  tree.RedistNodes();
  Profile<double>::Toc();

  Profile<double>::Toc();
}

template<typename TreeType,
         typename InputFunction>
void
init_tree(TreeType& tree,
          InputFunction input_fn,
          int data_dof) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  Profile<double>::Tic("init_tree",false,5);
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
  Profile<double>::Toc();
}

// in case of the multi-dimensional values
// maximum value of all dimensions together will be returend
// and not the maximum norm value.
template<typename TreeType>
void
max_tree_values(TreeType& tree,
                typename TreeType::Real_t& max_value,
                int& max_depth) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;
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
void
sync_tree_refinement(TreeType& tree_in,
                     TreeType& tree_out) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;

  int np, myrank;
  MPI_Comm_size(*tree_in.Comm(), &np);
  MPI_Comm_rank(*tree_in.Comm(), &myrank);

  std::vector<pvfmm::MortonId> mins_in;
  GetTreeMortonIdMins(tree_in, mins_in);

  std::vector<pvfmm::MortonId> mins_out;
  GetTreeMortonIdMins(tree_out, mins_out);

  for (int ncnt = 0; ncnt < mins_in.size(); ncnt++) {
    if (myrank != np-1) {
      if(mins_in[ncnt] >= mins_out[myrank] &&
         mins_in[ncnt] < mins_out[myrank+1]) {
        tree_out.FindNode(mins_in[ncnt], true, NULL);
      }
    } else {
      if(mins_in[ncnt] >= mins_out[myrank])
        tree_out.FindNode(mins_in[ncnt], true, NULL);
    }
  }

  // RedistNodes(MortionID_array)
  tree_out.RedistNodes(&mins_in[myrank]);
}

}  // namespace tbslas
#endif  // SRC_TREE_TREE_UTILS_H_
