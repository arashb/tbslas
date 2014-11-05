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

#include <vector>

#include "utils/profile.h"
#include "utils/common.h"

namespace tbslas {

template<typename real_t,
         typename NodeType,
         typename TreeType>
void clone_tree(TreeType& tree_in,
                TreeType& tree_out,
                int data_dof) {
  Profile<double>::Tic("clone_tree",false,5);
  typename NodeType::NodeData tree_data;
  tree_data.dim       = tree_in.Dim();
  tree_data.max_depth = MAX_DEPTH;
  tree_data.cheb_deg  = tree_in.RootNode()->ChebDeg();

  //Set input function pointer
  tree_data.input_fn = tbslas::dummy_fn<real_t>;
  tree_data.data_dof = data_dof;
  tree_data.tol      = 1;

  //Set source coordinates.
  std::vector<real_t> pt_coord;
  NodeType* node =
      static_cast< NodeType* >(tree_in.PreorderFirst());

  while(node!=NULL){
    if(node->IsLeaf() && !node->IsGhost()){
      real_t* c=node->Coord();
      real_t s=pow(0.5,node->Depth()+1);
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

template <typename real_t,
          typename NodeType,
          typename TreeType,
          typename InputFunction>
void construct_tree(const size_t N,
                    const size_t M,
                    const int cheb_deg,
                    const int depth,
                    const bool adap,
                    const real_t tol,
                    const MPI_Comm& comm,
                    const InputFunction input_fn,
                    const int data_dof,
                    TreeType& tree) {
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
  std::vector<real_t> pt_coord;
  pt_coord           = tbslas::point_distrib<real_t>(UnifGrid,N,comm);
  tree_data.max_pts  = M; // Points per octant.
  tree_data.pt_coord = pt_coord;

  //initialize with input data.
  Profile<double>::Tic("Initialize",false,5);
  tree.Initialize(&tree_data);
  Profile<double>::Toc();

  Profile<double>::Tic("Initialize",false,5);
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

template<typename real_t,
         typename NodeType,
         typename TreeType,
         typename InputFunction>
void
init_tree(TreeType& tree,
          InputFunction input_fn,
          int data_dof) {
  Profile<double>::Tic("init_tree",false,5);
  NodeType* n_curr = tree.PostorderFirst();
  int cheb_deg     = n_curr->ChebDeg();
  int sdim         = tree.Dim();

  // compute chebychev points positions on the fly
  std::vector<real_t> cheb_pos = pvfmm::cheb_nodes<real_t>(cheb_deg, sdim);
  int num_points               = cheb_pos.size()/sdim;

  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf()) break;
    n_curr = tree.PostorderNxt(n_curr);
  }

  while (n_curr != NULL) {
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
      input_fn(points_pos.data(), num_points, points_val.data());

      pvfmm::cheb_approx<real_t, real_t>(points_val.data(),
                                         cheb_deg,
                                         data_dof,
                                         &(n_curr->ChebData()[0])
                                         );
    }
    n_curr = tree.PostorderNxt(n_curr);
  }
  Profile<double>::Toc();
}

}  // namespace tbslas
#endif // SRC_TREE_TREE_UTILS_H_
