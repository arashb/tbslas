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

#include "tree/tree_common.h"

namespace tbslas {

template<typename real_t,
         class NodeType,
         class TreeType>
void clone_tree(TreeType& tree_in,
                TreeType& tree_out ) {
  typename NodeType::NodeData tree_data;
  tree_data.dim       = tree_in.Dim();
  tree_data.max_depth = MAX_DEPTH;
  tree_data.cheb_deg  = tree_in.RootNode()->ChebDeg();

  //Set input function pointer
  tree_data.input_fn = tbslas::dummy_fn<real_t>;
  tree_data.data_dof = tree_in.RootNode()->DataDOF();
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
  tree_out.Initialize(&tree_data);

  //2:1 Balancing
  tree_in.Balance21(pvfmm::FreeSpace);

  //Redistribute nodes.
  tree_in.RedistNodes();
}

template <typename real_t, typename InputFunction>
void construct_tree(const size_t N,
                    const size_t M,
                    const int cheb_deg,
                    const int depth,
                    const bool adap,
                    const real_t tol,
                    const MPI_Comm& comm,
                    const InputFunction input_fn,
                    const int data_dof,
                    Tree_t<real_t>& tree) {
  //Various parameters.
  typename Node_t<real_t>::NodeData tree_data;
  tree_data.dim       = COORD_DIM;
  tree_data.max_depth = depth;
  tree_data.cheb_deg  = cheb_deg;

  //Set input function pointer
  tree_data.input_fn  = input_fn;
  tree_data.data_dof  = data_dof;
  tree_data.tol       = tol;

  //Set source coordinates.
  std::vector<real_t> pt_coord;
  pt_coord= tbslas::point_distrib<real_t>(UnifGrid,N,comm);
  tree_data.max_pts  = M; // Points per octant.
  tree_data.pt_coord = pt_coord;

  //initialize with input data.
  tree.Initialize(&tree_data);
  tree.RefineTree();
  tree.Balance21(pvfmm::FreeSpace);
  tree.RedistNodes();
}

template<typename real_t,
         class NodeType,
         class TreeType,
         typename InputFunction>
void
init_tree(TreeType& tree,
          InputFunction input_fn) {
  NodeType* n_curr = tree.PostorderFirst();
  int data_dof = n_curr->DataDOF();
  int cheb_deg = n_curr->ChebDeg();
  int sdim     = tree.Dim();

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
}

template<class TreeType>
void swap_trees_pointers(TreeType** ta,
                         TreeType** tb) {
  TreeType* tmp = *ta;
  *ta = *tb;
  *tb = tmp;
}

template<typename real_t>
std::vector<int>
isOutside(Node_t<real_t>* n,
          const std::vector<real_t> x,
          const std::vector<real_t> y,
          const std::vector<real_t> z) {
  assert((x.size() == y.size()) && (y.size() == z.size()));

  real_t* node_coord = n->Coord();
  int depth          = n->Depth();
  real_t length      = static_cast<real_t>(std::pow(0.5, depth));

  real_t xmin = node_coord[0];
  real_t xmax = xmin + length;
  real_t ymin = node_coord[1];
  real_t ymax = ymin + length;
  real_t zmin = node_coord[2];
  real_t zmax = zmin + length;

  std::vector<int> out_index_list;
  for (int i = 0; i < x.size(); i++) {
    if ( x[i] < xmin || x[i] > xmax) {
      out_index_list.push_back(i);
      continue;
    }
    if ( y[i] < ymin || y[i] > ymax) {
      out_index_list.push_back(i);
      continue;
    }
    if ( z[i] < zmin || z[i] > zmax) {
      out_index_list.push_back(i);
      continue;
    }
  }
  return out_index_list;
}

}
#endif // SRC_TREE_TREE_UTILS_H_
