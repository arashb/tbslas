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
#include <string>

#include "utils/profile.h"
#include "utils/common.h"

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

  max_value = lcl_max_values[0];
  max_depth = static_cast<int>(lcl_max_values[1]);
}

template<typename TreeType>
int
tree2vec(TreeType& tree,
         std::vector<typename TreeType::Real_t>& vec) {
  typedef typename TreeType::Real_t RealType;
  typedef typename TreeType::Node_t NodeType;

  RealType scal_exp = 1.0;
  int sdim        = 3;
  std::vector<NodeType*> nlist;

  // Get non-ghost, leaf nodes.
  std::vector<NodeType*>& nlist_ = tree.GetNodeList();
  for (size_t i = 0; i < nlist_.size(); i++) {
    if (nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()) {
      nlist.push_back(nlist_[i]);
    }
  }
  assert(nlist.size()>0);

  // int omp_p       = omp_get_max_threads();
  int cheb_deg = nlist[0]->ChebDeg();
  int data_dof = nlist[0]->DataDOF();
  size_t n_coeff3 = (cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;

  // TODO: add OpenMP parallelization
  //#pragma omp parallel for
  //for (size_t tid=0;tid<omp_p;tid++) {
  // size_t i_start=(nlist.size()* tid   )/omp_p;
  // size_t i_end  =(nlist.size()*(tid+1))/omp_p;
  // vec.reserve(n_coeff3*data_dof*nlist.size());

  size_t i_start = 0;
  size_t i_end   = nlist.size();
  for (size_t i = i_start; i < i_end; i++) {
    pvfmm::Vector<RealType>& coeff_vec = nlist[i]->ChebData();
    double s = 1;//std::pow(0.5,sdim*nlist[i]->Depth()*0.5*scal_exp);

    size_t vec_offset = i*n_coeff3*data_dof;
    for (size_t j = 0; j < n_coeff3*data_dof; j++)
      // vec[j+vec_offset] = coeff_vec[j]*s;
      vec.push_back(coeff_vec[j]*s);
  }
  return 0;
}

// int vec2tree(Vec& Y, FMMData fmm_data){
//   PetscErrorCode ierr;
//   FMM_Tree_t* tree=fmm_data.tree;
//   const MPI_Comm* comm=tree->Comm();
//   int cheb_deg=fmm_data.fmm_mat->ChebDeg();

//   std::vector<FMMNode_t*> nlist;
//   { // Get non-ghost, leaf nodes.
//     std::vector<FMMNode_t*>& nlist_=tree->GetNodeList();
//     for(size_t i=0;i<nlist_.size();i++){
//       if(nlist_[i]->IsLeaf() && !nlist_[i]->IsGhost()){
//         nlist.push_back(nlist_[i]);
//       }
//     }
//   }

//   int omp_p=omp_get_max_threads();
//   size_t n_coeff3=(cheb_deg+1)*(cheb_deg+2)*(cheb_deg+3)/6;

//   {
//     PetscInt Y_size;
//     ierr = VecGetLocalSize(Y, &Y_size);
//     int data_dof=Y_size/(n_coeff3*nlist.size());

//     const PetscScalar *Y_ptr;
//     ierr = VecGetArrayRead(Y, &Y_ptr);

// #pragma omp parallel for
//     for(size_t tid=0;tid<omp_p;tid++){
//       size_t i_start=(nlist.size()* tid   )/omp_p;
//       size_t i_end  =(nlist.size()*(tid+1))/omp_p;
//       for(size_t i=i_start;i<i_end;i++){
//         pvfmm::Vector<double>& coeff_vec=nlist[i]->ChebData();
//         double s=std::pow(2.0,COORD_DIM*nlist[i]->Depth()*0.5*SCAL_EXP);

//         size_t Y_offset=i*n_coeff3*data_dof;
//         for(size_t j=0;j<n_coeff3*data_dof;j++) coeff_vec[j]=PetscRealPart(Y_ptr[j+Y_offset])*s;
//         nlist[i]->DataDOF()=data_dof;
//       }
//     }
//   }

//   return 0;
// }

template<typename TreeType>
void
save_tree(TreeType& tree, std::string file_name) {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

  std::vector<RealType> vec;
  tbslas::tree2vec<TreeType>(tree, vec);
  for (int i = 0; i < vec.size(); i++) {
    printf("%f ", vec[i]);
  }

}

template<typename TreeType>
void
load_tree() {
  typedef typename TreeType::Node_t NodeType;
  typedef typename TreeType::Real_t RealType;

}

}  // namespace tbslas
#endif  // SRC_TREE_TREE_UTILS_H_
