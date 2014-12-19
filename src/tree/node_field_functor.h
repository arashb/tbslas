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

#ifndef SRC_TREE_NODE_FIELD_FUNCTOR_H_
#define SRC_TREE_NODE_FIELD_FUNCTOR_H_

#include <vector>

#include <pvfmm_common.hpp>
#include <cheb_node.hpp>

#include "utils/profile.h"
#include "semilag/cubic_interp_policy.h"

namespace tbslas {

template <class Tree_t>
void EvalTree(Tree_t* tree,
          size_t N,
          typename Tree_t::Real_t* trg_coord_,
          typename Tree_t::Real_t* value) {

  typedef typename Tree_t::Node_t Node_t;
  typedef typename Tree_t::Real_t Real_t;

  size_t data_dof=0;
  pvfmm::MortonId min_mid;
  std::vector<Node_t*> nodes;
  { // Get list of leaf nodes.
    const std::vector<Node_t*>& all_nodes=tree->GetNodeList();
    for (size_t i=0; i< all_nodes.size(); i++) {
      if (all_nodes[i]->IsLeaf() && !all_nodes[i]->IsGhost()) {
        nodes.push_back(all_nodes[i]);
      }
    }
    assert(nodes.size());
    min_mid=nodes[0]->GetMortonId();
    data_dof=nodes[0]->DataDOF();
  }

  pvfmm::Vector<pvfmm::MortonId> trg_mid(N);
  #pragma omp parallel for
  for (size_t i = 0; i < N; i++) {
    trg_mid[i] = pvfmm::MortonId(&trg_coord_[i*COORD_DIM]);
  }

  // Compute scatter_index.
  pvfmm::Vector<size_t> scatter_index;
  pvfmm::par::SortScatterIndex(trg_mid  , scatter_index, *tree->Comm(), &min_mid);

  // Scatter coordinates and values.
  pvfmm::Vector<Real_t> trg_coord(N*COORD_DIM,trg_coord_);
  pvfmm::par::ScatterForward  (trg_mid  , scatter_index, *tree->Comm());
  pvfmm::par::ScatterForward  (trg_coord, scatter_index, *tree->Comm());

  pvfmm::Vector<Real_t> trg_value(trg_mid.Dim()*data_dof);
  { // Read tree data
    std::vector<size_t> part_indx(nodes.size()+1);
    part_indx[nodes.size()]=trg_mid.Dim();
    #pragma omp parallel for
    for (size_t j=0;j<nodes.size();j++) {
      part_indx[j]=std::lower_bound(&trg_mid[0],
                                    &trg_mid[0]+trg_mid.Dim(),
                                    nodes[j]->GetMortonId()) - &trg_mid[0];
    }

    size_t omp_p=omp_get_max_threads();
    #pragma omp parallel for
    for (size_t pid=0;pid<omp_p;pid++) {
      size_t a=((pid+0)*nodes.size())/omp_p;
      size_t b=((pid+1)*nodes.size())/omp_p;

      std::vector<Real_t> coord;
      pvfmm::Vector<Real_t> tmp_out;
      for (size_t j=a;j<b;j++) {
        size_t n_pts=part_indx[j+1]-part_indx[j];
        if(!n_pts) continue;

        Real_t* c=nodes[j]->Coord();
        size_t d=nodes[j]->Depth();
        Real_t s=(Real_t)(1ULL<<d);

        coord.resize(n_pts*COORD_DIM);
        Real_t* coord_ptr=&trg_coord[0]+part_indx[j]*COORD_DIM;
        for (size_t i=0;i<n_pts;i++) {
          coord[i*COORD_DIM+0]=(coord_ptr[i*COORD_DIM+0]-c[0])*s;
          coord[i*COORD_DIM+1]=(coord_ptr[i*COORD_DIM+1]-c[1])*s;
          coord[i*COORD_DIM+2]=(coord_ptr[i*COORD_DIM+2]-c[2])*s;
        }

        if (tmp_out.Dim()<n_pts*data_dof) {
          tmp_out.ReInit(n_pts*data_dof);
        }
        tmp_out.SetZero();
        // ************************************************************
        // CUBIC INTERPOLATION
        // ************************************************************
        // construct regular grid
        Profile<double>::Tic("ConstructRegGrid", false, 5);
        int reg_grid_resolution = nodes[j]->ChebDeg()*8;
        std::vector<Real_t> reg_grid_coord_1d(reg_grid_resolution);
        tbslas::get_reg_grid_points<Real_t, 1>(reg_grid_resolution,
                                               reg_grid_coord_1d.data());

        Real_t length = 1.0/s;
        std::vector<Real_t> x(reg_grid_resolution);
        std::vector<Real_t> y(reg_grid_resolution);
        std::vector<Real_t> z(reg_grid_resolution);
        for(size_t i=0;i<reg_grid_resolution;i++) {
          x[i]=nodes[j]->Coord()[0]+length*reg_grid_coord_1d[i];
          y[i]=nodes[j]->Coord()[1]+length*reg_grid_coord_1d[i];
          z[i]=nodes[j]->Coord()[2]+length*reg_grid_coord_1d[i];
        }
        Profile<double>::Toc();
        // for (int i = 0; i< reg_grid_resolution; i++)
        //   printf("point: [%f]\n",
        //          reg_grid_coord_1d[i]);
        // evaluate values at the regular grid points
        Profile<double>::Tic("ChebEvalRegGrid", false, 5);
        int reg_grid_num_points = std::pow(reg_grid_resolution, COORD_DIM);
        std::vector<Real_t> reg_grid_vals(reg_grid_num_points*data_dof);
        // void cheb_eval(const Vector<T>& coeff_, int cheb_deg, const std::vector<T>& in_x, const std::vector<T>& in_y, const std::vector<T>& in_z, Vector<T>& out, mem::MemoryManager* mem_mgr){

        nodes[j]->ReadVal(x, y, z, reg_grid_vals.data());
        Profile<double>::Toc();
        // cubic interpolatation at target points
        // pvfmm::Vector<Real_t> reg_grid_vals_tmp(reg_grid_num_points*data_dof);
        // Profile<double>::Tic("ConvertMemLayout", false, 1);
        // for (int i = 0; i < reg_grid_num_points; i++) {
        //   for (int data_dof_cnt = 0; data_dof_cnt < data_dof; data_dof_cnt++) {
        //     reg_grid_vals_tmp[i+reg_grid_num_points*data_dof_cnt] = reg_grid_vals[i*data_dof+data_dof_cnt];
        //   }
        // }
        // Profile<double>::Toc();
        Profile<double>::Tic("TriCubicInterp", false, 5);
        Real_t spacing = 1.0/(reg_grid_resolution-1);
        for (int data_dim_cnt = 0; data_dim_cnt < data_dof; data_dim_cnt++) {
          likely::TriCubicInterpolator tc_inerpolator(&reg_grid_vals[data_dim_cnt*reg_grid_num_points],
                                                      spacing,
                                                      reg_grid_resolution);
          for (int i = 0; i < n_pts; i++) {
            tmp_out[i*data_dof+data_dim_cnt] = tc_inerpolator(coord[i*COORD_DIM+0],
                                                              coord[i*COORD_DIM+1],
                                                              coord[i*COORD_DIM+2]);
          }
        }
        Profile<double>::Toc();

        // ************************************************************
        // CHEBYSHEV INTERPOLATION
        // ************************************************************
        // Profile<double>::Tic("ChebEval", false, 5);
        // pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
        // pvfmm::cheb_eval(coeff, nodes[j]->ChebDeg(), coord, tmp_out);
        // Profile<double>::Toc();
        memcpy(&trg_value[0]+part_indx[j]*data_dof, &tmp_out[0], n_pts*data_dof*sizeof(Real_t));
      }
    }
  }

  pvfmm::par::ScatterReverse(trg_value, scatter_index, *tree->Comm(), N);
  memcpy(value, &trg_value[0], trg_value.Dim()*sizeof(Real_t));
}

template<typename real_t,
         class NodeType>
class NodeFieldFunctor {

 public:
  explicit NodeFieldFunctor(NodeType* node): node_(node) {
  }

  virtual ~NodeFieldFunctor() {
  }

  void operator () (const real_t* points_pos,
                    int num_points,
                    real_t* out) {
    Profile<double>::Tic("EvalTree",false, 0);
    EvalTree(node_, num_points, const_cast<real_t*>(points_pos), out);
    Profile<double>::Toc();
  }

 private:
  NodeType* node_;
};

}      // namespace tbslas
#endif  // SRC_TREE_NODE_FIELD_FUNCTOR_H_
