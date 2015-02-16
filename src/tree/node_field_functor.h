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
#include <profile.hpp>

#include "utils/common.h"
#include "semilag/cubic_interp_policy.h"

namespace tbslas {

template <class Tree_t>
void EvalTree(Tree_t* tree,
              size_t N,
              typename Tree_t::Real_t* trg_coord_,
              typename Tree_t::Real_t* value,
              pvfmm::BoundaryType bc_type) {
  typedef typename Tree_t::Node_t Node_t;
  typedef typename Tree_t::Real_t Real_t;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

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

  if (bc_type == pvfmm::Periodic) {
    #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
      Real_t* c = &trg_coord_[i*COORD_DIM];
      if (c[0] < 0)
        c[0] = c[0] + 1.0;
      else if ( c[0] > 1.0)
        c[0] = c[0] - 1.0;

      if (c[1] < 0)
        c[1] = c[1] + 1.0;
      else if ( c[1] > 1.0)
        c[1] = c[1] - 1.0;

      if (c[2] < 0)
        c[2] = c[2] + 1.0;
      else if ( c[2] > 1.0)
        c[2] = c[2] - 1.0;
    }
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
      pvfmm::Vector<Real_t> tmp_out;    // buffer used in chebyshev evaluation
      std::vector<Real_t> query_values; // buffer used in cubic interpolation
      Real_t* output = NULL;

      for (size_t j=a;j<b;j++) {
        const size_t n_pts=part_indx[j+1]-part_indx[j];
        if(!n_pts) continue;

        Real_t* c=nodes[j]->Coord();
        size_t d=nodes[j]->Depth();
        Real_t s=(Real_t)(1ULL<<d);

        Real_t* coord_ptr=&trg_coord[0]+part_indx[j]*COORD_DIM;

        if (!sim_config->use_cubic) {
          //////////////////////////////////////////////////////////////
          // CHEBYSHEV INTERPOLATION
          //////////////////////////////////////////////////////////////
          if (tmp_out.Dim()<n_pts*data_dof) {
            tmp_out.ReInit(n_pts*data_dof);
          }
          tmp_out.SetZero();
          coord.resize(n_pts*COORD_DIM);
          for (size_t i=0;i<n_pts;i++) {
            // scale to [-1,1] -> used in cheb_eval
            coord[i*COORD_DIM+0]=(coord_ptr[i*COORD_DIM+0]-c[0])*2.0*s-1.0;
            coord[i*COORD_DIM+1]=(coord_ptr[i*COORD_DIM+1]-c[1])*2.0*s-1.0;
            coord[i*COORD_DIM+2]=(coord_ptr[i*COORD_DIM+2]-c[2])*2.0*s-1.0;
          }

          pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
          pvfmm::cheb_eval(coeff, nodes[j]->ChebDeg(), coord, tmp_out);
          output = &tmp_out[0];
        } else {
          //////////////////////////////////////////////////////////////
          // CUBIC INTERPOLATION
          //////////////////////////////////////////////////////////////
          query_values.resize(n_pts*data_dof);
          // ************************************************************
          // CONSTRUCT REGULAR GRID
          // ************************************************************
          int reg_grid_resolution =
              nodes[j]->ChebDeg()*sim_config->cubic_upsampling_factor;
          Real_t spacing = 1.0/(reg_grid_resolution-1);
          std::vector<Real_t> reg_grid_coord_1d(reg_grid_resolution);
          tbslas::get_reg_grid_points<Real_t, 1>(reg_grid_resolution,
                                                 reg_grid_coord_1d.data());
          // ************************************************************
          // EVALUATE AT THE REGULAR GRID
          // ************************************************************
          int reg_grid_num_points = std::pow(reg_grid_resolution, COORD_DIM);
          std::vector<Real_t> reg_grid_vals(reg_grid_num_points*data_dof);
          pvfmm::Vector<Real_t> reg_grid_vals_tmp(reg_grid_num_points*data_dof);
          // evaluate using the tree
          if(!sim_config->cubic_use_analytical) {
            // scale to [-1,1] -> used in cheb_eval
            std::vector<Real_t> x(reg_grid_resolution);
            for(size_t i=0;i<reg_grid_resolution;i++) {
              x[i] = -1.0+2.0*reg_grid_coord_1d[i];
            }
            pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
            pvfmm::cheb_eval(coeff, nodes[j]->ChebDeg(), x, x, x, reg_grid_vals_tmp);
          } else {   // evaluate using analytical function
            std::vector<Real_t> reg_grid_anal_coord(3);
            std::vector<Real_t> reg_grid_anal_vals(1*data_dof);
            int nx = reg_grid_resolution;
            for (int xi = 0; xi < nx; xi++) {
              for (int yi = 0; yi < nx; yi++) {
                for (int zi = 0; zi < nx; zi++) {
                  reg_grid_anal_coord[0] = c[0] + reg_grid_coord_1d[xi]/s;
                  reg_grid_anal_coord[1] = c[1] + reg_grid_coord_1d[yi]/s;
                  reg_grid_anal_coord[2] = c[2] + reg_grid_coord_1d[zi]/s;
                  assert(!nodes[j]->input_fn.IsEmpty());
                  nodes[j]->input_fn(reg_grid_anal_coord.data(),
                                     1,
                                     reg_grid_anal_vals.data());
                  for(int l=0;l<data_dof;l++)
                    reg_grid_vals_tmp[xi+(yi+(zi+l*nx)*nx)*nx] = reg_grid_anal_vals[l];
                }
              }
            }
          }
          for (int i = 0; i < reg_grid_vals.size(); i++) {
            reg_grid_vals[i] = reg_grid_vals_tmp[i];
          }
          // ************************************************************
          // 3D CUBIC INTERPOLATION
          // ************************************************************
          // scale to [0,1] in local node
          std::vector<Real_t> query_points(n_pts*COORD_DIM);
          for ( int pi = 0; pi < n_pts; pi++) {
            query_points[pi*COORD_DIM+0] = (coord_ptr[pi*COORD_DIM+0]-c[0])*s;
            query_points[pi*COORD_DIM+1] = (coord_ptr[pi*COORD_DIM+1]-c[1])*s;
            query_points[pi*COORD_DIM+2] = (coord_ptr[pi*COORD_DIM+2]-c[2])*s;
          }
          tbslas::CubicInterpPolicy<Real_t>::interp(reg_grid_vals,
                                                    data_dof,
                                                    spacing,
                                                    query_points,
                                                    query_values);
          output = &query_values[0];
        } // end of cubic interpolation

        memcpy(&trg_value[0]+part_indx[j]*data_dof, output, n_pts*data_dof*sizeof(Real_t));
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
  explicit NodeFieldFunctor(NodeType* node):
      node_(node) {
  }

  virtual ~NodeFieldFunctor() {
  }

  void operator () (const real_t* points_pos,
                    int num_points,
                    real_t* out) {
    tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Tic("EvalTree", &sim_config->comm, false, 5);
    EvalTree(node_, num_points, const_cast<real_t*>(points_pos), out,sim_config->bc);
    pvfmm::Profile::Toc();
  }

 private:
  NodeType* node_;
};

}      // namespace tbslas
#endif  // SRC_TREE_NODE_FIELD_FUNCTOR_H_
