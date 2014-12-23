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
        const size_t n_pts=part_indx[j+1]-part_indx[j];
        if(!n_pts) continue;

        Real_t* c=nodes[j]->Coord();
        size_t d=nodes[j]->Depth();
        Real_t s=(Real_t)(1ULL<<d);

        if (tmp_out.Dim()<n_pts*data_dof) {
          tmp_out.ReInit(n_pts*data_dof);
        }
        tmp_out.SetZero();

        Real_t* coord_ptr=&trg_coord[0]+part_indx[j]*COORD_DIM;

        //////////////////////////////////////////////////////////////
        // CHEBYSHEV INTERPOLATION
        //////////////////////////////////////////////////////////////
        // coord.resize(n_pts*COORD_DIM);

        // for (size_t i=0;i<n_pts;i++) {
        //   // scale to [-1,1] -> used in cheb_eval
        //   coord[i*COORD_DIM+0]=(coord_ptr[i*COORD_DIM+0]-c[0])*2.0*s-1.0;
        //   coord[i*COORD_DIM+1]=(coord_ptr[i*COORD_DIM+1]-c[1])*2.0*s-1.0;
        //   coord[i*COORD_DIM+2]=(coord_ptr[i*COORD_DIM+2]-c[2])*2.0*s-1.0;
        // }

        // Profile<double>::Tic("ChebEval", false, 5);
        // pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
        // pvfmm::cheb_eval(coeff, nodes[j]->ChebDeg(), coord, tmp_out);
        // Profile<double>::Toc();

        //////////////////////////////////////////////////////////////
        // CUBIC INTERPOLATION
        //////////////////////////////////////////////////////////////
        {
        // ************************************************************
        // CONSTRUCT REGULAR GRID
        // ************************************************************
        Profile<double>::Tic("ConstructRegGrid", false, 5);
        int reg_grid_resolution = nodes[j]->ChebDeg()*4;
        Real_t spacing = 1.0/(reg_grid_resolution-1);
        std::vector<Real_t> reg_grid_coord_1d(reg_grid_resolution);
        tbslas::get_reg_grid_points<Real_t, 1>(reg_grid_resolution,
                                               reg_grid_coord_1d.data());
        Profile<double>::Toc();

        // ************************************************************
        // EVALUATE AT THE REGULAR GRID
        // ************************************************************
        Profile<double>::Tic("EvalRegGrid", false, 5);
        int reg_grid_num_points = std::pow(reg_grid_resolution, COORD_DIM);
        std::vector<Real_t> reg_grid_vals(reg_grid_num_points*data_dof);
        // ************************************************************
        // CHEB_EVAL
        // ************************************************************
        // scale to [-1,1] -> used in cheb_eval
        // std::vector<Real_t> x(reg_grid_resolution);
        // std::vector<Real_t> y(reg_grid_resolution);
        // std::vector<Real_t> z(reg_grid_resolution);
        // for(size_t i=0;i<reg_grid_resolution;i++) {
        //   x[i] = -1.0+2.0*reg_grid_coord_1d[i];
        //   y[i] = -1.0+2.0*reg_grid_coord_1d[i];
        //   z[i] = -1.0+2.0*reg_grid_coord_1d[i];
        // }
        // pvfmm::Vector<Real_t> reg_grid_vals_tmp(reg_grid_num_points*data_dof);
        // pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
        // pvfmm::cheb_eval(coeff, nodes[j]->ChebDeg(), x, y, z, reg_grid_vals_tmp);

        // for(size_t i = 0; i < reg_grid_num_points;i++) {
        //   for (int data_dof_cnt = 0; data_dof_cnt < data_dof; data_dof_cnt++) {
        //     reg_grid_vals[i+data_dof_cnt*reg_grid_num_points] =
        //         reg_grid_vals_tmp[i*data_dof+data_dof_cnt];
        //   }
        // }
        // ************************************************************
        // READ_EVAL
        // ************************************************************
        {
          Real_t length = 1.0/s;
          std::vector<Real_t> x(reg_grid_resolution);
          std::vector<Real_t> y(reg_grid_resolution);
          std::vector<Real_t> z(reg_grid_resolution);
          // scale to global coordinates -> used in node->ReadVal
          for(size_t i = 0; i < reg_grid_resolution; i++) {
            x[i] = c[0] + length * reg_grid_coord_1d[i];
            y[i] = c[1] + length * reg_grid_coord_1d[i];
            z[i] = c[2] + length * reg_grid_coord_1d[i];
          }
          nodes[j]->ReadVal(x, y, z,
                            reg_grid_vals.data());
        }

        Profile<double>::Toc();

        // ************************************************************
        // TRICUBIC INTERPOLATION
        // ************************************************************
        // Profile<double>::Tic("InterpTriCubic", false, 5);
        // for (int data_dim_cnt = 0; data_dim_cnt < data_dof; data_dim_cnt++) {
        //   likely::TriCubicInterpolator tc_inerpolator(&reg_grid_vals[data_dim_cnt*reg_grid_num_points],
        //                                               spacing,
        //                                               reg_grid_resolution);
        //   for (int pi = 0; pi < n_pts; pi++) {
        //     Real_t xq =(coord_ptr[pi*COORD_DIM+0]-c[0])*s;
        //     Real_t yq =(coord_ptr[pi*COORD_DIM+1]-c[1])*s;
        //     Real_t zq =(coord_ptr[pi*COORD_DIM+2]-c[2])*s;
        //     tmp_out[pi*data_dof + data_dim_cnt] = tc_inerpolator(xq, yq, zq);
        //   }
        // }
        // Profile<double>::Toc();

        // ************************************************************
        // 3D CUBIC INTERPOLATION
        // ************************************************************
        Profile<double>::Tic("InterpCubic", false, 5);
        for(int data_dim_cnt = 0; data_dim_cnt < data_dof; data_dim_cnt++) {
          Real_t* reg_grid_vals_dim =
              &reg_grid_vals[data_dim_cnt*reg_grid_num_points];
          for ( int pi = 0; pi < n_pts; pi++) {

            //in case that target coordinate is outside the unit domain->set to 0
            if ( coord_ptr[pi*COORD_DIM+0] < 0 || coord_ptr[pi*COORD_DIM+0] > 1.0 ||
                 coord_ptr[pi*COORD_DIM+1] < 0 || coord_ptr[pi*COORD_DIM+1] > 1.0 ||
                 coord_ptr[pi*COORD_DIM+2] < 0 || coord_ptr[pi*COORD_DIM+2] > 1.0
                 ) {
              tmp_out[pi*data_dof + data_dim_cnt] = 0;
              continue;
            }

            // scale to [0,1] in local node
            Real_t xq =(coord_ptr[pi*COORD_DIM+0]-c[0])*s;
            Real_t yq =(coord_ptr[pi*COORD_DIM+1]-c[1])*s;
            Real_t zq =(coord_ptr[pi*COORD_DIM+2]-c[2])*s;

            Real_t dx(xq/spacing);
            Real_t dy(yq/spacing);
            Real_t dz(zq/spacing);

            // Calculate the corresponding lower-bound grid indices.
            int xi = static_cast<int>(std::floor(dx));
            int yi = static_cast<int>(std::floor(dy));
            int zi = static_cast<int>(std::floor(dz));
            // printf("TAR. POINTS: [%f, %f, %f] -> INDEX:[%d, %d, %d]\n",
            //        xq, yq, zq, xi, yi, zi);
            assert(xi >= 0 && yi >= 0  && zi >= 0);
            assert(xi < reg_grid_resolution &&
                   yi < reg_grid_resolution &&
                   zi < reg_grid_resolution);

            int xshift, yshift, zshift;
            if (xi == 0) xshift = 0;
            else if (xi == reg_grid_resolution -2)
              xshift = xi -2;
            else xshift = xi - 1;

            if (yi == 0) yshift = 0;
            else if (yi == reg_grid_resolution -2)
              yshift = yi -2;
            else yshift = yi - 1;

            if (zi == 0) zshift = 0;
            else if (zi == reg_grid_resolution -2)
              zshift = zi -2;
            else zshift = zi - 1;

            Real_t xx[4];
            Real_t yy[4];
            Real_t zz[4];
            for (int i = 0; i < 4; i++) {
              xx[i] = reg_grid_coord_1d[xshift + i];
              yy[i] = reg_grid_coord_1d[yshift + i];
              zz[i] = reg_grid_coord_1d[zshift + i];
            }

            int xli, yli, zli;
            int index;
            Real_t pp[4][4][4];
            for(int i = 0; i < 4; i++) {
              xli = xshift + i;
              for ( int j = 0; j < 4; j++) {
                yli = yshift + j;
                for ( int k = 0; k < 4; k++) {
                  zli = zshift + k;
                  index = xli + yli*reg_grid_resolution + zli*reg_grid_resolution*reg_grid_resolution;
                  pp[i][j][k] = reg_grid_vals_dim[index];
                }
              }
            }
            tmp_out[pi*data_dof + data_dim_cnt] =
                tbslas::CubicInterpPolicy<Real_t>::InterpCubic3D(xq, yq, zq,
                                                                 xx, yy, zz,
                                                                 pp);
          }  // for target point
        }  // for data_dof
        Profile<double>::Toc();
        }
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
    Profile<double>::Tic("EvalTree", false, 0);
    EvalTree(node_, num_points, const_cast<real_t*>(points_pos), out);
    Profile<double>::Toc();
  }

 private:
  NodeType* node_;
};

}      // namespace tbslas
#endif  // SRC_TREE_NODE_FIELD_FUNCTOR_H_
