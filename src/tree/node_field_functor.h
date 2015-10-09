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

#include <ompUtils.h>

#include "utils/common.h"

namespace tbslas {


template <class Real_t>
void fast_interp(const std::vector<Real_t>& reg_grid_vals, int data_dof,
              int N_reg, const std::vector<Real_t>& query_points,
              std::vector<Real_t>& query_values){

  Real_t lagr_denom[4];
  for(int i=0;i<4;i++){
    lagr_denom[i]=1;
    for(int j=0;j<4;j++){
      if(i!=j) lagr_denom[i]/=(Real_t)(i-j);
    }
  }

  int N_reg3=N_reg*N_reg*N_reg;
  int N_pts=query_points.size()/COORD_DIM;
  query_values.resize(N_pts*data_dof);

  for(int i=0;i<N_pts;i++){
    if( query_points[COORD_DIM*i+0] < 0 || query_points[COORD_DIM*i+0] > 1.0 ||
        query_points[COORD_DIM*i+1] < 0 || query_points[COORD_DIM*i+1] > 1.0 ||
        query_points[COORD_DIM*i+2] < 0 || query_points[COORD_DIM*i+2] > 1.0 ){
      for(int k=0;k<data_dof;k++) query_values[i*data_dof+k]=0;
      continue;
    }

    Real_t point[COORD_DIM];
    int grid_indx[COORD_DIM];
    for(int j=0;j<COORD_DIM;j++){
      point[j]=query_points[COORD_DIM*i+j]*(N_reg-1);
      grid_indx[j]=((int)point[j])-1;
      if(grid_indx[j]<      0) grid_indx[j]=      0;
      if(grid_indx[j]>N_reg-4) grid_indx[j]=N_reg-4;
      point[j]-=grid_indx[j];
    }

    Real_t M[3][4];
    for(int j=0;j<COORD_DIM;j++){
      Real_t x=point[j];
      for(int k=0;k<4;k++){
        M[j][k]=lagr_denom[k];
        for(int l=0;l<4;l++){
          if(k!=l) M[j][k]*=(x-l);
        }
      }
    }

    for(int k=0;k<data_dof;k++){
      Real_t val=0;
      for(int j2=0;j2<4;j2++){
        for(int j1=0;j1<4;j1++){
          for(int j0=0;j0<4;j0++){
            int indx = (grid_indx[0]+j0) + N_reg*(grid_indx[1]+j1) + N_reg*N_reg*(grid_indx[2]+j2);
            val += M[0][j0]*M[1][j1]*M[2][j2] * reg_grid_vals[indx+k*N_reg3];
          }
        }
      }
      query_values[i*data_dof+k]=val;
    }
  }
}


template <class Tree_t>
void EvalTree(Tree_t* tree,
              size_t N,
              typename Tree_t::Real_t* trg_coord_,
              typename Tree_t::Real_t* value,
              pvfmm::BoundaryType bc_type) {
  size_t omp_p=omp_get_max_threads();
  typedef typename Tree_t::Node_t Node_t;
  typedef typename Tree_t::Real_t Real_t;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  int myrank;
  int np;
  MPI_Comm_rank(sim_config->comm, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  //////////////////////////////////////////////////////////////
  // GET LEAF NODES AND MINIMUM MORTON ID OF THE CURRENT PROCESS
  //////////////////////////////////////////////////////////////
  pvfmm::Profile::Tic("MortonId", &sim_config->comm, false, 5);
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
  //std::cout<<"Leaf Nodes = "<<nodes.size()<<'\n';

  //////////////////////////////////////////////////////////////
  // APPLY PERIODIC BOUNDARY CONDITION
  //////////////////////////////////////////////////////////////
  if (bc_type == pvfmm::Periodic) {
    #pragma omp parallel for
    for (size_t i = 0; i < N*COORD_DIM; i++) {
      Real_t& c = trg_coord_[i];
      if(c <  0.0) c = c + 1.0;
      if(c >= 1.0) c = c - 1.0;
    }
  }

  //////////////////////////////////////////////////////////////
  // COMPUTE MORTON ID OF THE TARGET POINTS
  //////////////////////////////////////////////////////////////
  static pvfmm::Vector<pvfmm::MortonId> trg_mid; trg_mid.Resize(N);
  #pragma omp parallel for
  for (size_t i = 0; i < N; i++) {
    trg_mid[i] = pvfmm::MortonId(&trg_coord_[i*COORD_DIM]);
  }
  pvfmm::Profile::Toc();

  //////////////////////////////////////////////////////////////
  // GATHER ALL THE MINIMUM MORTION IDS
  //////////////////////////////////////////////////////////////
  // TODO: do not need allgather
  std::vector<pvfmm::MortonId> glb_min_mid(np);
  MPI_Allgather(&min_mid, 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                &glb_min_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), sim_config->comm);

  //////////////////////////////////////////////////////////////
  // LOCAL SORT
  //////////////////////////////////////////////////////////////
  // size_t lcl_start, lcl_end, lcl_trg_cnt;

  // {
  //   // LOCAL SORT WITHOUT TRACKING THE INDICES
  //   pvfmm::omp_par::merge_sort(&trg_mid[0], &trg_mid[0]+trg_mid.Dim());
  //   lcl_start = std::lower_bound(&trg_mid[0],
  //                                &trg_mid[0]+trg_mid.Dim(),
  //                                glb_min_mid[myrank]) - &trg_mid[0];
  //   if (myrank+1 < np)
  //     lcl_end = std::lower_bound(&trg_mid[0],
  //                                &trg_mid[0]+trg_mid.Dim(),
  //                                glb_min_mid[myrank+1]) - &trg_mid[0];
  //   else
  //     lcl_end = trg_mid.Dim();
  //   lcl_trg_cnt = lcl_end - lcl_start;
  // }

//   {
//     // LOCAL SORT WITH TRACKING THE INDICES
//     typedef pvfmm::par::SortPair<pvfmm::MortonId,size_t> Pair_t;
//     pvfmm::Vector<Pair_t> parray(trg_mid.Dim());
//     {
//       long long loc_size=trg_mid.Dim();
// #pragma omp parallel for
//       for(size_t i=0; i < loc_size; i++) {
//       parray[i].key  = trg_mid[i];
//       parray[i].data = i;
//     }
//     }

//     pvfmm::Vector<Pair_t> parray_sorted;
//     pvfmm::par::HyperQuickSort(parray, parray_sorted, MPI_COMM_SELF);

//     Pair_t p1; p1.key = glb_min_mid[myrank];
//     lcl_start = std::lower_bound(&parray_sorted[0],
//           &parray_sorted[0]+parray_sorted.Dim(),
//           p1,
//           std::less<Pair_t>()) - &parray_sorted[0];

//     if (myrank+1 < np) {
//       Pair_t p2; p2.key = glb_min_mid[myrank+1];
//       lcl_end = std::lower_bound(&parray_sorted[0],
//           &parray_sorted[0]+parray_sorted.Dim(),
//           p2,
//           std::less<Pair_t>()) - &parray_sorted[0];
//     } else {
//       lcl_end = parray_sorted.Dim();
//     }

//     static pvfmm::Vector<size_t> scatter_index;
//     scatter_index.Resize(parray_sorted.Dim());
// #pragma omp parallel for
//     for(size_t i=0; i<parray_sorted.Dim(); i++) {
//       scatter_index[i] = parray_sorted[i].data;
//     }

//     lcl_trg_cnt = lcl_end - lcl_start;
//   }

//   /* print number of departure points in current process */
//   int total_trg_cnt = trg_mid.Dim();
//   int* rbuf_total_trg_cnt = (int *)malloc(np*sizeof(int));
//   int* rbuf_lcl_trg_cnt = (int *)malloc(np*sizeof(int));
//   MPI_Gather(&total_trg_cnt, 1, MPI_INT, rbuf_total_trg_cnt, 1, MPI_INT, 0, *tree->Comm());
//   MPI_Gather(&lcl_trg_cnt, 1, MPI_INT, rbuf_lcl_trg_cnt, 1, MPI_INT, 0, *tree->Comm());
//   if (!myrank) {
//     std::cout << "LCL_TRG_CNT/TOT_TRG_CNT: ";
//     for (int i = 0 ; i < np; i++)
//       std::cout << " " << rbuf_lcl_trg_cnt[i]<<"/" << rbuf_total_trg_cnt[i];
//     std::cout << std::endl;
//   }
//   delete rbuf_total_trg_cnt;
//   delete rbuf_lcl_trg_cnt;

  //////////////////////////////////////////////////////////////
  // LOCAL EVALUATION
  //////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////
  // GLOBAL EVALUATION
  //////////////////////////////////////////////////////////////
  pvfmm::Profile::Tic("ScatterIndex", &sim_config->comm, true, 5);
  static pvfmm::Vector<size_t> scatter_index;
  pvfmm::par::SortScatterIndex(trg_mid  , scatter_index, *tree->Comm(), &min_mid);
  pvfmm::Profile::Toc();

  // Scatter coordinates and values.
  static pvfmm::Vector<Real_t> trg_coord;
  pvfmm::Profile::Tic("ScatterForward", &sim_config->comm, true, 5);
  {
    trg_coord.Resize(N*COORD_DIM);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=N*COORD_DIM*(tid+0)/omp_p;
      size_t b=N*COORD_DIM*(tid+1)/omp_p;
      if(b-a) memcpy(&trg_coord[0]+a, &trg_coord_[0]+a, (b-a)*sizeof(Real_t));
    }
    pvfmm::par::ScatterForward  (trg_coord, scatter_index, *tree->Comm());

    trg_mid.Resize(trg_coord.Dim()/COORD_DIM);
    #pragma omp parallel for
    for(size_t i=0;i<trg_mid.Dim();i++){
      trg_mid[i] = pvfmm::MortonId(&trg_coord[i*COORD_DIM]);
    }
  }
  pvfmm::Profile::Toc();

  /* print number of departure points in current process */
  int trg_cnt = trg_mid.Dim();
  int* rbuf_trg_cnt = (int *)malloc(np*sizeof(int));
  MPI_Gather(&trg_cnt, 1, MPI_INT, rbuf_trg_cnt, 1, MPI_INT, 0, *tree->Comm());
  if (!myrank) {
    std::cout << "TRG_CNT: ";
    for (int i = 0 ; i < np; i++)
      std::cout << " " << rbuf_trg_cnt[i];
    std::cout << std::endl;
  }
  delete rbuf_trg_cnt;

  static pvfmm::Vector<Real_t> trg_value;
  pvfmm::Profile::Tic("Evaluation", &sim_config->comm, false, 5);
  trg_value.Resize(trg_mid.Dim()*data_dof);
  { // Read tree data
    std::vector<size_t> part_indx(nodes.size()+1);
    part_indx[nodes.size()] = trg_mid.Dim();
    #pragma omp parallel for
    for (size_t j=0;j<nodes.size();j++) {
      part_indx[j]=std::lower_bound(&trg_mid[0],
                                    &trg_mid[0]+trg_mid.Dim(),
                                    nodes[j]->GetMortonId()) - &trg_mid[0];
    }

    #pragma omp parallel for
    for (size_t pid=0;pid<omp_p;pid++) {
      size_t a=((pid+0)*nodes.size())/omp_p;
      size_t b=((pid+1)*nodes.size())/omp_p;

      std::vector<Real_t> coord;
      pvfmm::Vector<Real_t> tmp_out;    // buffer used in chebyshev evaluation
      std::vector<Real_t> query_values; // buffer used in cubic interpolation
      std::vector<Real_t> query_points;
      Real_t* output = NULL;

      // ************************************************************
      // CONSTRUCT REGULAR GRID
      // ************************************************************
      int reg_grid_resolution =
          nodes[0]->ChebDeg()*sim_config->cubic_upsampling_factor;
      Real_t spacing = 1.0/(reg_grid_resolution-1);
      std::vector<Real_t> reg_grid_coord_1d(reg_grid_resolution);
      tbslas::get_reg_grid_points<Real_t, 1>(reg_grid_resolution,
                                             reg_grid_coord_1d.data());

      // ************************************************************
      // EVALUATE AT THE REGULAR GRID
      // ************************************************************
      int reg_grid_num_points = std::pow(reg_grid_resolution, COORD_DIM);
      std::vector<Real_t> reg_grid_vals(reg_grid_num_points*data_dof);

      // scale to [-1,1] -> used in cheb_eval
      std::vector<Real_t> x(reg_grid_resolution);
      for(size_t i=0;i<reg_grid_resolution;i++) {
        x[i] = -1.0+2.0*reg_grid_coord_1d[i];
      }

      pvfmm::Matrix<Real_t> Mp1;
      pvfmm::Vector<Real_t> v1, v2;
      { // Precomputation
        int cheb_deg=nodes[0]->ChebDeg();
        std::vector<Real_t> p1(reg_grid_resolution*(cheb_deg+1));
        pvfmm::cheb_poly(cheb_deg,&x[0],reg_grid_resolution,&p1[0]);
        Mp1.ReInit(cheb_deg+1,reg_grid_resolution,&p1[0]);

        // Create work buffers
        size_t buff_size=std::max(cheb_deg+1,reg_grid_resolution)*std::max(cheb_deg+1,reg_grid_resolution)*std::max(cheb_deg+1,reg_grid_resolution)*data_dof;
        v1.Resize(buff_size);
        v2.Resize(buff_size);
      }

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
            tmp_out.Resize(n_pts*data_dof);
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
          query_points.resize(n_pts*COORD_DIM);
          // ************************************************************
          // EVALUATE AT THE REGULAR GRID
          // ************************************************************
          if(!sim_config->cubic_use_analytical) {
            pvfmm::Vector<Real_t>& coeff_=nodes[j]->ChebData();
            pvfmm::Vector<Real_t> reg_grid_vals_tmp(reg_grid_num_points*data_dof, &reg_grid_vals[0], false);
            { // cheb_eval
              int cheb_deg=nodes[0]->ChebDeg();
              size_t d=(size_t)cheb_deg+1;
              size_t n_coeff=(d*(d+1)*(d+2))/6;
              size_t dof=coeff_.Dim()/n_coeff;
              assert(coeff_.Dim()==dof*n_coeff);

              size_t n1=x.size();
              assert(reg_grid_vals_tmp.Dim()==n1*n1*n1*dof);

              { // Rearrange coefficients into a tensor.
                pvfmm::Vector<Real_t> coeff(d*d*d*dof,&v1[0],false);
                coeff.SetZero();
                size_t indx=0;
                for(size_t l=0;l<dof;l++){
                  for(size_t i=0;i<d;i++){
                    for(size_t j=0;j<d-i;j++){
                      Real_t* coeff_ptr=&coeff[(j+(i+l*d)*d)*d];
                      for(size_t k=0;k<d-i-j;k++){
                        coeff_ptr[k]=coeff_[indx];
                        indx++;
                      }
                    }
                  }
                }
              }

              { // Apply Mp1
                pvfmm::Matrix<Real_t> Mi  ( d* d*dof, d,&v1[0],false);
                pvfmm::Matrix<Real_t> Mo  ( d* d*dof,n1,&v2[0],false);
                pvfmm::Matrix<Real_t>::GEMM(Mo, Mi, Mp1);

                pvfmm::Matrix<Real_t> Mo_t(n1, d* d*dof,&v1[0],false);
                for(size_t i=0;i<Mo.Dim(0);i++)
                for(size_t j=0;j<Mo.Dim(1);j++){
                  Mo_t[j][i]=Mo[i][j];
                }
              }
              { // Apply Mp1
                pvfmm::Matrix<Real_t> Mi  (n1* d*dof, d,&v1[0],false);
                pvfmm::Matrix<Real_t> Mo  (n1* d*dof,n1,&v2[0],false);
                pvfmm::Matrix<Real_t>::GEMM(Mo, Mi, Mp1);

                pvfmm::Matrix<Real_t> Mo_t(n1,n1* d*dof,&v1[0],false);
                for(size_t i=0;i<Mo.Dim(0);i++)
                for(size_t j=0;j<Mo.Dim(1);j++){
                  Mo_t[j][i]=Mo[i][j];
                }
              }
              { // Apply Mp1
                pvfmm::Matrix<Real_t> Mi  (n1*n1*dof, d,&v1[0],false);
                pvfmm::Matrix<Real_t> Mo  (n1*n1*dof,n1,&v2[0],false);
                pvfmm::Matrix<Real_t>::GEMM(Mo, Mi, Mp1);

                pvfmm::Matrix<Real_t> Mo_t(n1,n1*n1*dof,&v1[0],false);
                for(size_t i=0;i<Mo.Dim(0);i++)
                for(size_t j=0;j<Mo.Dim(1);j++){
                  Mo_t[j][i]=Mo[i][j];
                }
              }

              { // Copy to reg_grid_vals_tmp
                pvfmm::Matrix<Real_t> Mo  ( n1*n1*n1,dof,&v1[0],false);
                pvfmm::Matrix<Real_t> Mo_t(dof,n1*n1*n1,&reg_grid_vals_tmp[0],false);
                for(size_t i=0;i<Mo.Dim(0);i++)
                for(size_t j=0;j<Mo.Dim(1);j++){
                  Mo_t[j][i]=Mo[i][j];
                }
              }
            }
          } else {   // evaluate using analytical function
            std::vector<Real_t> reg_grid_anal_coord(3);
            std::vector<Real_t> reg_grid_anal_vals(1*data_dof);
            pvfmm::Vector<Real_t> reg_grid_vals_tmp(reg_grid_num_points*data_dof, &reg_grid_vals[0], false);
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
          // ************************************************************
          // 3D CUBIC INTERPOLATION
          // ************************************************************
          // scale to [0,1] in local node
          for ( int pi = 0; pi < n_pts; pi++) {
            query_points[pi*COORD_DIM+0] = (coord_ptr[pi*COORD_DIM+0]-c[0])*s;
            query_points[pi*COORD_DIM+1] = (coord_ptr[pi*COORD_DIM+1]-c[1])*s;
            query_points[pi*COORD_DIM+2] = (coord_ptr[pi*COORD_DIM+2]-c[2])*s;
          }
          fast_interp(reg_grid_vals, data_dof, reg_grid_resolution, query_points, query_values);
          output = &query_values[0];
        } // end of cubic interpolation

        memcpy(&trg_value[0]+part_indx[j]*data_dof, output, n_pts*data_dof*sizeof(Real_t));
      }
    }
    pvfmm::Profile::Add_FLOP(trg_coord.Dim()/COORD_DIM * (COORD_DIM*16 + data_dof*256)); // cubic interpolation
  }
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("ScatterReverse", &sim_config->comm, true, 5);
  pvfmm::par::ScatterReverse(trg_value, scatter_index, *tree->Comm(), N);
  memcpy(value, &trg_value[0], trg_value.Dim()*sizeof(Real_t));
  pvfmm::Profile::Toc();
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
    pvfmm::Profile::Tic("EvalTree", &sim_config->comm, true, 5);
    EvalTree(node_, num_points, const_cast<real_t*>(points_pos), out,sim_config->bc);
    pvfmm::Profile::Toc();
  }

 private:
  NodeType* node_;
};

}      // namespace tbslas
#endif  // SRC_TREE_NODE_FIELD_FUNCTOR_H_
