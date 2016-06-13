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
#include <cmath>
#include <pvfmm_common.hpp>
#include <cheb_node.hpp>
#include <profile.hpp>

#include <ompUtils.h>

#include "utils/common.h"

namespace pvfmm {
template <class Real, class Vec = Real>
void vec_eval(int n, int d0, int d, int dof, Real *px, Real *py, Real *pz,
              Real *coeff, Real *tmp_out) {
  const int VecLen = sizeof(Vec) / sizeof(Real);

  for (int l0 = 0; l0 < dof; l0++) {
    for (int l1 = 0; l1 < n; l1 += 4*VecLen) {
      Vec u0 = zero_intrin<Vec>();
      Vec u1 = zero_intrin<Vec>();
      Vec u2 = zero_intrin<Vec>();
      Vec u3 = zero_intrin<Vec>();
      for (int i = 0; i < d; i++) {
        Vec v0 = zero_intrin<Vec>();
        Vec v1 = zero_intrin<Vec>();
        Vec v2 = zero_intrin<Vec>();
        Vec v3 = zero_intrin<Vec>();
        Vec pz0_ = load_intrin<Vec>(&pz[i * n + l1+0*VecLen]);
        Vec pz1_ = load_intrin<Vec>(&pz[i * n + l1+1*VecLen]);
        Vec pz2_ = load_intrin<Vec>(&pz[i * n + l1+2*VecLen]);
        Vec pz3_ = load_intrin<Vec>(&pz[i * n + l1+3*VecLen]);
        for (int j = 0; i + j < d; j++) {
          Vec w0 = zero_intrin<Vec>();
          Vec w1 = zero_intrin<Vec>();
          Vec w2 = zero_intrin<Vec>();
          Vec w3 = zero_intrin<Vec>();
          Vec py0_ = load_intrin<Vec>(&py[j * n + l1+0*VecLen]);
          Vec py1_ = load_intrin<Vec>(&py[j * n + l1+1*VecLen]);
          Vec py2_ = load_intrin<Vec>(&py[j * n + l1+2*VecLen]);
          Vec py3_ = load_intrin<Vec>(&py[j * n + l1+3*VecLen]);
          for (int k = 0; i + j + k < d; k++) {
            Vec px0_ = load_intrin<Vec>(&px[k * n + l1+0*VecLen]);
            Vec px1_ = load_intrin<Vec>(&px[k * n + l1+1*VecLen]);
            Vec px2_ = load_intrin<Vec>(&px[k * n + l1+2*VecLen]);
            Vec px3_ = load_intrin<Vec>(&px[k * n + l1+3*VecLen]);
            Vec c = set_intrin<Vec, Real>(coeff[k + d0 * (j + d0 * (i + d0 * l0))]);
            w0 = add_intrin(w0, mul_intrin(px0_, c));
            w1 = add_intrin(w1, mul_intrin(px1_, c));
            w2 = add_intrin(w2, mul_intrin(px2_, c));
            w3 = add_intrin(w3, mul_intrin(px3_, c));
          }
          v0 = add_intrin(v0, mul_intrin(py0_, w0));
          v1 = add_intrin(v1, mul_intrin(py1_, w1));
          v2 = add_intrin(v2, mul_intrin(py2_, w2));
          v3 = add_intrin(v3, mul_intrin(py3_, w3));
        }
        u0 = add_intrin(u0, mul_intrin(pz0_, v0));
        u1 = add_intrin(u1, mul_intrin(pz1_, v1));
        u2 = add_intrin(u2, mul_intrin(pz2_, v2));
        u3 = add_intrin(u3, mul_intrin(pz3_, v3));
      }
      store_intrin(&tmp_out[l0 * n + l1+0*VecLen], u0);
      store_intrin(&tmp_out[l0 * n + l1+1*VecLen], u1);
      store_intrin(&tmp_out[l0 * n + l1+2*VecLen], u2);
      store_intrin(&tmp_out[l0 * n + l1+3*VecLen], u3);
    }
  }
}
}

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
          Real_t M1M2=M[1][j1]*M[2][j2];
          long indx_j1j2 = N_reg*( (grid_indx[1]+j1) + N_reg*(grid_indx[2]+j2) );
          for(int j0=0;j0<4;j0++){
            long indx = (grid_indx[0]+j0) + indx_j1j2;
            val += M[0][j0]*M1M2 * reg_grid_vals[indx+k*N_reg3];
          }
        }
      }
      query_values[i*data_dof+k]=val;
    }
  }
}

template <class Real_t, class Tree_t>
void EvalNodesLocal(std::vector<typename Tree_t::Node_t*>& nodes,
                    pvfmm::Vector<Real_t>& trg_coord,
                    pvfmm::Vector<Real_t>& trg_value) { // Read nodes data
  size_t omp_p=omp_get_max_threads();
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  size_t data_dof=nodes[0]->DataDOF();

  //////////////////////////////////////////////////////////////
  // COMPUTE MORTON IDs OF TARGET POINTS
  //////////////////////////////////////////////////////////////
  static pvfmm::Vector<pvfmm::MortonId> trg_mid;
  trg_mid.Resize(trg_coord.Dim()/COORD_DIM);

#pragma omp parallel for
  for(size_t i=0; i < trg_mid.Dim(); i++) {
    // trg_mid[i] = pvfmm::MortonId(&trg_coord[i*COORD_DIM]);
    Real_t shift = 1.0/(1UL<<MAX_DEPTH);
    trg_mid[i] = pvfmm::MortonId(( (trg_coord[i*COORD_DIM+0] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord[i*COORD_DIM+0]-shift:trg_coord[i*COORD_DIM+0],
                                 ( (trg_coord[i*COORD_DIM+1] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord[i*COORD_DIM+1]-shift:trg_coord[i*COORD_DIM+1],
                                 ( (trg_coord[i*COORD_DIM+2] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord[i*COORD_DIM+2]-shift:trg_coord[i*COORD_DIM+2]);

  }

  //////////////////////////////////////////////////////////////
  // DETERMINE THE CORRESPONDING NODE FOR EACH TARGET POINT
  //////////////////////////////////////////////////////////////
  std::vector<size_t> part_indx(nodes.size()+1);
  part_indx[nodes.size()] = trg_mid.Dim();

#pragma omp parallel for
  for (size_t j=0;j<nodes.size();j++) {
    part_indx[j]=std::lower_bound(&trg_mid[0],
                                  &trg_mid[0]+trg_mid.Dim(),
                                  nodes[j]->GetMortonId()) - &trg_mid[0];
  }

  // std::cout << "PART_INDX: " ;
  // for (int i = 0; i < part_indx.size(); i++) {
  //   std::cout << " " << part_indx[i];
  // }
  // std::cout << std::endl;

  // std::cout << "PART_INDX_LB: " ;
  // for (int i = 0; i < part_indx.size()-1; i++) {
  //   std::cout << " " << part_indx[i+1] - part_indx[i];
  // }
  // std::cout << std::endl;

#pragma omp parallel for schedule(static)
  for (size_t pid=0;pid<omp_p;pid++) {
    std::vector<Real_t> coord;
    pvfmm::Vector<Real_t> tmp_out;    // buffer used in chebyshev evaluation

    std::vector<Real_t> cx;
    std::vector<Real_t> cy;
    std::vector<Real_t> cz;
    pvfmm::Matrix<Real_t> px;
    pvfmm::Matrix<Real_t> py;
    pvfmm::Matrix<Real_t> pz;

    pvfmm::Vector<Real_t> coeff_;
    { // Init coeff_
      int d0=nodes[0]->ChebDeg()+1+8;
      coeff_.ReInit(d0*d0*d0*data_dof);
      coeff_.SetZero();
    }

    size_t pt_start0=((pid+0)*trg_mid.Dim())/omp_p;
    size_t pt_end0=((pid+1)*trg_mid.Dim())/omp_p;
    long node_start=std::lower_bound(&part_indx[0], &part_indx[0]+part_indx.size(), pt_start0)-&part_indx[0]-1;
    long node_end=std::lower_bound(&part_indx[0], &part_indx[0]+part_indx.size(), pt_end0)-&part_indx[0];
    node_start=std::max<long>(0,std::min<long>(node_start, nodes.size()));
    node_end=std::max<long>(0,std::min<long>(node_end, nodes.size()));
    for (size_t j=node_start;j<node_end;j++) {
      size_t pt_start_=std::max<long>(pt_start0, part_indx[j]);
      size_t pt_end_=std::min<long>(pt_end0, part_indx[j+1]);
      if(pt_start_>=pt_end_) continue;

      { // set coeff_
        pvfmm::Vector<Real_t>& coeff=nodes[j]->ChebData();
        int cheb_deg=nodes[j]->ChebDeg();
        int d=cheb_deg+1;
        int d0;
        for(d0=d;d0%4;d0++);

        assert(coeff.Dim()==(size_t)(d*(d+1)*(d+2)*data_dof)/6);
        long indx=0;
        for(int l0=0;l0<data_dof;l0++){
          for(int i=0;i<d;i++){
            for(int j=0;i+j<d;j++){
              for(int k=0;i+j+k<d;k++){
                coeff_[k+d0*(j+d0*(i+d0*l0))]=coeff[indx];
                indx++;
              }
            }
          }
        }
      }

      const size_t blk_size=288;
      for(size_t pt_start=pt_start_;pt_start<pt_end_;pt_start+=blk_size){
        size_t pt_end=std::min(pt_end_, pt_start+blk_size);
        const size_t n_pts=pt_end-pt_start;

        Real_t* coord_ptr = &trg_coord[0]+pt_start*COORD_DIM;
        {
          //////////////////////////////////////////////////////////////
          // CHEBYSHEV INTERPOLATION
          //////////////////////////////////////////////////////////////
          coord.resize((n_pts+16)*COORD_DIM);
          { // Set coord
            Real_t* c = nodes[j]->Coord();
            size_t  d = nodes[j]->Depth();
            Real_t  s = (Real_t)(1ULL<<d);
            for (size_t i=0;i<n_pts;i++) {
              // scale to [-1,1] -> used in cheb_eval
              coord[i*COORD_DIM+0]=(coord_ptr[i*COORD_DIM+0]-c[0])*2.0*s-1.0;
              coord[i*COORD_DIM+1]=(coord_ptr[i*COORD_DIM+1]-c[1])*2.0*s-1.0;
              coord[i*COORD_DIM+2]=(coord_ptr[i*COORD_DIM+2]-c[2])*2.0*s-1.0;
            }
          }
          for (size_t i=n_pts;i<n_pts+16;i++) {
            // scale to [-1,1] -> used in cheb_eval
            coord[i*COORD_DIM+0]=0.0;
            coord[i*COORD_DIM+1]=0.0;
            coord[i*COORD_DIM+2]=0.0;
          }

          if(coord.size()){
            int cheb_deg=nodes[j]->ChebDeg();
            int d=cheb_deg+1;
            int d0;
            for(d0=d;d0%4;d0++);
            int n=n_pts;
            for(;n%16;n++);

            cx.resize(n);
            cy.resize(n);
            cz.resize(n);
            for(long i=0;i<n;i++){
              cx[i]=coord[i*COORD_DIM+0];
              cy[i]=coord[i*COORD_DIM+1];
              cz[i]=coord[i*COORD_DIM+2];
            }

            px.Resize(d0,n); px.SetZero();
            py.Resize(d0,n); py.SetZero();
            pz.Resize(d0,n); pz.SetZero();
            pvfmm::cheb_poly(cheb_deg,&(cx[0]),n,&(px[0][0]));
            pvfmm::cheb_poly(cheb_deg,&(cy[0]),n,&(py[0][0]));
            pvfmm::cheb_poly(cheb_deg,&(cz[0]),n,&(pz[0][0]));

            if (tmp_out.Dim()<n*data_dof) {
              tmp_out.Resize(n*data_dof);
            }
            if(pvfmm::mem::TypeTraits<Real_t>::ID()==pvfmm::mem::TypeTraits<float>::ID()){
              typedef float Real;
              #if defined __MIC__
                #define Vec_t Real_t
              #elif defined __AVX__
                #define Vec_t __m256
              #elif defined __SSE3__
                #define Vec_t __m128
              #else
                #define Vec_t Real
              #endif
              ::pvfmm::vec_eval<Real, Vec_t>(n,d0,d,data_dof, (Real*)&px[0][0], (Real*)&py[0][0], (Real*)&pz[0][0], (Real*)&coeff_[0], (Real*)&tmp_out[0]);
              #undef Vec_t
            }else if(pvfmm::mem::TypeTraits<Real_t>::ID()==pvfmm::mem::TypeTraits<double>::ID()){
              typedef double Real;
              #if defined __MIC__
                #define Vec_t Real_t
              #elif defined __AVX__
                #define Vec_t __m256d
              #elif defined __SSE3__
                #define Vec_t __m128d
              #else
                #define Vec_t Real
              #endif
              ::pvfmm::vec_eval<Real, Vec_t>(n,d0,d,data_dof, (Real*)&px[0][0], (Real*)&py[0][0], (Real*)&pz[0][0], (Real*)&coeff_[0], (Real*)&tmp_out[0]);
              #undef Vec_t
            }else{
              typedef Real_t Real;
              #define Vec_t Real
              ::pvfmm::vec_eval<Real, Vec_t>(n,d0,d,data_dof, (Real*)&px[0][0], (Real*)&py[0][0], (Real*)&pz[0][0], (Real*)&coeff_[0], (Real*)&tmp_out[0]);
              #undef Vec_t
            }

            { // Copy to trg_value
              Real_t *trg_value_ = &trg_value[0] + pt_start * data_dof;
              for (long i = 0; i < n_pts; i++) {
                for (int j = 0; j < data_dof; j++) {
                  trg_value_[i * data_dof + j] = tmp_out[j * n + i];
                }
              }
            }
          }
        }
      }
    }
  }
  { // Add FLOP
    long d=nodes[0]->ChebDeg()+1;
    pvfmm::Profile::Add_FLOP(trg_coord.Dim()/COORD_DIM * ( COORD_DIM*d*3 + ((d*(d+1)*(d+2))/6) * data_dof * 2 ) );
  }
}

template <class Tree_t>
void EvalTree(Tree_t* tree,
              typename Tree_t::Real_t* trg_coord_,
              size_t N,
              typename Tree_t::Real_t* value,
              pvfmm::BoundaryType bc_type) {
  size_t omp_p=omp_get_max_threads();
  typedef typename Tree_t::Node_t Node_t;
  typedef typename Tree_t::Real_t Real_t;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  int myrank;
  int np;
  MPI_Comm_rank(sim_config->comm, &myrank);
  MPI_Comm_size(sim_config->comm, &np);

  //////////////////////////////////////////////////////////////
  // GET LEAF NODES AND MINIMUM MORTON ID OF THE CURRENT PROCESS
  //////////////////////////////////////////////////////////////
  // pvfmm::Profile::Tic("MinMortonId", &sim_config->comm, false, 5);
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
  // pvfmm::Profile::Toc();

  //////////////////////////////////////////////////////////////
  // GATHER MINIMUM MORTON IDS OF ALL PARTITIONS
  //////////////////////////////////////////////////////////////
  std::vector<pvfmm::MortonId> glb_min_mid(np);
  MPI_Allgather(&min_mid, 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                &glb_min_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                sim_config->comm);

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
  // LOCAL SORT
  //////////////////////////////////////////////////////////////
  typedef pvfmm::par::SortPair<pvfmm::MortonId,size_t> Pair_t;
  //pvfmm::Vector<Pair_t> iarray_trg_mid(N);
  pvfmm::Vector<Pair_t> iarray_trg_mid_sorted(N);
  size_t lcl_start, lcl_end, trg_cnt_inside, trg_cnt_outside;

  // pvfmm::Profile::Tic("LclSort", &sim_config->comm, false, 5);
  //////////////////////////////////////////////////
  // LOCAL SORT WITH TRACKING THE INDICES
  //////////////////////////////////////////////////
  pvfmm::Profile::Tic("LclHQSort", &sim_config->comm, false, 5);
  Real_t shift = 1.0/(1UL<<MAX_DEPTH);
#pragma omp parallel for
  for(size_t i = 0; i < N; i++) {
    iarray_trg_mid_sorted[i].key  = pvfmm::MortonId(( (trg_coord_[i*COORD_DIM+0] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_[i*COORD_DIM+0]-shift:trg_coord_[i*COORD_DIM+0],
                                                    ( (trg_coord_[i*COORD_DIM+1] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_[i*COORD_DIM+1]-shift:trg_coord_[i*COORD_DIM+1],
                                                    ( (trg_coord_[i*COORD_DIM+2] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_[i*COORD_DIM+2]-shift:trg_coord_[i*COORD_DIM+2]);
    // iarray_trg_mid_sorted[i].key  = pvfmm::MortonId(&trg_coord_[i*COORD_DIM]);
    iarray_trg_mid_sorted[i].data = i;
  }

  //pvfmm::par::HyperQuickSort(iarray_trg_mid, iarray_trg_mid_sorted, MPI_COMM_SELF);
  pvfmm::omp_par::merge_sort(&iarray_trg_mid_sorted[0], &iarray_trg_mid_sorted[0]+iarray_trg_mid_sorted.Dim());

  Pair_t p1;
  p1.key = glb_min_mid[myrank];
  lcl_start = std::lower_bound(&iarray_trg_mid_sorted[0],
                               &iarray_trg_mid_sorted[0]+iarray_trg_mid_sorted.Dim(),
                               p1,
                               std::less<Pair_t>()) - &iarray_trg_mid_sorted[0];

  if (myrank+1 < np) {
    Pair_t p2; p2.key = glb_min_mid[myrank+1];
    lcl_end = std::lower_bound(&iarray_trg_mid_sorted[0],
                               &iarray_trg_mid_sorted[0]+iarray_trg_mid_sorted.Dim(),
                               p2,
                               std::less<Pair_t>()) - &iarray_trg_mid_sorted[0];
  } else {
    lcl_end = iarray_trg_mid_sorted.Dim();
  }

  // [lcl_start, lcl_end[
  trg_cnt_inside  = lcl_end - lcl_start;
  trg_cnt_outside = N - trg_cnt_inside;
  pvfmm::Profile::Toc();  // Sort

  //////////////////////////////////////////////////
  // COMMINUCATE THE OUTSIDER POINTS
  //////////////////////////////////////////////////
  static pvfmm::Vector<size_t> out_scatter_index;
  static pvfmm::Vector<Real_t> trg_coord_outside;

  // pvfmm::Profile::Tic("OutCpyCoord", &sim_config->comm, true, 5);
  trg_coord_outside.Resize(trg_cnt_outside*COORD_DIM);
#pragma omp parallel for
  for (int i = 0 ; i < lcl_start; i++) {
    size_t src_indx = iarray_trg_mid_sorted[i].data;
    size_t dst_indx = i;
    for(size_t j = 0; j < COORD_DIM; j++) {
      trg_coord_outside[dst_indx*COORD_DIM+j] = trg_coord_[src_indx*COORD_DIM+j];
    }
  }

#pragma omp parallel for
  for (int i = 0 ; i < (N-lcl_end); i++) {
    size_t src_indx = iarray_trg_mid_sorted[lcl_end+i].data;
    size_t dst_indx = lcl_start + i;
    for(size_t j = 0; j < COORD_DIM; j++) {
      trg_coord_outside[dst_indx*COORD_DIM+j] = trg_coord_[src_indx*COORD_DIM+j];
    }
  }
  // pvfmm::Profile::Toc();  // OUT COPY COORDINATES

  // pvfmm::Profile::Tic("OutMortonID", &sim_config->comm, true, 5);
  static pvfmm::Vector<pvfmm::MortonId> trg_mid_outside;
  trg_mid_outside.Resize(trg_cnt_outside);
#pragma omp parallel for
  for (size_t i = 0; i < trg_cnt_outside; i++) {
    trg_mid_outside[i] = pvfmm::MortonId(( (trg_coord_outside[i*COORD_DIM+0] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_outside[i*COORD_DIM+0]-shift:trg_coord_outside[i*COORD_DIM+0],
                                         ( (trg_coord_outside[i*COORD_DIM+1] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_outside[i*COORD_DIM+1]-shift:trg_coord_outside[i*COORD_DIM+1],
                                         ( (trg_coord_outside[i*COORD_DIM+2] == 1.0) && (sim_config->bc != pvfmm::Periodic) )?trg_coord_outside[i*COORD_DIM+2]-shift:trg_coord_outside[i*COORD_DIM+2]);

    // trg_mid_outside[i] = pvfmm::MortonId(&trg_coord_outside[i*COORD_DIM]);
  }
  // pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("OutScatterIndex", &sim_config->comm, true, 5);
  pvfmm::par::SortScatterIndex(trg_mid_outside, out_scatter_index, *tree->Comm(), &min_mid);
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("OutScatterForward", &sim_config->comm, true, 5);
  pvfmm::par::ScatterForward(trg_coord_outside, out_scatter_index, *tree->Comm());
  pvfmm::Profile::Toc();

  size_t trg_cnt_others = trg_coord_outside.Dim()/COORD_DIM;
  size_t trg_cnt_total  = trg_cnt_inside + trg_cnt_others;

  //////////////////////////////////////////////////
  // EVALUATE THE OUTSIDER POINTS
  //////////////////////////////////////////////////
  static pvfmm::Vector<Real_t> trg_value_outsider;
  pvfmm::Profile::Tic("OutEvaluation", &sim_config->comm, false, 5);
  trg_value_outsider.Resize(trg_cnt_others*data_dof);
  EvalNodesLocal<Real_t, Tree_t>(nodes, trg_coord_outside, trg_value_outsider);
  pvfmm::Profile::Toc();

  //////////////////////////////////////////////////
  // SCATTER REVERSE THE OUTSIDER POINT'S VALUES
  //////////////////////////////////////////////////
  pvfmm::Profile::Tic("OutScatterReverse", &sim_config->comm, true, 5);
  pvfmm::par::ScatterReverse(trg_value_outsider, out_scatter_index, *tree->Comm(), trg_cnt_outside);
  pvfmm::Profile::Toc();

  //////////////////////////////////////////////////
  // SET OUTSIDER POINTS EVALUATION VALUES
  //////////////////////////////////////////////////
  // pvfmm::Profile::Tic("OutSetVal", &sim_config->comm, true, 5);
#pragma omp parallel for
  for (int i = 0 ; i < lcl_start; i++) {
    size_t src_indx = i;
    size_t dst_indx = iarray_trg_mid_sorted[i].data;
    for(size_t j = 0; j < data_dof; j++) {
      value[dst_indx*data_dof+j] = trg_value_outsider[src_indx*data_dof+j];
    }
  }
#pragma omp parallel for
  for (int i = 0 ; i < (N-lcl_end); i++) {
    size_t src_indx = lcl_start + i;
    size_t dst_indx = iarray_trg_mid_sorted[lcl_end+i].data;
    for(size_t j = 0; j < data_dof; j++) {
      value[dst_indx*data_dof+j] = trg_value_outsider[src_indx*data_dof+j];
    }
  }
  // pvfmm::Profile::Toc();  // OUT SET VALUES

  //////////////////////////////////////////////////
  // PRINT TARGET COUNTS INFO
  //////////////////////////////////////////////////
  /* if (sim_config->profile) { */
  /*   size_t sbuff[4] = {trg_cnt_inside, */
  /*                      trg_cnt_outside, */
  /*                      trg_value_outsider.Dim()/data_dof, */
  /*                      trg_cnt_others}; */
  /*   size_t* rbuff = (size_t *)malloc(np*4*sizeof(size_t)); */
  /*   MPI_Gather(sbuff, 4, pvfmm::par::Mpi_datatype<size_t>::value(), */
  /*              rbuff, 4, pvfmm::par::Mpi_datatype<size_t>::value(), */
  /*              0, *tree->Comm()); */
  /*   if (myrank == 0) { */
  /*     std::ostringstream os; */
  /*     os << "TRG_CNT_IN_TOT: "; */
  /*     for (int i = 0 ; i < np; i++) { */
  /*       size_t* data = &rbuff[i*4]; */
  /*       std::cout */
  /*           << " PROC: " << i */
  /*           << " TRG_CNT_IN: "     << data[0] */
  /*           << " TRG_CNT_OUT: "    << data[1] */
  /*           << " TRG_CNT_RCV: "    << data[2] */
  /*           << " TRG_CNT_OTHRS: "  << data[3] */
  /*           << std::endl; */
  /*       os << data[0] + data[3] << " "; */
  /*     } */
  /*     std::cout << os.str() << std::endl; */
  /*   } */
  /*   delete rbuff; */
  /* } */

  //////////////////////////////////////////////////
  // COLLECT THE COORDINATE VALUES
  //////////////////////////////////////////////////
  static  pvfmm::Vector<Real_t> trg_coord_inside;
  // pvfmm::Profile::Tic("InCpyCoord", &sim_config->comm, true, 5);
  trg_coord_inside.Resize(trg_cnt_inside*COORD_DIM);
#pragma omp parallel for
  for (size_t i = 0; i < trg_cnt_inside; i++) {
    size_t src_indx = iarray_trg_mid_sorted[lcl_start+i].data;
    size_t dst_indx = i;
    for(size_t j = 0; j < COORD_DIM;j++) {
      trg_coord_inside[dst_indx*COORD_DIM+j] = trg_coord_[src_indx*COORD_DIM+j];
    }
  }
  // pvfmm::Profile::Toc();

  //////////////////////////////////////////////////
  // EVALUATE THE LOCAL POINTS
  //////////////////////////////////////////////////
  static pvfmm::Vector<Real_t> trg_value_insider;
  pvfmm::Profile::Tic("InEvaluation", &sim_config->comm, false, 5);
  trg_value_insider.Resize(trg_cnt_inside*data_dof);
  EvalNodesLocal<Real_t, Tree_t>(nodes, trg_coord_inside, trg_value_insider);
  pvfmm::Profile::Toc();

  //////////////////////////////////////////////////
  // SET INSIDER POINTS EVALUATION VALUES
  //////////////////////////////////////////////////
  // pvfmm::Profile::Tic("InSetVal", &sim_config->comm, false, 5);
#pragma omp parallel for
  for (size_t i = 0; i < trg_cnt_inside; i++) {
    size_t src_indx = i;
    size_t dst_indx = iarray_trg_mid_sorted[lcl_start+i].data;
    for (int j = 0; j < data_dof; j++) {
      value[dst_indx*data_dof+j] = trg_value_insider[i*data_dof+j];
    }
  }
  // pvfmm::Profile::Toc();

  // pvfmm::Profile::Toc();  // LOCAL SORT

  //////////////////////////////////////////////////
  // PRINT TARGET COUNTS INFO
  //////////////////////////////////////////////////
  // if (sim_config->profile) {
  //   size_t sbuff[4] = {trg_cnt_inside,
  //                      trg_cnt_outside,
  //                      trg_value_outsider.Dim()/data_dof,
  //                      trg_cnt_others};
  //   size_t* rbuff = (size_t *)malloc(np*4*sizeof(size_t));
  //   MPI_Gather(sbuff, 4, pvfmm::par::Mpi_datatype<size_t>::value(),
  //              rbuff, 4, pvfmm::par::Mpi_datatype<size_t>::value(),
  //              0, *tree->Comm());
  //   if (myrank == 0) {
  //     std::ostringstream os;
  //     os << "TRG_CNT_IN_TOT: ";
  //     for (int i = 0 ; i < np; i++) {
  //       size_t* data = &rbuff[i*4];
  //       std::cout
  //           << " PROC: " << i
  //           << " TRG_CNT_IN: "     << data[0]
  //           << " TRG_CNT_OUT: "    << data[1]
  //           << " TRG_CNT_OTHRS_RCV: "  << data[2]
  //           << " TRG_CNT_OTHRS_SND: "  << data[3]
  //           << std::endl;
  //       os << data[0] + data[3] << " ";
  //     }
  //     std::cout << os.str() << std::endl;
  //   }
  //   delete rbuff;
  // }

  //////////////////////////////////////////////////////////////
  // GLOBAL SORT
  //////////////////////////////////////////////////////////////
//   pvfmm::Profile::Tic("GlobalSort", &sim_config->comm, false, 5);
//   {
//     //////////////////////////////////////////////////////////////
//     // COMPUTE MORTON ID OF THE TARGET POINTS
//     //////////////////////////////////////////////////////////////
//     pvfmm::Profile::Tic("MortonId", &sim_config->comm, true, 5);
//     static pvfmm::Vector<pvfmm::MortonId> trg_mid; trg_mid.Resize(N);
// #pragma omp parallel for
//     for (size_t i = 0; i < N; i++) {
//       trg_mid[i] = pvfmm::MortonId(&trg_coord_[i*COORD_DIM]);
//     }
//     pvfmm::Profile::Toc();

//     //////////////////////////////////////////////////////////////
//     // SCATTER THE COORDINATES
//     //////////////////////////////////////////////////////////////
//     pvfmm::Profile::Tic("ScatterIndex", &sim_config->comm, true, 5);
//     static pvfmm::Vector<size_t> scatter_index;
//     pvfmm::par::SortScatterIndex(trg_mid, scatter_index, *tree->Comm(), &min_mid);
//     pvfmm::Profile::Toc();

//     static pvfmm::Vector<Real_t> trg_coord;
//     pvfmm::Profile::Tic("ScatterForward", &sim_config->comm, true, 5);
//     {
//       trg_coord.Resize(N*COORD_DIM);
// #pragma omp parallel for
//       for(size_t tid=0;tid<omp_p;tid++){
//         size_t a=N*COORD_DIM*(tid+0)/omp_p;
//         size_t b=N*COORD_DIM*(tid+1)/omp_p;
//         if(b-a) memcpy(&trg_coord[0]+a, &trg_coord_[0]+a, (b-a)*sizeof(Real_t));
//       }
//       pvfmm::par::ScatterForward(trg_coord, scatter_index, *tree->Comm());
//     }
//     pvfmm::Profile::Toc();

//     // std::cout << "P" << myrank << " TRG_CNT: " << trg_coord.Dim()/COORD_DIM << std::endl;

//     //////////////////////////////////////////////////////////////
//     // LOCAL POINTS EVALUATION
//     //////////////////////////////////////////////////////////////
//     size_t num_trg_points = trg_coord.Dim()/COORD_DIM;
//     static pvfmm::Vector<Real_t> trg_value;
//     pvfmm::Profile::Tic("Evaluation", &sim_config->comm, false, 5);
//     trg_value.Resize(num_trg_points*data_dof);
//     EvalNodesLocal<Real_t, Tree_t>(nodes, trg_coord, trg_value);
//     pvfmm::Profile::Toc();

//     //////////////////////////////////////////////////////////////
//     // GATHERING GLOBAL POINTS VALUES
//     //////////////////////////////////////////////////////////////
//     pvfmm::Profile::Tic("ScatterReverse", &sim_config->comm, true, 5);
//     pvfmm::par::ScatterReverse(trg_value, scatter_index, *tree->Comm(), N);
//     pvfmm::Profile::Toc();

//     //////////////////////////////////////////////////////////////
//     // SETTING EVALUATION VALUES
//     //////////////////////////////////////////////////////////////
//     // memcpy(value, &trg_value[0], trg_value.Dim()*sizeof(Real_t));
//   }  // GLOBAL SORT
  // pvfmm::Profile::Toc();
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
    EvalTree(node_, const_cast<real_t*>(points_pos), num_points, out,sim_config->bc);
    pvfmm::Profile::Toc();
  }

  void operator () (const real_t* points_pos,
                    int num_points,
                    real_t time,
                    real_t* out) {
    (*this)(points_pos, num_points, out);
  }
 private:
  NodeType* node_;
};

}      // namespace tbslas
#endif  // SRC_TREE_NODE_FIELD_FUNCTOR_H_
