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
// SYSTEM
#include <mpi.h>
#include <pvfmm_common.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <utility>
// PVFMM
#include <profile.hpp>
#include <fmm_cheb.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>
#include <cheb_node.hpp>
// LOCAL
#include <utils.hpp>
#include <field_wrappers.h>
// TBSLAS
#include <utils/common.h>
#include <utils/reporter.h>
#include <utils/fields.h>

#include <tree/tree_semilag.h>
#include <tree/tree_utils.h>

#include <diffusion/kernel.h>

int NUM_TIME_STEPS = 1;

double TBSLAS_DT;
double TBSLAS_DIFF_COEFF;
double TBSLAS_ALPHA;
double EXP_ALPHA;

// current simulation time
// double tcurr = 0.25;
double tcurr_init = 25;
double tcurr = 25;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

// template <class Real_t>
// void get_exp_alpha_field_wrapper(const Real_t* coord,
//                                   int n,
//                                   Real_t* out) {
//   const Real_t xc = 0.5;
//   const Real_t yc = 0.5;
//   const Real_t zc = 0.55;
//   const Real_t R = 0.3;
//   const Real_t alpha = EXP_ALPHA;
//   tbslas::get_exp_alpha_field(coord, n, out, xc, yc, zc, R, alpha);
// }

// template <class Real_t>
// void get_taylor_green_field_wrapper(const Real_t* coord,
//                             int n,
//                             Real_t* out) {
//   tbslas::get_taylor_green_field(coord, n, out);
// }

// template <class Real_t>
// void get_multiple_guassian_kernel_wraper(const Real_t* coord,
//                                          int n,
//                                          Real_t* out) {
//   // FIRST GAUSSIAN
//   const Real_t xc1  = 0.6;
//   const Real_t yc1  = 0.6;
//   const Real_t zc1  = 0.6;
//   std::vector<Real_t> out1(n);
//   tbslas::gaussian_kernel(coord, n, out1.data(), xc1, yc1, zc1);
//   // FIRST GAUSSIAN
//   const Real_t xc2  = 0.4;
//   const Real_t yc2  = 0.4;
//   const Real_t zc2  = 0.4;
//   std::vector<Real_t> out2(n);
//   tbslas::gaussian_kernel(coord, n, out2.data(), xc2, yc2, zc2);
//   // FIRST GAUSSIAN
//   const Real_t xc3  = 0.3;
//   const Real_t yc3  = 0.3;
//   const Real_t zc3  = 0.7;
//   std::vector<Real_t> out3(n);
//   tbslas::gaussian_kernel(coord, n, out3.data(), xc3, yc3, zc3);
//   for (int i = 0; i < n; i++) {
//     out[i] = out1[i] + out2[i] + out3[i];
//   }
// }

// template <class Real_t>
// void get_guassian_kernel_wraper(const Real_t* coord,
//                      int n,
//                      Real_t* out) {
//   const Real_t xc  = 0.7;
//   const Real_t yc  = 0.7;
//   const Real_t zc  = 0.7;
//   tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
// }

// template <class Real_t>
// void get_hopf_field_wrapper(const Real_t* coord,
//                             int n,
//                             Real_t* out) {
//   const Real_t xc = 0.5;
//   const Real_t yc = 0.5;
//   const Real_t zc = 0.5;
//   tbslas::get_hopf_field(coord, n, out, xc, yc, zc);
// }

// template <class Real_t>
// void get_diffusion_kernel_atT(const Real_t* coord,
//                               int n,
//                               Real_t* out) {
//   const Real_t amp = 1e-2;
//   const Real_t xc = 0.5;
//   const Real_t yc = 0.5;
//   const Real_t zc = 0.5;
//   tbslas::diffusion_kernel(coord,
//                            n,
//                            out,
//                            TBSLAS_DIFF_COEFF,
//                            tcurr,
//                            amp,
//                            xc,
//                            yc,
//                            zc);
// }

// template <class Real_t>
// void get_diffusion_kernel_hopf(const Real_t* coord,
//                                int n,
//                                Real_t* out) {
//   const Real_t amp = 1e-2;
//   const Real_t xc = 0.6;
//   const Real_t yc = 0.6;
//   const Real_t zc = 0.6;
//   double time_curr = 2.5;
//   tbslas::diffusion_kernel(coord,
//                            n,
//                            out,
//                            TBSLAS_DIFF_COEFF,
//                            time_curr,
//                            amp,
//                            xc,
//                            yc,
//                            zc);
// }

// template <class Real_t>
// void get_gaussian_kernel_wrapper(const Real_t* coord,
//                      int n,
//                      Real_t* out) {
//   const Real_t xc  = 0.5;
//   const Real_t yc  = 0.5;
//   const Real_t zc  = 0.55;
//   tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
// }

// template <class Real_t>
// void get_diffusion_kernel_atT_hom(const Real_t* coord,
//                                   int n,
//                                   Real_t* out) {
//   const Real_t amp = 1e-2;
//   const Real_t xc = 0.5+(tcurr-tcurr_init)*-0.5;
//   const Real_t yc = 0.5;
//   const Real_t zc = 0.5;
//   tbslas::diffusion_kernel(coord,
//                            n,
//                            out,
//                            TBSLAS_DIFF_COEFF,
//                            tcurr,
//                            amp,
//                            xc,
//                            yc,
//                            zc);
// }

// template <class Real_t>
// void fn_input_t2(const Real_t* coord,
//                  int n,
//                  Real_t* out) {
//   tbslas::gaussian_kernel_diffusion_input(coord,
//                                           n,
//                                           out,
//                                           TBSLAS_ALPHA);
// }

// template <class Real_t>
// void fn_poten_t2(const Real_t* coord,
//                  int n,
//                  Real_t* out) {
//   tbslas::gaussian_kernel(coord,
//                           n,
//                           out);
// }

template <class Real_t>
void RunAdvectDiff(int test, size_t N, size_t M, bool unif, int mult_order,
                   int cheb_deg, int depth, bool adap, Real_t tol, int merge, MPI_Comm comm) {
  typedef double RealType;
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<Real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  typedef typename FMM_Tree_t::Node_t NodeType;

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  // Find out number of OMP thereads.
  int omp_p=omp_get_max_threads();
  // **********************************************************************
  // SETUP FMM KERNEL
  // **********************************************************************
  const pvfmm::Kernel<Real_t>* mykernel=NULL;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<double> modified_laplace_kernel_d =
      pvfmm::BuildKernel<double, tbslas::modified_laplace_poten>
      (tbslas::GetModfiedLaplaceKernelName<double>(TBSLAS_ALPHA), 3, std::pair<int,int>(1,1));

  // **********************************************************************
  // SETUP TEST CASE
  // **********************************************************************
  void (*fn_1)(const double* , int , double*)=NULL;
  fn_1 = tbslas::get_linear_field_y<double,3>;
  void (*fn_input_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_poten_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_veloc_)(const Real_t* , int , Real_t*)=NULL;

  switch (test) {
    case 1:
      fn_input_ = get_diffusion_kernel_atT<Real_t>;
      fn_poten_ = get_diffusion_kernel_atT<Real_t>;
      fn_veloc_ = tbslas::get_vorticity_field<double,3>;
      mykernel  = &modified_laplace_kernel_d;
      // bndry = pvfmm::FreeSpace;
      bndry = pvfmm::Periodic;
      break;
    case 2:
      fn_input_ = fn_input_t2<Real_t>;
      fn_poten_ = fn_poten_t2<Real_t>;
      fn_veloc_ = tbslas::get_vorticity_field<double,3>;
      mykernel  = &modified_laplace_kernel_d;
      // bndry = pvfmm::FreeSpace;
      bndry = pvfmm::Periodic;
      break;
    case 3:
      fn_input_ = get_diffusion_kernel_atT_hom<Real_t>;
      fn_poten_ = get_diffusion_kernel_atT_hom<Real_t>;
      fn_veloc_ = tbslas::get_vel_field_hom_x<double,3>;//tbslas::get_vorticity_field<double,3>,
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 4:
      // tcurr = 0;
      fn_input_ = fn_input_t2<Real_t>;//get_slotted_cylinder<double,3>;
      fn_poten_ = fn_input_t2<Real_t>;//get_slotted_cylinder<double,3>;
      fn_veloc_ = tbslas::get_vel_field_hom_x<double,3>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 5:
      fn_input_ = get_diffusion_kernel_hopf<Real_t>;
      fn_poten_ = get_diffusion_kernel_hopf<Real_t>;
      fn_veloc_ = get_hopf_field_wrapper<double>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 6:
      fn_input_ = get_exp_alpha_field_wrapper<Real_t>;
      fn_poten_ = get_exp_alpha_field_wrapper<Real_t>;
      fn_veloc_ = get_hopf_field_wrapper<double>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 7:
      fn_input_ = get_exp_alpha_field_wrapper<Real_t>;
      fn_poten_ = get_exp_alpha_field_wrapper<Real_t>;
      fn_veloc_ = get_taylor_green_field_wrapper<Real_t>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 8:
      fn_input_ = get_diffusion_kernel_hopf<Real_t>;
      fn_poten_ = get_diffusion_kernel_hopf<Real_t>;
      fn_veloc_ = get_taylor_green_field_wrapper<Real_t>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 9:
      fn_input_ = get_guassian_kernel_wraper<Real_t>;
      fn_poten_ = get_guassian_kernel_wraper<Real_t>;
      fn_veloc_ = get_hopf_field_wrapper<double>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    case 10:
      fn_input_ = get_gaussian_kernel_wrapper<Real_t>;
      fn_poten_ = get_gaussian_kernel_wrapper<Real_t>;
      fn_veloc_ = get_taylor_green_field_wrapper<Real_t>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::Periodic;
      break;
    default:
      fn_input_=NULL;
      fn_poten_=NULL;
      break;
  }
  // Find out my identity in the default communicator
  int myrank, np;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &np);

  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  sim_config->bc = bndry;

  // **********************************************************************
  // SETUP TREE FOR PREVIOUS TIME STEP
  // **********************************************************************
  FMM_Tree_t* treep = new FMM_Tree_t(comm);
  tbslas::ConstructTree<FMM_Tree_t>(N,
                                    M,
                                    cheb_deg,
                                    depth,
                                    adap,
                                    tol,
                                    comm,
                                    fn_input_,
                                    1,
                                    *treep);
  // **********************************************************************
  // SETUP TREE FOR CURRENT TIMESTEP
  // **********************************************************************
  tcurr += TBSLAS_DT;
  FMM_Tree_t* treec = new FMM_Tree_t(comm);
  tbslas::ConstructTree<FMM_Tree_t>(N,
                                    M,
                                    cheb_deg,
                                    depth,
                                    adap,
                                    tol,
                                    comm,
                                    fn_input_,
                                    1,
                                    *treec);

  if (sim_config->vtk_save_rate) {
    //     treep->Write2File(tbslas::GetVTKFileName(-1, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
    treec->Write2File(tbslas::GetVTKFileName(0, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
  }

  // **********************************************************************
  // SETUP VELOCITY FIELD TREE
  // **********************************************************************
  FMM_Tree_t* tvel = new FMM_Tree_t(comm);
  tbslas::ConstructTree<FMM_Tree_t>(N,
                                    M,
                                    cheb_deg,
                                    depth,
                                    adap,
                                    tol,
                                    comm,
                                    fn_veloc_,
                                    3,
                                    *tvel);
  if (sim_config->vtk_save_rate) {
    tvel->Write2File(tbslas::GetVTKFileName(0, "velocity").c_str(), sim_config->vtk_order);
  }

  double in_al2,in_rl2,in_ali,in_rli;
  CheckChebOutput<FMM_Tree_t>(treec,
                              fn_poten_,
                              mykernel->ker_dim[1],
                              in_al2,
                              in_rl2,
                              in_ali,
                              in_rli,
                              std::string("Input"));

  int con_noct_sum = 0;
  int con_noct_max = 0;
  int con_noct_min = 0;

  int vel_noct_sum = 0;
  int vel_noct_max = 0;
  int vel_noct_min = 0;

  if (sim_config->profile) {
    int con_noct = tbslas::CountNumLeafNodes(*treep);
    con_noct_sum += con_noct;
    con_noct_max = con_noct;
    con_noct_min = con_noct;

    int vel_noct = tbslas::CountNumLeafNodes(*tvel);
    vel_noct_sum += vel_noct;
    vel_noct_max = vel_noct;
    vel_noct_min = vel_noct;
  }

  // set the input_fn to NULL -> needed for adaptive refinement
  {
    std::vector<NodeType*>  nlist = treep->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const Real_t* , int , Real_t*))NULL;
    }

    nlist = treec->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const Real_t* , int , Real_t*))NULL;
    }
  }

  // GET THE TREE PARAMETERS FROM CURRENT TREE
  FMMNode_t* n_curr = treec->PostorderFirst();
  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf())
      break;
    n_curr = treec->PostorderNxt(n_curr);
  }
  int data_dof = n_curr->DataDOF();
  int sdim     = treec->Dim();

  int timestep = 1;
  std::vector<double> dprts_points_pos;
  std::vector<double> treep_points_val;
  std::vector<double> treec_points_val;
  std::vector<double> treen_points_val;

  // **********************************************************************
  // SETUP FMM
  // **********************************************************************
  //Initialize FMM_Mat.
  FMM_Mat_t* fmm_mat = NULL;
  {
    fmm_mat = new FMM_Mat_t;
    fmm_mat->Initialize(mult_order,
                        cheb_deg,
                        comm,
                        mykernel);
  }

  for (; timestep < NUM_TIME_STEPS+1; timestep +=1) {
    // =====================================================================
    // (SEMI) MERGE TO FIX IMBALANCE
    // =====================================================================
    switch(merge) {
      case 2:
        pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
        tbslas::MergeTree(*tvel, *treep);
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
        tbslas::MergeTree(*tvel, *treec);
        pvfmm::Profile::Toc();
        break;
      case 3:
        pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
        tbslas::SemiMergeTree(*tvel, *treep);
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
        tbslas::SemiMergeTree(*tvel, *treec);
        pvfmm::Profile::Toc();
        break;
    }

    // use previous time step's tree for the next time step
    FMM_Tree_t* treen = treep;

    if (sim_config->profile) {
      int con_noct = tbslas::CountNumLeafNodes(*treen);
      con_noct_sum += con_noct;
      if (con_noct > con_noct_max) con_noct_max = con_noct;
      if (con_noct < con_noct_min) con_noct_min = con_noct;

      int vel_noct = tbslas::CountNumLeafNodes(*tvel);
      vel_noct_sum += vel_noct;
      if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
      if (vel_noct < vel_noct_min) vel_noct_min = vel_noct;
    }

    // UPDATE THE SIMULATION CURRENT TIME
    tcurr += TBSLAS_DT;

    tbslas::NodeFieldFunctor<double,FMM_Tree_t> vel_evaluator(tvel);
    tbslas::NodeFieldFunctor<double,FMM_Tree_t> trc_evaluator(treec);
    tbslas::NodeFieldFunctor<double,FMM_Tree_t> trp_evaluator(treep);

    pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);
    {
      // =====================================================================
      // SOLVE SEMILAG
      // =====================================================================
      // COLLECT THE MERGED TREE POINTS
      tbslas::CollectChebTreeGridPoints(*treen, dprts_points_pos);
      int treen_num_points = dprts_points_pos.size()/COORD_DIM;
      treep_points_val.resize(treen_num_points*data_dof);
      treec_points_val.resize(treen_num_points*data_dof);
      treen_points_val.resize(treen_num_points*data_dof);

      pvfmm::Profile::Tic("SLM", &sim_config->comm, false, 5);
      {
        // ===================================
        // FIRST STEP BACKWARD TRAJ COMPUTATION
        // ===================================
        ComputeTrajRK2(vel_evaluator,
                       dprts_points_pos,
                       tcurr,
                       tcurr - TBSLAS_DT,
                       sim_config->num_rk_step,
                       dprts_points_pos);
        trc_evaluator(dprts_points_pos.data(), treen_num_points, treec_points_val.data());
        // ===================================
        // SECOND STEP BACKWARD TRAJ COMPUTATION
        // ===================================
        ComputeTrajRK2(vel_evaluator,
                       dprts_points_pos,
                       tcurr - TBSLAS_DT,
                       tcurr - TBSLAS_DT*2,
                       sim_config->num_rk_step,
                       dprts_points_pos);
        trp_evaluator(dprts_points_pos.data(), treen_num_points,  treep_points_val.data());
        // ===================================
        // COMBINE AND STORE THE SEMILAG VALUES
        // ===================================
        double ccoeff = 4.0/3;
        double pcoeff = 1.0/3;

#pragma omp parallel for
        for (int i = 0; i < treen_points_val.size(); i++) {
          treen_points_val[i] = ccoeff*treec_points_val[i] - pcoeff*treep_points_val[i] ;
        }

        tbslas::SetTreeGridValues(*treen,
                                  cheb_deg,
                                  data_dof,
                                  treen_points_val);

        // FMMNode_t* n_next = treen->PostorderFirst();
        // while (n_next != NULL) {
        //   if(!n_next->IsGhost() && n_next->IsLeaf()) break;
        //   n_next = treen->PostorderNxt(n_next);
        // }
        // std::vector<NodeType*> nodes;
        // while (n_next != NULL) {
        //   if (n_next->IsLeaf() && !n_next->IsGhost()) nodes.push_back(n_next);
        //   n_next = treen->PostorderNxt(n_next);
        // }

        //int num_points_per_node = (cheb_deg+1)*(cheb_deg+1)*(cheb_deg+1);
        //int tree_next_node_counter = 0;
        //while (n_next != NULL) {
        //  if (n_next->IsLeaf() && !n_next->IsGhost()) {
        //    tbslas::NewPt2ChebPt<double>(&treen_points_val[tree_next_node_counter*num_points_per_node*data_dof],
        //                               cheb_deg, data_dof);
        //    pvfmm::cheb_approx<double, double>(&treen_points_val[tree_next_node_counter*num_points_per_node*data_dof],
        //                                     cheb_deg,
        //                                     data_dof,
        //                                     &(n_next->ChebData()[0]));
        //    tree_next_node_counter++;
        //  }
        //  n_next = treen->PostorderNxt(n_next);
        //}

//         int omp_p=omp_get_max_threads();
//         static pvfmm::Matrix<RealType> M;
//         tbslas::GetPt2CoeffMatrix<RealType>(cheb_deg, M);
//         int num_points_per_node=M.Dim(0);
//         pvfmm::Matrix<RealType> Mvalue(treen_points_val.size()/num_points_per_node,M.Dim(0),&treen_points_val[0],false);
//         pvfmm::Matrix<RealType> Mcoeff(treen_points_val.size()/num_points_per_node,M.Dim(1));
// #pragma omp parallel for schedule(static)
//         for(int pid=0;pid<omp_p;pid++){
//           long a=(pid+0)*nodes.size()/omp_p;
//           long b=(pid+1)*nodes.size()/omp_p;
//           pvfmm::Matrix<RealType> Mi((b-a)*data_dof, Mvalue.Dim(1), &Mvalue[a*data_dof][0], false);
//           pvfmm::Matrix<RealType> Mo((b-a)*data_dof, Mcoeff.Dim(1), &Mcoeff[a*data_dof][0], false);
//           pvfmm::Matrix<RealType>::GEMM(Mo, Mi, M);
//           for(long j=0;j<b-a;j++){
//             memcpy(&(nodes[a+j]->ChebData()[0]), &Mo[j*data_dof][0], M.Dim(1)*data_dof*sizeof(RealType));
//           }
//         }

        pvfmm::Profile::Add_FLOP(3*treen_points_val.size()); // for combining the two previous time steps values
      }
      pvfmm::Profile::Toc();  // SL

      // =========================================================================
      // RUN FMM
      // =========================================================================
      pvfmm::Profile::Tic("FMM",&comm,true);
      treen->InitFMM_Tree(false,bndry);
      treen->SetupFMM(fmm_mat);
      treen->RunFMM();
      treen->Copy_FMMOutput(); //Copy FMM output to tree Data.
      pvfmm::Profile::Toc();

    }
    pvfmm::Profile::Toc();        // solve

    // =====================================================================
    // REFINE TREE
    // =====================================================================
    pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
    treen->RefineTree();
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("Balance21", &sim_config->comm, false, 5);
    treen->Balance21(sim_config->bc);
    pvfmm::Profile::Toc();

    //TODO: ONLY FOR STEADY VELOCITY TREES
    tvel->RefineTree();

    // ======================================================================
    // Write2File
    // ======================================================================
    if (sim_config->vtk_save_rate) {
      if ( timestep % sim_config->vtk_save_rate == 0) {
        treen->Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
        double al2,rl2,ali,rli;
        CheckChebOutput<FMM_Tree_t>(treen,
                                    fn_poten_,
                                    mykernel->ker_dim[1],
                                    al2,rl2,ali,rli,
                                    std::string("Output_TN" + tbslas::ToString(static_cast<long long>(timestep))));
      }
    }
    treep = treec;
    treec = treen;
  }

  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  double al2,rl2,ali,rli;
  CheckChebOutput<FMM_Tree_t>(treec,
                              fn_poten_,
                              mykernel->ker_dim[1],
                              al2,rl2,ali,rli,
                              std::string("Output"));

  int tcon_max_depth=0;
  int tvel_max_depth=0;
  tbslas::GetTreeMaxDepth<FMM_Tree_t>(*treec, tcon_max_depth);
  tbslas::GetTreeMaxDepth<FMM_Tree_t>(*tvel, tvel_max_depth);

  typedef tbslas::Reporter<Real_t> Rep;
  if(!myrank) {
    Rep::AddData("NP", np, tbslas::REP_INT);
    Rep::AddData("OMP", sim_config->num_omp_threads, tbslas::REP_INT);

    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

    Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
    Rep::AddData("CMaxD", tcon_max_depth, tbslas::REP_INT);
    Rep::AddData("VMaxD", tvel_max_depth, tbslas::REP_INT);

    Rep::AddData("CBC", sim_config->use_cubic?1:0, tbslas::REP_INT);
    Rep::AddData("CUF", sim_config->cubic_upsampling_factor, tbslas::REP_INT);

    Rep::AddData("DT", sim_config->dt);
    Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);

    Rep::AddData("TEST", test, tbslas::REP_INT);
    Rep::AddData("MERGE", merge, tbslas::REP_INT);

    Rep::AddData("DIFF", TBSLAS_DIFF_COEFF);
    Rep::AddData("ALPHA", TBSLAS_ALPHA);

    Rep::AddData("InAL2", in_al2);
    Rep::AddData("OutAL2", al2);

    Rep::AddData("InRL2", in_rl2);
    Rep::AddData("OutRL2", rl2);

    Rep::AddData("InALINF", in_ali);
    Rep::AddData("OutALINF", ali);

    Rep::AddData("InRLINF", in_rli);
    Rep::AddData("OutRLINF", rli);

    Rep::AddData("CMinNOCT", con_noct_min, tbslas::REP_INT);
    Rep::AddData("CAvgNOCT", con_noct_sum/(sim_config->total_num_timestep+1),
                 tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
    Rep::AddData("CMaxNOCT", con_noct_max, tbslas::REP_INT);

    Rep::AddData("VMinNOCT", vel_noct_min, tbslas::REP_INT);
    Rep::AddData("VAvgNOCT", vel_noct_sum/(sim_config->total_num_timestep+1),
                 tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
    Rep::AddData("VMaxNOCT", vel_noct_max, tbslas::REP_INT);

    Rep::Report();
  }

  // **********************************************************************
  // CLEAN UP MEMORY
  // **********************************************************************
  //Delete matrices.
  if(fmm_mat)
    delete fmm_mat;

  //Delete the tree.
  delete treep;
  delete treec;
  delete tvel;
}

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  parse_command_line_options(argc, argv);

  bool  unif=              (commandline_option(argc, argv, "-unif",    NULL, false, "-unif                : Uniform point distribution."        )!=NULL);
  int      m=       strtoul(commandline_option(argc, argv,    "-m",    "10", false, "-m    <int> = (10)   : Multipole order (+ve even integer)."),NULL,10);
  int   test=       strtoul(commandline_option(argc, argv, "-test",     "1", false,
                                               "-test <int> = (1)    : 1) Laplace, Smooth Gaussian, Periodic Boundary"),NULL,10);
  int   merge = strtoul(commandline_option(argc, argv, "-merge",     "1", false,
                                           "-merge <int> = (1)    : 1) no merge 2) complete merge 3) Semi-Merge"),NULL,10);
  double exp_alpha = strtod(commandline_option(argc, argv,  "-ea",  "10", false, "-ea <real> = (10) : diffusivity" ), NULL);

  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  tbslas::SimConfig* sim_config     = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Enable(sim_config->profile);
  sim_config->vtk_filename_variable = "conc";

  NUM_TIME_STEPS    = sim_config->total_num_timestep;
  TBSLAS_DT         = sim_config->dt;
  TBSLAS_DIFF_COEFF = sim_config->diff;
  TBSLAS_ALPHA      = 3.0/2.0/TBSLAS_DT/TBSLAS_DIFF_COEFF;
  EXP_ALPHA         = exp_alpha;
  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }
  // =========================================================================
  // RUN
  // =========================================================================
  pvfmm::Profile::Tic("AdvDif",&comm,true);
  RunAdvectDiff<double>(test,
                        sim_config->tree_num_point_sources,
                        sim_config->tree_num_points_per_octanct,
                        unif,
                        m,
                        sim_config->tree_chebyshev_order,
                        sim_config->tree_max_depth,
                        sim_config->tree_adap,
                        sim_config->tree_tolerance,
                        merge,
                        comm);
  pvfmm::Profile::Toc();
  //Output Profiling results.
  pvfmm::Profile::print(&comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
