// *************************************************************************
// Copyright (C) 2016 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_TREE_NS_H_
#define SRC_TREE_NS_H_

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
#include <tree/tree_set_functor.h>
#include <tree/tree_extrap_functor.h>

#include <kernels/mod_stokes.h>

extern double TBSLAS_MOD_STOKES_DIFF_COEFF;
extern double TBSLAS_MOD_STOKES_ALPHA;

template <class real_t, typename TreeType >
void SolveNS1O(TreeType* tvelp,
	       TreeType* tvelc,
	       const real_t tcurr_init,
	       const int num_timestep,
	       const real_t dt) {
  typedef typename TreeType::Node_t NodeType;
  typedef pvfmm::FMM_Cheb<NodeType> FMM_Mat_t;
  typedef TreeType FMM_Tree_t;

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  double tcurr = 0;

  // ======================================================================
  // SETUP FMM KERNEL
  // ======================================================================
  TBSLAS_MOD_STOKES_DIFF_COEFF = sim_config->diff;
  TBSLAS_MOD_STOKES_ALPHA      = 1.0/dt;

  const pvfmm::Kernel<real_t>* mykernel = NULL;
  const pvfmm::Kernel<real_t> modified_stokes_kernel_d =
    pvfmm::BuildKernel<real_t, tbslas::modified_stokes_vel>
    (tbslas::GetModfiedStokesKernelName<real_t>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3));
  mykernel  = &modified_stokes_kernel_d;

  {
    std::vector<NodeType*>  nlist = tvelc->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const real_t* , int , real_t*))NULL;
    }

    nlist = tvelp->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const real_t* , int , real_t*))NULL;
    }
  }

  // ======================================================================
  // GET THE TREE PARAMETERS FROM CURRENT TREE
  // ======================================================================
  NodeType* n_curr = tvelc->PostorderFirst();
  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf())
      break;
    n_curr = tvelc->PostorderNxt(n_curr);
  }
  int data_dof = n_curr->DataDOF();

  // ======================================================================
  // SETUP FMM
  // ======================================================================
  //Initialize FMM_Mat.
  FMM_Mat_t* fmm_mat = NULL;
  {
    fmm_mat = new FMM_Mat_t;
    fmm_mat->Initialize(sim_config->mult_order,
                        sim_config->tree_chebyshev_order,
                        sim_config->comm,
                        mykernel);
  }

  // ======================================================================
  // RUN
  // ======================================================================
  std::vector<real_t> arrvl_points_pos;
  std::vector<real_t> dprts_points_pos;
  std::vector<real_t> tconp_points_val;
  std::vector<real_t> tconc_points_val;
  std::vector<real_t> treen_points_val;

  for (int timestep = 1; timestep < num_timestep+1; timestep +=1) {
    // =====================================================================
    // (SEMI) MERGE TO FIX IMBALANCE
    // =====================================================================
    switch(sim_config->merge) {
    case 2:
      pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
      tbslas::MergeTree(*tvelc, *tvelp);
      pvfmm::Profile::Toc();
      break;
    case 3:
      pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
      tbslas::SemiMergeTree(*tvelc, *tvelp);
      pvfmm::Profile::Toc();
      break;
    }

    // use previous time step's tree for the next time step
    FMM_Tree_t* treen = tvelp;

    // if (sim_config->profile) {
    //   // DETERMINE THE NUM OF LEAF OCTANTS
    //   int vel_noct = tbslas::CountNumLeafNodes(*treen);
    //   vel_noct_sum += vel_noct;
    //   if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
    //   if (vel_noct < vel_noct_min) vel_noct_min = vel_noct;

    //   // DETERMINE THE MAX DEPTH OF THREE
    //   tbslas::GetTreeMaxDepth<FMM_Tree_t>(*treen, tvel_depth);
    //   if (tvel_depth > tvel_depth_max) tvel_depth_max = tvel_depth;

    //   real_t cfl_val = tbslas::compute_cheb_CFL(dt,
    // 						tvel_depth,
    // 						sim_config->tree_chebyshev_order);
    //   cfl_sum += cfl_val;
    //   if (cfl_val > cfl_max) cfl_max = cfl_val;
    //   if (cfl_val < cfl_min) cfl_min = cfl_val;
    // }

    // UPDATE THE SIMULATION CURRENT TIME
    tcurr = tcurr_init + timestep*dt;

    // CONSTRUCT THE VELOCITY TREES FUNCTORS
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelp_functor(tvelp);
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelc_functor(tvelc);
    tbslas::FieldExtrapFunctor<real_t,FMM_Tree_t> tvele_functor(tvelp, tvelc);

    pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &sim_config->comm, true);
    {
      // =====================================================================
      // SOLVE SEMILAG
      // =====================================================================
      // COLLECT THE MERGED TREE POINTS
      int num_leaf = tbslas::CollectChebTreeGridPoints(*treen, arrvl_points_pos);
      int treen_num_points = arrvl_points_pos.size()/COORD_DIM;
      dprts_points_pos.resize(arrvl_points_pos.size());
      tconp_points_val.resize(treen_num_points*data_dof);
      tconc_points_val.resize(treen_num_points*data_dof);
      treen_points_val.resize(treen_num_points*data_dof);

      pvfmm::Profile::Tic("SLM", &sim_config->comm, false, 5);
      {
	// ===================================
	// FIRST STEP BACKWARD TRAJ COMPUTATION
	// ===================================
	ComputeTrajRK2(tvelc_functor,
		       tvele_functor,
		       arrvl_points_pos,
		       tcurr,
		       tcurr - dt,
		       sim_config->num_rk_step,
		       dprts_points_pos);
	tvelc_functor(dprts_points_pos.data(), treen_num_points, tconc_points_val.data());

	/* // =================================== */
	/* // SECOND STEP BACKWARD TRAJ COMPUTATION */
	/* // =================================== */
	/* ComputeTrajRK2(tvelp_functor, */
	/* 	       tvelc_functor, */
	/* 	       arrvl_points_pos, */
	/* 	       tcurr, */
	/* 	       tcurr - dt*2, */
	/* 	       sim_config->num_rk_step, */
	/* 	       dprts_points_pos); */
	/* tvelp_functor(dprts_points_pos.data(), treen_num_points,  tconp_points_val.data()); */
      }
      pvfmm::Profile::Toc();  // SL

      // ===================================
      // COMBINE AND STORE THE SEMILAG VALUES
      // ===================================
      real_t ccoeff = 1.0/dt;
      /* real_t pcoeff = -0.5/dt; */

#pragma omp parallel for
      for (int i = 0; i < treen_points_val.size(); i++) {
	treen_points_val[i] = ccoeff*tconc_points_val[i];// + pcoeff*tconp_points_val[i] ;
      }

      // =========================================================================
      // FIX THE VALUES MEMORY LAYOUT
      // =========================================================================
      int d = sim_config->tree_chebyshev_order+1;
      int num_pnts_per_node = d*d*d;
      std::vector<real_t> mt_pnts_val_ml(treen_num_points*data_dof);
      for (int nindx = 0; nindx < num_leaf; nindx++) {
	int input_shift = nindx*num_pnts_per_node*data_dof;
	for (int j = 0; j < num_pnts_per_node; j++) {
	  for (int i = 0 ; i < data_dof; i++) {
	    mt_pnts_val_ml[input_shift+j+i*num_pnts_per_node] = treen_points_val[input_shift+j*data_dof+i];
	  }
	}
      }

      // ===================================
      // STORE THE VALUES IN TREE
      // ===================================
      tbslas::SetTreeGridValues(*treen,
				sim_config->tree_chebyshev_order,
				data_dof,
				mt_pnts_val_ml);
      pvfmm::Profile::Add_FLOP(3*treen_points_val.size()); // for combining the two previous time steps values

      // =========================================================================
      // RUN FMM
      // =========================================================================
      pvfmm::Profile::Tic("FMM",&sim_config->comm,true);
      treen->InitFMM_Tree(false,sim_config->bc);
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

    // ======================================================================
    // Write2File
    // ======================================================================
    if (sim_config->vtk_save_rate) {
      if ( timestep % sim_config->vtk_save_rate == 0) {
        treen->Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
	tcurr = tcurr_init + timestep*dt;
        real_t al2,rl2,ali,rli;
        // CheckChebOutput<FMM_Tree_t>(treen,
        //                             fn_veloc_,
        //                             mykernel->ker_dim[1],
        //                             al2,rl2,ali,rli,
        //                             std::string("Output_TN" + tbslas::ToString(static_cast<long long>(timestep))));
	// if (test == 3) {
	//   // =========================================================================
	//   // CREATE THE VORTICITY TREE FOR TEST CASE 3
	//   // =========================================================================
	//   FMM_Tree_t* tvort = new FMM_Tree_t(comm);
	//   tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
	// 				    sim_config->tree_num_points_per_octanct,
	// 				    sim_config->tree_chebyshev_order,
	// 				    sim_config->tree_max_depth,
	// 				    sim_config->tree_adap,
	// 				    sim_config->tree_tolerance,
	// 				    comm,
	// 				    get_vorticity_field_wrapper<double>, // used as a dummy function
	// 				    3,
	// 				    *tvort);

	//   tbslas::SyncTreeRefinement(*treen, *tvort);

	//   // =========================================================================
	//   // COMPUTE THE VORTICITY
	//   // =========================================================================
	//   tbslas::ComputeTreeCurl<FMM_Tree_t>(*treen, *tvort);

	//   // =========================================================================
	//   // SAVE THE VORTICITY TREE
	//   // =========================================================================
	//   tvort->Write2File(tbslas::GetVTKFileName(timestep, "vort").c_str(), sim_config->vtk_order);
	// }
      }
    }
    tvelp = tvelc;
    tvelc = treen;
  }

  //Delete matrices.
  if(fmm_mat)
    delete fmm_mat;
}

template <class real_t, typename TreeType >
void SolveNS2O(TreeType* tvelp,
	       TreeType* tvelc,
	       const real_t tcurr_init,
	       const int num_timestep,
	       const real_t dt) {
  typedef typename TreeType::Node_t NodeType;
  typedef pvfmm::FMM_Cheb<NodeType> FMM_Mat_t;
  typedef TreeType FMM_Tree_t;

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  double tcurr = 0;

  // ======================================================================
  // SETUP FMM KERNEL
  // ======================================================================
  TBSLAS_MOD_STOKES_DIFF_COEFF = sim_config->diff;
  TBSLAS_MOD_STOKES_ALPHA      = 3.0/2.0/dt;

  const pvfmm::Kernel<real_t>* mykernel = NULL;
  const pvfmm::Kernel<real_t> modified_stokes_kernel_d =
    pvfmm::BuildKernel<real_t, tbslas::modified_stokes_vel>
    (tbslas::GetModfiedStokesKernelName<real_t>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3));
  mykernel  = &modified_stokes_kernel_d;

  {
    std::vector<NodeType*>  nlist = tvelc->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const real_t* , int , real_t*))NULL;
    }

    nlist = tvelp->GetNodeList();
    for(int i = 0; i < nlist.size(); i++) {
      nlist[i]->input_fn = (void (*)(const real_t* , int , real_t*))NULL;
    }
  }

  // ======================================================================
  // GET THE TREE PARAMETERS FROM CURRENT TREE
  // ======================================================================
  NodeType* n_curr = tvelc->PostorderFirst();
  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf())
      break;
    n_curr = tvelc->PostorderNxt(n_curr);
  }
  int data_dof = n_curr->DataDOF();

  // ======================================================================
  // SETUP FMM
  // ======================================================================
  //Initialize FMM_Mat.
  FMM_Mat_t* fmm_mat = NULL;
  {
    fmm_mat = new FMM_Mat_t;
    fmm_mat->Initialize(sim_config->mult_order,
                        sim_config->tree_chebyshev_order,
                        sim_config->comm,
                        mykernel);
  }

  // ======================================================================
  // RUN
  // ======================================================================
  std::vector<real_t> arrvl_points_pos;
  std::vector<real_t> dprts_points_pos;
  std::vector<real_t> tconp_points_val;
  std::vector<real_t> tconc_points_val;
  std::vector<real_t> treen_points_val;

  for (int timestep = 1; timestep < num_timestep+1; timestep +=1) {
    // =====================================================================
    // (SEMI) MERGE TO FIX IMBALANCE
    // =====================================================================
    switch(sim_config->merge) {
    case 2:
      pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
      tbslas::MergeTree(*tvelc, *tvelp);
      pvfmm::Profile::Toc();
      break;
    case 3:
      pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
      tbslas::SemiMergeTree(*tvelc, *tvelp);
      pvfmm::Profile::Toc();
      break;
    }

    // use previous time step's tree for the next time step
    FMM_Tree_t* treen = tvelp;

    // if (sim_config->profile) {
    //   // DETERMINE THE NUM OF LEAF OCTANTS
    //   int vel_noct = tbslas::CountNumLeafNodes(*treen);
    //   vel_noct_sum += vel_noct;
    //   if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
    //   if (vel_noct < vel_noct_min) vel_noct_min = vel_noct;

    //   // DETERMINE THE MAX DEPTH OF THREE
    //   tbslas::GetTreeMaxDepth<FMM_Tree_t>(*treen, tvel_depth);
    //   if (tvel_depth > tvel_depth_max) tvel_depth_max = tvel_depth;

    //   real_t cfl_val = tbslas::compute_cheb_CFL(dt,
    // 						tvel_depth,
    // 						sim_config->tree_chebyshev_order);
    //   cfl_sum += cfl_val;
    //   if (cfl_val > cfl_max) cfl_max = cfl_val;
    //   if (cfl_val < cfl_min) cfl_min = cfl_val;
    // }

    // UPDATE THE SIMULATION CURRENT TIME
    tcurr = tcurr_init + timestep*dt;

    // CONSTRUCT THE VELOCITY TREES FUNCTORS
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelp_functor(tvelp);
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelc_functor(tvelc);
    tbslas::FieldExtrapFunctor<real_t,FMM_Tree_t> tvele_functor(tvelp, tvelc);

    pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &sim_config->comm, true);
    {
      // =====================================================================
      // SOLVE SEMILAG
      // =====================================================================
      // COLLECT THE MERGED TREE POINTS
      int num_leaf = tbslas::CollectChebTreeGridPoints(*treen, arrvl_points_pos);
      int treen_num_points = arrvl_points_pos.size()/COORD_DIM;
      dprts_points_pos.resize(arrvl_points_pos.size());
      tconp_points_val.resize(treen_num_points*data_dof);
      tconc_points_val.resize(treen_num_points*data_dof);
      treen_points_val.resize(treen_num_points*data_dof);

      pvfmm::Profile::Tic("SLM", &sim_config->comm, false, 5);
      {
	// ===================================
	// FIRST STEP BACKWARD TRAJ COMPUTATION
	// ===================================
	ComputeTrajRK2(tvelc_functor,
		       tvele_functor,
		       arrvl_points_pos,
		       tcurr,
		       tcurr - dt,
		       sim_config->num_rk_step,
		       dprts_points_pos);
	tvelc_functor(dprts_points_pos.data(), treen_num_points, tconc_points_val.data());

	// ===================================
	// SECOND STEP BACKWARD TRAJ COMPUTATION
	// ===================================
	ComputeTrajRK2(tvelp_functor,
		       tvelc_functor,
		       arrvl_points_pos,
		       tcurr,
		       tcurr - dt*2,
		       sim_config->num_rk_step,
		       dprts_points_pos);
	tvelp_functor(dprts_points_pos.data(), treen_num_points,  tconp_points_val.data());
      }
      pvfmm::Profile::Toc();  // SL

      // ===================================
      // COMBINE AND STORE THE SEMILAG VALUES
      // ===================================
      real_t ccoeff = 2.0/dt;
      real_t pcoeff = -0.5/dt;

#pragma omp parallel for
      for (int i = 0; i < treen_points_val.size(); i++) {
	treen_points_val[i] = ccoeff*tconc_points_val[i] + pcoeff*tconp_points_val[i] ;
      }

      // =========================================================================
      // FIX THE VALUES MEMORY LAYOUT
      // =========================================================================
      int d = sim_config->tree_chebyshev_order+1;
      int num_pnts_per_node = d*d*d;
      std::vector<real_t> mt_pnts_val_ml(treen_num_points*data_dof);
      for (int nindx = 0; nindx < num_leaf; nindx++) {
	int input_shift = nindx*num_pnts_per_node*data_dof;
	for (int j = 0; j < num_pnts_per_node; j++) {
	  for (int i = 0 ; i < data_dof; i++) {
	    mt_pnts_val_ml[input_shift+j+i*num_pnts_per_node] = treen_points_val[input_shift+j*data_dof+i];
	  }
	}
      }

      // ===================================
      // STORE THE VALUES IN TREE
      // ===================================
      tbslas::SetTreeGridValues(*treen,
				sim_config->tree_chebyshev_order,
				data_dof,
				mt_pnts_val_ml);
      pvfmm::Profile::Add_FLOP(3*treen_points_val.size()); // for combining the two previous time steps values

      // =========================================================================
      // RUN FMM
      // =========================================================================
      pvfmm::Profile::Tic("FMM",&sim_config->comm,true);
      treen->InitFMM_Tree(false,sim_config->bc);
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

    // ======================================================================
    // Write2File
    // ======================================================================
    if (sim_config->vtk_save_rate) {
      if ( timestep % sim_config->vtk_save_rate == 0) {
        treen->Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
	tcurr = tcurr_init + timestep*dt;
        real_t al2,rl2,ali,rli;
        // CheckChebOutput<FMM_Tree_t>(treen,
        //                             fn_veloc_,
        //                             mykernel->ker_dim[1],
        //                             al2,rl2,ali,rli,
        //                             std::string("Output_TN" + tbslas::ToString(static_cast<long long>(timestep))));
	if (sim_config->test == 3) {
	  // =========================================================================
	  // CREATE THE VORTICITY TREE FOR TEST CASE 3
	  // =========================================================================
	  FMM_Tree_t* tvort = new FMM_Tree_t(sim_config->comm);
	  tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
					    sim_config->tree_num_points_per_octanct,
					    sim_config->tree_chebyshev_order,
					    sim_config->tree_max_depth,
					    sim_config->tree_adap,
					    sim_config->tree_tolerance,
					    sim_config->comm,
					    get_vorticity_field_wrapper<double>, // used as a dummy function
					    3,
					    *tvort);

	  tbslas::SyncTreeRefinement(*treen, *tvort);

	  // =========================================================================
	  // COMPUTE THE VORTICITY
	  // =========================================================================
	  tbslas::ComputeTreeCurl<FMM_Tree_t>(*treen, *tvort);

	  // =========================================================================
	  // SAVE THE VORTICITY TREE
	  // =========================================================================
	  tvort->Write2File(tbslas::GetVTKFileName(timestep, "vort").c_str(), sim_config->vtk_order);
	}
      }
    }
    tvelp = tvelc;
    tvelc = treen;
  }

  //Delete matrices.
  if(fmm_mat)
    delete fmm_mat;
}

#endif  // SRC_UTILS_SIM_CONFIG_H_
