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
#include <tree/tree_set_functor.h>
#include <tree/tree_extrap_functor.h>

#include <kernels/mod_stokes.h>

double TBSLAS_MOD_STOKES_DIFF_COEFF;
double TBSLAS_MOD_STOKES_ALPHA;
double EXP_ALPHA;

// current simulation time
// double tcurr = 0.25;
double tcurr_init = 0;
double tcurr = 0;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;


template <class real_t>
void RunNS(int test, size_t N, size_t M, bool unif, int mult_order,
	   int cheb_deg, int depth, bool adap, real_t tol, int merge, MPI_Comm comm) {
  typedef double RealType;
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  typedef typename FMM_Tree_t::Node_t NodeType;

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();

  // ======================================================================
  // SETUP FMM KERNEL
  // ======================================================================
  const pvfmm::Kernel<real_t>* mykernel=NULL;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<real_t> modified_stokes_kernel_d =
    pvfmm::BuildKernel<real_t, tbslas::modified_stokes_vel>
    (tbslas::GetModfiedStokesKernelName<real_t>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3));

  // ======================================================================
  // SETUP TEST CASE
  // ======================================================================
  void (*fn_veloc_)(const real_t* , int , real_t*)=NULL;

  switch (test) {
  case 1:
    fn_veloc_ = get_taylor_green_field_ns<real_t>;
    mykernel  = &modified_stokes_kernel_d;
    // bndry = pvfmm::FreeSpace;
    bndry = pvfmm::Periodic;
    break;
  case 2:
    fn_veloc_ = get_two_vortex_tubes_vorticity_field_ns<real_t>;
    mykernel  = &modified_stokes_kernel_d;
    // bndry = pvfmm::FreeSpace;
    bndry = pvfmm::Periodic;
    break;
  case 3:
      fn_veloc_ = get_two_vortex_tubes_curl_vorticity_field_ns<real_t>;
      mykernel  = &modified_stokes_kernel_d;
      // bndry = pvfmm::FreeSpace;
      bndry = pvfmm::Periodic;
    break;
  default:
    fn_veloc_ = NULL;
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

  // =========================================================================
  // INIT THE VELOCITY TREES
  // =========================================================================
  FMM_Tree_t* tvelp = new FMM_Tree_t(comm);
  FMM_Tree_t* tvelc = new FMM_Tree_t(comm);
  {
    tcurr = tcurr_init - sim_config->dt;
    tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				      sim_config->tree_num_points_per_octanct,
				      sim_config->tree_chebyshev_order,
				      sim_config->tree_max_depth,
				      sim_config->tree_adap,
				      sim_config->tree_tolerance,
				      comm,
				      fn_veloc_,
				      3,
				      *tvelp);
    tcurr = tcurr_init;
    tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				      sim_config->tree_num_points_per_octanct,
				      sim_config->tree_chebyshev_order,
				      sim_config->tree_max_depth,
				      sim_config->tree_adap,
				      sim_config->tree_tolerance,
				      comm,
				      fn_veloc_,
				      3,
				      *tvelc);
  }
  if (test == 3){
    // ======================================================================
    // SETUP FMM
    // ======================================================================
    //Initialize FMM_Mat.
    const pvfmm::Kernel<real_t>* mykernel= &pvfmm::StokesKernel<real_t>::velocity();
    FMM_Mat_t* fmm_mat = NULL;
    {
      fmm_mat = new FMM_Mat_t;
      fmm_mat->Initialize(mult_order,
			  cheb_deg,
			  comm,
			  mykernel);
    }
    // =========================================================================
    // RUN FMM
    // =========================================================================
    pvfmm::Profile::Tic("FMM",&comm,true);
    tvelp->InitFMM_Tree(false,bndry);
    tvelp->SetupFMM(fmm_mat);
    tvelp->RunFMM();
    tvelp->Copy_FMMOutput(); //Copy FMM output to tree Data.
    pvfmm::Profile::Toc();

   // =========================================================================
    // RUN FMM
    // =========================================================================
    pvfmm::Profile::Tic("FMM",&comm,true);
    tvelc->InitFMM_Tree(false,bndry);
    tvelc->SetupFMM(fmm_mat);
    tvelc->RunFMM();
    tvelc->Copy_FMMOutput(); //Copy FMM output to tree Data.
    pvfmm::Profile::Toc();
 
  }

  // =========================================================================
  // STORING THE INITIAL VELOCITY TREES
  // =========================================================================
  if (sim_config->vtk_save_rate) {
    tvelp->Write2File(tbslas::GetVTKFileName(0, "velp").c_str(),
		      sim_config->vtk_order);
    tvelc->Write2File(tbslas::GetVTKFileName(0, "vel").c_str(),
		      sim_config->vtk_order);
  }


  real_t in_al2,in_rl2,in_ali,in_rli;
  CheckChebOutput<FMM_Tree_t>(tvelc,
			      fn_veloc_,
			      mykernel->ker_dim[1],
			      in_al2,
			      in_rl2,
			      in_ali,
			      in_rli,
			      std::string("Input"));

  int vel_noct_sum = 0;
  int vel_noct_max = 0;
  int vel_noct_min = 0;
  int tvel_depth     = 0;
  int tvel_depth_max = 0;

  int cfl_sum = 0;
  int cfl_max = 0;
  int cfl_min = 0;

  if (sim_config->profile) {
    // DETERMINE THE NUM OF LEAF OCTANTS
    int vel_noct = tbslas::CountNumLeafNodes(*tvelc);
    vel_noct_sum += vel_noct;
    vel_noct_max = vel_noct;
    vel_noct_min = vel_noct;

    // DETERMINE THE MAX DEPTH OF THREE
    tbslas::GetTreeMaxDepth<FMM_Tree_t>(*tvelc, tvel_depth);
    if (tvel_depth > tvel_depth_max) tvel_depth_max = tvel_depth;

    real_t cfl_val = tbslas::compute_cheb_CFL(sim_config->dt,
					      tvel_depth,
					      sim_config->tree_chebyshev_order);
    cfl_sum += cfl_val;
    cfl_max  = cfl_val;
    cfl_min  = cfl_val;
  }

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
  FMMNode_t* n_curr = tvelc->PostorderFirst();
  while (n_curr != NULL) {
    if(!n_curr->IsGhost() && n_curr->IsLeaf())
      break;
    n_curr = tvelc->PostorderNxt(n_curr);
  }
  int data_dof = n_curr->DataDOF();
  int sdim     = tvelc->Dim();

  // ======================================================================
  // SETUP FMM
  // ======================================================================
  //Initialize FMM_Mat.
  FMM_Mat_t* fmm_mat = NULL;
  {
    fmm_mat = new FMM_Mat_t;
    fmm_mat->Initialize(mult_order,
                        cheb_deg,
                        comm,
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

  for (int timestep = 1; timestep < sim_config->total_num_timestep+1; timestep +=1) {
    // =====================================================================
    // (SEMI) MERGE TO FIX IMBALANCE
    // =====================================================================
    switch(merge) {
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

    if (sim_config->profile) {
      // DETERMINE THE NUM OF LEAF OCTANTS
      int vel_noct = tbslas::CountNumLeafNodes(*treen);
      vel_noct_sum += vel_noct;
      if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
      if (vel_noct < vel_noct_min) vel_noct_min = vel_noct;

      // DETERMINE THE MAX DEPTH OF THREE
      tbslas::GetTreeMaxDepth<FMM_Tree_t>(*treen, tvel_depth);
      if (tvel_depth > tvel_depth_max) tvel_depth_max = tvel_depth;

      real_t cfl_val = tbslas::compute_cheb_CFL(sim_config->dt,
						tvel_depth,
						sim_config->tree_chebyshev_order);
      cfl_sum += cfl_val;
      if (cfl_val > cfl_max) cfl_max = cfl_val;
      if (cfl_val < cfl_min) cfl_min = cfl_val;
    }

    // UPDATE THE SIMULATION CURRENT TIME
    tcurr = tcurr_init + timestep*sim_config->dt;

    // CONSTRUCT THE VELOCITY TREES FUNCTORS
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelp_functor(tvelp);
    tbslas::NodeFieldFunctor<real_t,FMM_Tree_t>   tvelc_functor(tvelc);
    tbslas::FieldExtrapFunctor<real_t,FMM_Tree_t> tvele_functor(tvelp, tvelc);

    pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);
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
		       tcurr - sim_config->dt,
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
		       tcurr - sim_config->dt*2,
		       sim_config->num_rk_step,
		       dprts_points_pos);
	tvelp_functor(dprts_points_pos.data(), treen_num_points,  tconp_points_val.data());
      }
      pvfmm::Profile::Toc();  // SL

      // ===================================
      // COMBINE AND STORE THE SEMILAG VALUES
      // ===================================
      real_t ccoeff = 2.0/sim_config->dt;
      real_t pcoeff = -0.5/sim_config->dt;

#pragma omp parallel for
      for (int i = 0; i < treen_points_val.size(); i++) {
	treen_points_val[i] = ccoeff*tconc_points_val[i] + pcoeff*tconp_points_val[i] ;
      }

      // =========================================================================
      // FIX THE VALUES MEMORY LAYOUT
      // =========================================================================
      int d = cheb_deg+1;
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
				cheb_deg,
				data_dof,
				mt_pnts_val_ml);
      // treen->Write2File(tbslas::GetVTKFileName(timestep, "interm").c_str(), sim_config->vtk_order);
      pvfmm::Profile::Add_FLOP(3*treen_points_val.size()); // for combining the two previous time steps values

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

    // ======================================================================
    // Write2File
    // ======================================================================
    if (sim_config->vtk_save_rate) {
      if ( timestep % sim_config->vtk_save_rate == 0) {
        treen->Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
	tcurr = tcurr_init + timestep*sim_config->dt;
        real_t al2,rl2,ali,rli;
        CheckChebOutput<FMM_Tree_t>(treen,
                                    fn_veloc_,
                                    mykernel->ker_dim[1],
                                    al2,rl2,ali,rli,
                                    std::string("Output_TN" + tbslas::ToString(static_cast<long long>(timestep))));
      }
    }
    tvelp = tvelc;
    tvelc = treen;
  }

  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  tcurr = tcurr_init + sim_config->total_num_timestep*sim_config->dt;
  real_t al2,rl2,ali,rli;
  CheckChebOutput<FMM_Tree_t>(tvelc,
                              fn_veloc_,
                              mykernel->ker_dim[1],
                              al2,rl2,ali,rli,
                              std::string("Output"));


  typedef tbslas::Reporter<real_t> Rep;
  if(!myrank) {
    Rep::AddData("NP", np, tbslas::REP_INT);
    Rep::AddData("OMP", sim_config->num_omp_threads, tbslas::REP_INT);

    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

    Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
    Rep::AddData("VMaxD", tvel_depth_max, tbslas::REP_INT);

    Rep::AddData("DT", sim_config->dt);
    Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);

    Rep::AddData("TEST", test, tbslas::REP_INT);
    Rep::AddData("MERGE", merge, tbslas::REP_INT);

    Rep::AddData("DIFF", sim_config->diff);
    Rep::AddData("ALPHA", TBSLAS_MOD_STOKES_ALPHA);

    Rep::AddData("InAL2", in_al2);
    Rep::AddData("OutAL2", al2);

    Rep::AddData("InRL2", in_rl2);
    Rep::AddData("OutRL2", rl2);

    Rep::AddData("InALINF", in_ali);
    Rep::AddData("OutALINF", ali);

    Rep::AddData("InRLINF", in_rli);
    Rep::AddData("OutRLINF", rli);

    Rep::AddData("VMinNOCT", vel_noct_min, tbslas::REP_INT);
    Rep::AddData("VAvgNOCT", vel_noct_sum/(sim_config->total_num_timestep+1),
                 tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
    Rep::AddData("VMaxNOCT", vel_noct_max, tbslas::REP_INT);

    Rep::AddData("VMinCFL", cfl_min);
    Rep::AddData("VAvgCFL", cfl_sum/(sim_config->total_num_timestep+1)); // NUMBER OF TIMESTEPS + INITIAL TREE
    Rep::AddData("VMaxCFL", cfl_max);

    Rep::Report();
  }

  // ======================================================================
  // CLEAN UP MEMORY
  // ======================================================================
  //Delete matrices.
  if(fmm_mat)
    delete fmm_mat;

  //Delete the tree.
  delete tvelp;
  delete tvelc;
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
  sim_config->vtk_filename_variable = "vel";

  TBSLAS_MOD_STOKES_DIFF_COEFF = sim_config->diff;
  TBSLAS_MOD_STOKES_ALPHA      = 3.0/2.0/sim_config->dt;
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
  pvfmm::Profile::Tic("NS",&comm,true);
  RunNS<double>(test,
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
