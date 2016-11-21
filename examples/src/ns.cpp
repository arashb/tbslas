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

#include <tree/tree_ns.h>

double TBSLAS_MOD_STOKES_DIFF_COEFF;
double TBSLAS_MOD_STOKES_ALPHA;
double EXP_ALPHA;

// current simulation time
double tcurr_init = 0;
double tcurr = 0;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;


template <class real_t>
void RunNS() {
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  MPI_Comm comm = sim_config->comm;

  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  typedef typename FMM_Tree_t::Node_t NodeType;
  
  // ======================================================================
  // SETUP TEST CASE
  // ======================================================================
  void (*fn_veloc_)(const real_t* , int , real_t*)=NULL;
  pvfmm::BoundaryType bndry;

  switch (sim_config->test) {
  case 1:
    fn_veloc_ = get_taylor_green_field_ns<real_t>;
    // bndry = pvfmm::FreeSpace;
    // mykernel  = &modified_stokes_kernel_d;
    bndry = pvfmm::Periodic;
    break;
  case 2:
    fn_veloc_ = get_two_vortex_tubes_vorticity_field_ns<real_t>;
    // mykernel  = &modified_stokes_kernel_d;
    // bndry = pvfmm::FreeSpace;
    bndry = pvfmm::Periodic;
    break;
  case 3:
      fn_veloc_ = get_two_vortex_tubes_curl_vorticity_field_ns<real_t>;
      // mykernel  = &modified_stokes_kernel_d;
      // bndry = pvfmm::FreeSpace;
      bndry = pvfmm::Periodic;
    break;
  default:
    fn_veloc_ = NULL;
    break;
  }
  // Find out my identity in the default communicator
  int myrank, np;
  MPI_Comm_rank(sim_config->comm, &myrank);
  MPI_Comm_size(sim_config->comm, &np);

  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  sim_config->bc = bndry;

  // =========================================================================
  // INIT THE VELOCITY TREES
  // =========================================================================
  FMM_Tree_t* tvelp = new FMM_Tree_t(sim_config->comm);
  FMM_Tree_t* tvelc = new FMM_Tree_t(sim_config->comm);
  {
    tcurr = tcurr_init - sim_config->dt;
    tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				      sim_config->tree_num_points_per_octanct,
				      sim_config->tree_chebyshev_order,
				      sim_config->tree_max_depth,
				      sim_config->tree_adap,
				      sim_config->tree_tolerance,
				      sim_config->comm,
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
				      sim_config->comm,
				      fn_veloc_,
				      3,
				      *tvelc);
  }
  if (sim_config->test == 3){
    // ======================================================================
    // SETUP FMM
    // ======================================================================
    //Initialize FMM_Mat.
    const pvfmm::Kernel<real_t>* mykernel= &pvfmm::StokesKernel<real_t>::velocity();
    FMM_Mat_t* fmm_mat = NULL;
    {
      fmm_mat = new FMM_Mat_t;
      fmm_mat->Initialize(sim_config->mult_order,
			  sim_config->tree_chebyshev_order,
			  sim_config->comm,
			  mykernel);
    }
    // =========================================================================
    // RUN FMM
    // =========================================================================
    pvfmm::Profile::Tic("FMM",&sim_config->comm,true);
    tvelp->InitFMM_Tree(false,bndry);
    tvelp->SetupFMM(fmm_mat);
    tvelp->RunFMM();
    tvelp->Copy_FMMOutput(); //Copy FMM output to tree Data.
    pvfmm::Profile::Toc();

   // =========================================================================
    // RUN FMM
    // =========================================================================
    pvfmm::Profile::Tic("FMM",&sim_config->comm,true);
    tvelc->InitFMM_Tree(false,bndry);
    tvelc->SetupFMM(fmm_mat);
    tvelc->RunFMM();
    tvelc->Copy_FMMOutput(); //Copy FMM output to tree Data.
    pvfmm::Profile::Toc();

    delete fmm_mat;
 
  }

    // =========================================================================
    // STORING THE INITIAL VELOCITY TREES
    // =========================================================================
    if (sim_config->vtk_save_rate) {
      tvelp->Write2File(tbslas::GetVTKFileName(0, "velp").c_str(),
			sim_config->vtk_order);
      tvelc->Write2File(tbslas::GetVTKFileName(0, "vel").c_str(),
			sim_config->vtk_order);
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

	tbslas::SyncTreeRefinement(*tvelc, *tvort);

	// =========================================================================
	// COMPUTE THE VORTICITY
	// =========================================================================
	tbslas::ComputeTreeCurl<FMM_Tree_t>(*tvelc, *tvort);

	// =========================================================================
	// SAVE THE VORTICITY TREE
	// =========================================================================
	tvort->Write2File(tbslas::GetVTKFileName(0, "vort").c_str(), sim_config->vtk_order);

	delete tvort;
      }
    }


  real_t in_al2,in_rl2,in_ali,in_rli;
  // CheckChebOutput<FMM_Tree_t>(tvelc,
  // 			      fn_veloc_,
  // 			      mykernel->ker_dim[1],
  // 			      in_al2,
  // 			      in_rl2,
  // 			      in_ali,
  // 			      in_rli,
  // 			      std::string("Input"));

  int vel_noct_sum = 0;
  int vel_noct_max = 0;
  int vel_noct_min = 0;
  int tvel_depth     = 0;
  int tvel_depth_max = 0;

  int cfl_sum = 0;
  int cfl_max = 0;
  int cfl_min = 0;

  // if (sim_config->profile) {
  //   // DETERMINE THE NUM OF LEAF OCTANTS
  //   int vel_noct = tbslas::CountNumLeafNodes(*tvelc);
  //   vel_noct_sum += vel_noct;
  //   vel_noct_max = vel_noct;
  //   vel_noct_min = vel_noct;

  //   // DETERMINE THE MAX DEPTH OF THREE
  //   tbslas::GetTreeMaxDepth<FMM_Tree_t>(*tvelc, tvel_depth);
  //   if (tvel_depth > tvel_depth_max) tvel_depth_max = tvel_depth;

  //   real_t cfl_val = tbslas::compute_cheb_CFL(sim_config->dt,
  // 					      tvel_depth,
  // 					      sim_config->tree_chebyshev_order);
  //   cfl_sum += cfl_val;
  //   cfl_max  = cfl_val;
  //   cfl_min  = cfl_val;
  // }


  // SolveNS1O<real_t>(tvelp,
  // 		    tvelc, 
  // 		    tcurr_init,
  // 		    sim_config->total_num_timestep,
  // 		    sim_config->dt);

  SolveNS2O<real_t>(tvelp,
  		    tvelc, 
  		    tcurr_init,
  		    sim_config->total_num_timestep,
  		    sim_config->dt);

  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  tcurr = tcurr_init + sim_config->total_num_timestep*sim_config->dt;
  real_t al2,rl2,ali,rli;
  // CheckChebOutput<FMM_Tree_t>(tvelc,
  //                             fn_veloc_,
  //                             mykernel->ker_dim[1],
  //                             al2,rl2,ali,rli,
  //                             std::string("Output"));


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

    Rep::AddData("TEST", sim_config->test, tbslas::REP_INT);
    Rep::AddData("MERGE", sim_config->merge, tbslas::REP_INT);

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

  //Delete the tree.
  delete tvelp;
  delete tvelc;
}

template<class real_t>
void SetupFMMPrecomp() {
  tbslas::SimConfig* sim_config     = tbslas::SimConfigSingleton::Instance();

  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  typedef typename FMM_Tree_t::Node_t NodeType;

  // ======================================================================
  // SETUP FMM KERNEL
  // ======================================================================
  {
    TBSLAS_MOD_STOKES_DIFF_COEFF = sim_config->diff;
    TBSLAS_MOD_STOKES_ALPHA      = 3.0/2.0/sim_config->dt;

    const pvfmm::Kernel<real_t>* mykernel = NULL;
    const pvfmm::Kernel<real_t> modified_stokes_kernel_d =
      pvfmm::BuildKernel<real_t, tbslas::modified_stokes_vel>
      (tbslas::GetModfiedStokesKernelName<real_t>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3),
       NULL, NULL, NULL,
       NULL, NULL, NULL,
       NULL, NULL, NULL,
       false);
    mykernel  = &modified_stokes_kernel_d;

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
    if (fmm_mat)
      delete fmm_mat;
  }
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
  sim_config->mult_order = m;
  sim_config->merge = merge;
  sim_config->test = test;

  EXP_ALPHA         = exp_alpha;

  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }

  // SetupFMMPrecomp<double>();

  // =========================================================================
  // RUN
  // =========================================================================
  // pvfmm::Profile::Tic("NS",&comm,true);
  // RunNSFirstOrder<double>(test,
  // 			  sim_config->tree_num_point_sources,
  // 			  sim_config->tree_num_points_per_octanct,
  // 			  unif,
  // 			  m,
  // 			  sim_config->tree_chebyshev_order,
  // 			  sim_config->tree_max_depth,
  // 			  sim_config->tree_adap,
  // 			  sim_config->tree_tolerance,
  // 			  merge,
  // 			  comm);
  // pvfmm::Profile::Toc();

  // =========================================================================
  // RUN
  // =========================================================================
  pvfmm::Profile::Tic("NS",&comm,true);
  RunNS<double>();
  // sim_config->test,
  // 		sim_config->tree_num_point_sources,
  // 		sim_config->tree_num_points_per_octanct,
  // 		unif,
  // 		sim_config->mult_order,
  // 		sim_config->tree_chebyshev_order,
  // 		sim_config->tree_max_depth,
  // 		sim_config->tree_adap,
  // 		sim_config->tree_tolerance,
  // 		sim_config->merge,
  // 		comm);
  pvfmm::Profile::Toc();
  //Output Profiling results.
  // pvfmm::Profile::print(&comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
