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

// template <class real_t>
// void get_taylor_green_field_ns(const real_t* points_pos,
// 			       int num_points,
// 			       real_t* points_values) {

//   real_t A = 2.0/sqrt(3.0);
//   real_t B = 2.0/sqrt(3.0);
//   real_t C = 2.0/sqrt(3.0);
//   real_t a = 2*PI;
//   real_t b = 2*PI;
//   real_t c = 2*PI;
//   real_t theta = 0.0;
//   real_t m = +2*PI/3.0;
//   real_t n = -2*PI/3.0;
//   real_t q = 0;
//   const real_t* p;
//   real_t x,y,z;
//   for (int i = 0; i < num_points; i++) {
//     p = &points_pos[i*COORD_DIM];
//     x = p[0];
//     y = p[1];
//     z = p[2];
//     points_values[i*3+0] = A*sin(theta+m)*sin(a*x)*cos(b*y)*cos(c*z);
//     points_values[i*3+1] = B*sin(theta+n)*cos(a*x)*sin(b*y)*cos(c*z);
//     points_values[i*3+2] = C*sin(theta+q)*cos(a*x)*cos(b*y)*sin(c*z);
//   }
// }


template <class real_t>
void fn_input_t1(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] =  p[1];
    out[i*3+1] = -p[0];
    out[i*3+2] = 0;
  }
}

template <class real_t>
void fn_poten_t1(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] = 0;
    out[i*3+1] = 0;
    out[i*3+2] = -2;
  }
}

template <class real_t>
void fn_input_t3(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] =  0;
    out[i*3+1] = -p[0]*p[0];
    out[i*3+2] = 0;
  }
}

template <class real_t>
void fn_poten_t3(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] = 0;
    out[i*3+1] = 0;
    out[i*3+2] = -2*p[0];
  }
}

template <class real_t>
void fn_input_t4(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] = -p[1];
    out[i*3+1] = p[0]*p[1];
    out[i*3+2] = p[2];
  }
}

template <class real_t>
void fn_poten_t4(const real_t* coord,
                 int num_points,
                 real_t* out) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &coord[i*COORD_DIM];
    out[i*3+0] = 0;
    out[i*3+1] = 0;
    out[i*3+2] = p[1]+1;
  }
}

template <class real_t>
void RunNS(int test, size_t N, size_t M, bool unif, int mult_order,
	   int cheb_deg, int depth, bool adap, real_t tol, int merge, MPI_Comm comm) {
  typedef double RealType;
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  typedef typename FMM_Tree_t::Node_t NodeType;

  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  pvfmm::BoundaryType bndry;

  // // ======================================================================
  // // SETUP FMM KERNEL
  // // ======================================================================
  // const pvfmm::Kernel<real_t>* mykernel=NULL;
  // const pvfmm::Kernel<real_t> modified_stokes_kernel_d =
  //   pvfmm::BuildKernel<real_t, tbslas::modified_stokes_vel>
  //   (tbslas::GetModfiedStokesKernelName<real_t>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3));

  // ======================================================================
  // SETUP TEST CASE
  // ======================================================================
  void (*fn_input_)(const real_t* , int , real_t*)=NULL;
  void (*fn_poten_)(const real_t* , int , real_t*)=NULL;
  void (*fn_veloc_)(const real_t* , int , real_t*)=NULL;

  switch (test) {
  case 1:
    fn_input_ = fn_input_t1<real_t>;
    fn_poten_ = fn_poten_t1<real_t>;
    bndry = pvfmm::Periodic;
    break;
  case 2:
    fn_input_ = fn_input_t3<real_t>;
    fn_poten_ = fn_poten_t3<real_t>;
    bndry = pvfmm::Periodic;
    break;
  case 3:
    fn_input_ = fn_input_t4<real_t>;
    fn_poten_ = fn_poten_t4<real_t>;
    bndry = pvfmm::Periodic;
    break;
  case 4:
    fn_input_ = get_two_vortex_tubes_vorticity_field_ns<real_t>;
    fn_poten_ = get_two_vortex_tubes_curl_vorticity_field_ns<real_t>;
    bndry = pvfmm::Periodic;
    break;
  //   fn_veloc_ = get_two_vortex_tubes_vorticity_field_ns<real_t>;
  //   mykernel  = &modified_stokes_kernel_d;
  //   // bndry = pvfmm::FreeSpace;
  //   bndry = pvfmm::Periodic;
  //   break;
  // case 3:
  //     fn_veloc_ = get_two_vortex_tubes_curl_vorticity_field_ns<real_t>;
  //     mykernel  = &modified_stokes_kernel_d;
  //     // bndry = pvfmm::FreeSpace;
  //     bndry = pvfmm::Periodic;
  //   break;
  default:
    // fn_input__ = NULL;
    // fn_poten__ = NULL;
    // fn_veloc_ = NULL;
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
  // INIT THE TREE
  // =========================================================================
  FMM_Tree_t* ti = new FMM_Tree_t(comm);
  FMM_Tree_t* tf = new FMM_Tree_t(comm);
  FMM_Tree_t* tfa = new FMM_Tree_t(comm);

  tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				    sim_config->tree_num_points_per_octanct,
				    sim_config->tree_chebyshev_order,
				    sim_config->tree_max_depth,
				    sim_config->tree_adap,
				    sim_config->tree_tolerance,
				    comm,
				    fn_input_,
				    3,
				    *ti);

  tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				    sim_config->tree_num_points_per_octanct,
				    sim_config->tree_chebyshev_order,
				    sim_config->tree_max_depth,
				    sim_config->tree_adap,
				    sim_config->tree_tolerance,
				    comm,
				    fn_input_,
				    3,
				    *tf);

  tbslas::ConstructTree<FMM_Tree_t>(sim_config->tree_num_point_sources,
				    sim_config->tree_num_points_per_octanct,
				    sim_config->tree_chebyshev_order,
				    sim_config->tree_max_depth,
				    sim_config->tree_adap,
				    sim_config->tree_tolerance,
				    comm,
				    fn_poten_,
				    3,
				    *tfa);

  // =========================================================================
  // COMPUTE INPUT ERROR
  // =========================================================================
  real_t in_al2,in_rl2,in_ali,in_rli;
  CheckChebOutput<FMM_Tree_t>(ti,
  			      fn_input_,
  			      3,
  			      in_al2,
  			      in_rl2,
  			      in_ali,
  			      in_rli,
  			      std::string("Input"));

  // =========================================================================
  // COMPUTE CURL
  // =========================================================================
  tbslas::ComputeTreeCurl<FMM_Tree_t>(*ti, *tf);

  // =========================================================================
  // COMPUTE OUTPUT ERROR
  // =========================================================================
  real_t al2,rl2,ali,rli;
  CheckChebOutput<FMM_Tree_t>(tf,
                              fn_poten_,
			      3,
                              al2,rl2,ali,rli,
                              std::string("Output"));

  // =========================================================================
  // STORING THE INITIAL VELOCITY TREES
  // =========================================================================
  if (sim_config->vtk_save_rate) {
    ti->Write2File(tbslas::GetVTKFileName(0, "ti").c_str(),
		      sim_config->vtk_order);
    tf->Write2File(tbslas::GetVTKFileName(0, "tf").c_str(),
		      sim_config->vtk_order);
    tfa->Write2File(tbslas::GetVTKFileName(0, "tfa").c_str(),
		      sim_config->vtk_order);

  }


  typedef tbslas::Reporter<real_t> Rep;
  if(!myrank) {
    Rep::AddData("NP", np, tbslas::REP_INT);
    Rep::AddData("OMP", sim_config->num_omp_threads, tbslas::REP_INT);

    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

    Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
    // Rep::AddData("VMaxD", tvel_depth_max, tbslas::REP_INT);

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

    // Rep::AddData("VMinNOCT", vel_noct_min, tbslas::REP_INT);
    // Rep::AddData("VAvgNOCT", vel_noct_sum/(sim_config->total_num_timestep+1),
    //              tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
    // Rep::AddData("VMaxNOCT", vel_noct_max, tbslas::REP_INT);

    // Rep::AddData("VMinCFL", cfl_min);
    // Rep::AddData("VAvgCFL", cfl_sum/(sim_config->total_num_timestep+1)); // NUMBER OF TIMESTEPS + INITIAL TREE
    // Rep::AddData("VMaxCFL", cfl_max);

    Rep::Report();
  }

  // ======================================================================
  // CLEAN UP MEMORY
  // ======================================================================
  //Delete matrices.
  // if(fmm_mat)
  //   delete fmm_mat;

  //Delete the tree.
  delete ti;
  delete tf;
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

  // TBSLAS_MOD_STOKES_DIFF_COEFF = sim_config->diff;
  // TBSLAS_MOD_STOKES_ALPHA      = 3.0/2.0/sim_config->dt;
  // EXP_ALPHA         = exp_alpha;
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
  // pvfmm::Profile::print(&comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
