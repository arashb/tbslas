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

#include <mpi.h>
#include <pvfmm_common.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <utility>

#include <profile.hpp>
#include <fmm_cheb.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>

// TBSLAS
#include <utils/common.h>
#include <diffusion/kernel.h>
#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

int NUM_TIME_STEPS = 1;

double TBSLAS_DT;
double TBSLAS_DIFF_COEFF;
double TBSLAS_ALPHA;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

//////////////////////////////////////////////////////////////////////////////
// Test1: Laplace problem, Smooth Gaussian (exp(-beta*r^2)), FreeSpace Boundary
///////////////////////////////////////////////////////////////////////////////
const double xc = 0.5;
const double yc = 0.5;
const double zc = 0.5;
const double a = -160;  // -beta
template <class Real_t>
void fn_input_t1(const Real_t* coord, int n, Real_t* out) { //Input function
  int dof=1;
  Real_t alpha = TBSLAS_ALPHA;
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-xc)*(c[0]-xc)+(c[1]-yc)*(c[1]-yc)+(c[2]-zc)*(c[2]-zc);
      out[i*dof+0]=-(2*a*r_2+3)*2*a*exp(a*r_2)/alpha+exp(a*r_2);
    }
  }
}

template <class Real_t>
void fn_poten_t1(const Real_t* coord, int n, Real_t* out) { //Output potential
  int dof=1;
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-xc)*(c[0]-xc)+(c[1]-yc)*(c[1]-yc)+(c[2]-yc)*(c[2]-yc);
      out[i*dof+0]=exp(a*r_2);
    }
  }
}

template<class _FMM_Tree_Type,
         class _FMM_Mat_Type>
void
RunFMM(_FMM_Tree_Type* tree,
       _FMM_Mat_Type* fmm_mat,
       pvfmm::BoundaryType bndry) {
  tree->InitFMM_Tree(false,bndry);
  tree->SetupFMM(fmm_mat);
  tree->RunFMM();
  tree->Copy_FMMOutput(); //Copy FMM output to tree Data.
}

void
GetVTKFileName(char (*buffer)[300], int timestep) {
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  snprintf(*buffer,
           sizeof(*buffer),
           sim_config->vtk_filename_format.c_str(),
           tbslas::get_result_dir().c_str(),
           sim_config->vtk_filename_prefix.c_str(),
           sim_config->vtk_filename_variable.c_str(),
           timestep);
}
///////////////////////////////////////////////////////////////////////////////
template <class Real_t>
void RunAdvectDiff(int test_case, size_t N, size_t M, bool unif, int mult_order,
                   int cheb_deg, int depth, bool adap, Real_t tol, MPI_Comm comm) {
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<Real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  tbslas::SimConfig* sim_config = tbslas::SimConfigSingleton::Instance();
  char out_name_buffer[300];
  // Find out number of OMP thereads.
  int omp_p=omp_get_max_threads();
  // **********************************************************************
  // SETUP FMM KERNEL
  // **********************************************************************
  void (*fn_input_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_poten_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_grad_ )(const Real_t* , int , Real_t*)=NULL;
  const pvfmm::Kernel<Real_t>* mykernel=NULL;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<double> modified_laplace_kernel_d =
      pvfmm::BuildKernel<double, tbslas::modified_laplace_poten>
      (tbslas::GetModfiedLaplaceKernelName<double>(TBSLAS_ALPHA), 3, std::pair<int,int>(1,1));

  switch (test_case) {
    case 1:
      fn_input_ = fn_input_t1<Real_t>;
      fn_poten_ = fn_poten_t1<Real_t>;
      mykernel  = &modified_laplace_kernel_d;
      bndry = pvfmm::FreeSpace;
      break;
    default:
      fn_input_=NULL;
      fn_poten_=NULL;
      fn_grad_ =NULL;
      break;
  }
  // Find out my identity in the default communicator
  int myrank, p;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm,&p);
  // **********************************************************************
  // SETUP TREE
  // **********************************************************************
  typename FMMNode_t::NodeData tree_data;
  tree_data.dim       = COORD_DIM;
  tree_data.max_depth = depth;
  tree_data.cheb_deg  = cheb_deg;
  tree_data.input_fn  = fn_input_;
  tree_data.data_dof  = mykernel->ker_dim[0];
  tree_data.tol       = tol;
  //Set source coordinates.
  std::vector<Real_t> pt_coord;
  if(unif) pt_coord = tbslas::point_distrib<Real_t>(tbslas::UnifGrid,N,comm);
  else pt_coord = tbslas::point_distrib<Real_t>(tbslas::RandElps,N,comm); //RandElps, RandGaus
  tree_data.max_pts  = M; // Points per octant.
  tree_data.pt_coord = pt_coord;
  //Print various parameters.
  if (!myrank) {
    std::cout<<std::setprecision(2)<<std::scientific;
    std::cout<<"Number of MPI processes: "<<p<<'\n';
    std::cout<<"Number of OpenMP threads: "<<omp_p<<'\n';
    std::cout<<"Order of multipole expansions: "<<mult_order<<'\n';
    std::cout<<"Order of Chebyshev polynomials: "<<tree_data.cheb_deg<<'\n';
    std::cout<<"FMM Kernel name: "<<mykernel->ker_name<<'\n';
    std::cout<<"Number of point samples: "<<N<<'\n';
    std::cout<<"Uniform distribution: "<<(unif?"true":"false")<<'\n';
    std::cout<<"Maximum points per octant: "<<tree_data.max_pts<<'\n';
    std::cout<<"Chebyshev Tolerance: "<<tree_data.tol<<'\n';
    std::cout<<"Maximum Tree Depth: "<<depth<<'\n';
    std::cout<<"BoundaryType: "<<(bndry==pvfmm::Periodic?"Periodic":"FreeSpace")<<'\n';
  }
  //Create Tree and initialize with input data.
  FMM_Tree_t* tree = new FMM_Tree_t(comm);
  tree->Initialize(&tree_data);
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
  // **********************************************************************
  // VELOCITY FIELD
  // **********************************************************************
  FMM_Tree_t tvel_curr(comm);
  tbslas::ConstructTree<FMM_Tree_t>(N,
                                    M,
                                    cheb_deg,
                                    depth,
                                    adap,
                                    tol,
                                    comm,
                                    tbslas::get_vorticity_field<double,3>,
                                    3,
                                    tvel_curr);
  snprintf(out_name_buffer,
           sizeof(out_name_buffer),
           sim_config->vtk_filename_format.c_str(),
           tbslas::get_result_dir().c_str(),
           sim_config->vtk_filename_prefix.c_str(),
           "velocity",
           0);
  tvel_curr.Write2File(out_name_buffer, cheb_deg);
  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  int timestep = 1;
  double al2,rl2,ali,rli;
  for (; timestep < 2*NUM_TIME_STEPS+1; timestep +=2) {
    // **********************************************************************
    // SOLVE DIFFUSION: FMM
    // **********************************************************************
    RunFMM(tree, fmm_mat, bndry);
    // Write2File
    GetVTKFileName(&out_name_buffer, timestep);
    tree->Write2File(out_name_buffer, sim_config->vtk_order);
    // **********************************************************************
    // SOLVE ADVECTION: SEMI-LAGRANGIAN
    // **********************************************************************
    tbslas::RunSemilagSimulationInSitu(&tvel_curr,
                                       tree,
                                       1,
                                       sim_config->dt,
                                       sim_config->num_rk_step,
                                       true,
                                       false);
    //Find error in FMM output.
    if (timestep == 1) {
      CheckChebOutput<FMM_Tree_t>(tree,
                                  fn_poten_,
                                  mykernel->ker_dim[1],
                                  al2,rl2,ali,rli,
                                  std::string("Output"));
    }
    //Write2File
    GetVTKFileName(&out_name_buffer, timestep+1);
    tree->Write2File(out_name_buffer, sim_config->vtk_order);
  }
  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  int num_leaves = tbslas::CountNumLeafNodes(*tree);
  if(!myrank) {
    printf("#TBSLAS-HEADER: %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s %-15s\n",
           "TOL",
           "DT",
           "TN",
           "DIFF",
           "ALPHA",
           "AL2",
           "RL2",
           "ALINF",
           "RLINF",
           "NOCT");
    printf("#TBSLAS-RESULT: %-15.5e %-15.5e %-15d %-15.5e %-15.5e %-15.5e %-15.5e %-15.5e %-15.5e %-15d\n",
           sim_config->tree_tolerance,
           sim_config->dt,
           sim_config->total_num_timestep,
           TBSLAS_DIFF_COEFF,
           TBSLAS_ALPHA,
           al2,
           rl2,
           ali,
           rli,
           num_leaves
           );
  }
  // **********************************************************************
  // CLEAN UP MEMORY
  // **********************************************************************
  //Delete matrices.
  if(fmm_mat)
    delete fmm_mat;

  //Delete the tree.
  delete tree;
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

  // =========================================================================
  // SIMULATION PARAMETERS
  // =========================================================================
  tbslas::SimConfig* sim_config     = tbslas::SimConfigSingleton::Instance();

  pvfmm::Profile::Enable(sim_config->profile);

  sim_config->vtk_filename_prefix   = "advdiff";
  sim_config->vtk_filename_variable = "conc";

  NUM_TIME_STEPS    = sim_config->total_num_timestep;
  TBSLAS_DT         = sim_config->dt;
  TBSLAS_DIFF_COEFF = 0.01;
  TBSLAS_ALPHA      = (1.0)/TBSLAS_DT/TBSLAS_DIFF_COEFF;
  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }
  // =========================================================================
  // RUN
  // =========================================================================
  pvfmm::Profile::Tic("RunAdvectDiff",&comm,true);
  RunAdvectDiff<double>(test,
                        sim_config->tree_num_point_sources,
                        sim_config->tree_num_points_per_octanct,
                        unif,
                        m,
                        sim_config->tree_chebyshev_order,
                        sim_config->tree_max_depth,
                        sim_config->tree_adap,
                        sim_config->tree_tolerance,
                        comm);
  pvfmm::Profile::Toc();

  //Output Profiling results.
  pvfmm::Profile::print(&comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
