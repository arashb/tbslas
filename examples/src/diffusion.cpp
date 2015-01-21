#include <mpi.h>
#include <pvfmm_common.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

#include <profile.hpp>
#include <fmm_cheb.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>

// TBSLAS
#include <utils/common.h>
#include <diffusion/kernel.h>
// #include <tree/>

int NUM_TIME_STEPS = 1;

double TBSLAS_DT;
double TBSLAS_DIFF_COEFF;
double TBSLAS_ALPHA;

const char* OUTPUT_FILE_FORMAT = "%s/%s-VAR_%s-TS_%04d-RNK";
const char* OUTPUT_FILE_PREFIX = "diffusion";

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

//////////////////////////////////////////////////////////////////////////////
// Test1: Laplace problem, Smooth Gaussian, Free Space Boundary
///////////////////////////////////////////////////////////////////////////////
const double a = -160;  // -beta
template <class Real_t>
void fn_input_t1(const Real_t* coord, int n, Real_t* out) { //Input function
  int dof=1;
  Real_t alpha = TBSLAS_ALPHA;
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
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
      Real_t r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
      out[i*dof+0]=exp(a*r_2);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
template <class Real_t>
void RunDiffusion(int test_case, size_t N, size_t M, bool unif, int mult_order,
              int cheb_deg, int depth, bool adap, Real_t tol, MPI_Comm comm) {
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<Real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;

  void (*fn_input_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_poten_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_grad_ )(const Real_t* , int , Real_t*)=NULL;
  const pvfmm::Kernel<Real_t>* mykernel=NULL;
  const pvfmm::Kernel<Real_t>* mykernel_grad=NULL;;
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

  // Find out number of OMP thereads.
  int omp_p=omp_get_max_threads();

  // Find out my identity in the default communicator
  int myrank, p;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm,&p);

  // ======================================================================
  // TREE DATA
  // ======================================================================
  typename FMMNode_t::NodeData tree_data;
  tree_data.dim=COORD_DIM;
  tree_data.max_depth=depth;
  tree_data.cheb_deg=cheb_deg;
  //Set input function pointer
  tree_data.input_fn=fn_input_;
  tree_data.data_dof=mykernel->ker_dim[0];
  tree_data.tol=tol;
  //Set source coordinates.
  std::vector<Real_t> pt_coord;
  if(unif) pt_coord = tbslas::point_distrib<Real_t>(tbslas::UnifGrid,N,comm);
  else pt_coord = tbslas::point_distrib<Real_t>(tbslas::RandElps,N,comm); //RandElps, RandGaus
  tree_data.max_pts=M; // Points per octant.
  tree_data.pt_coord=pt_coord;
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
  // ======================================================================
  // SETUP TREE
  // ======================================================================
  FMM_Tree_t* tree=new FMM_Tree_t(comm);
  tree->Initialize(&tree_data);
  // ======================================================================
  // SETUP FMM
  // ======================================================================
  FMM_Mat_t* fmm_mat=NULL;
  {
    fmm_mat=new FMM_Mat_t;
    fmm_mat->Initialize(mult_order,
                        tree_data.cheb_deg,
                        comm,
                        mykernel);
  }
  // ======================================================================
  // DIFFUSUION SOLVER
  // ======================================================================
  char out_name_buffer[300];
  double al2,rl2,ali,rli;
  for (int ts_counter = 1; ts_counter < NUM_TIME_STEPS+1; ts_counter++) {
    // Setup FMM
    tree->InitFMM_Tree(false,bndry);
    tree->SetupFMM(fmm_mat);
    tree->RunFMM();
    tree->Copy_FMMOutput(); //Copy FMM output to tree Data.

    //Find error in FMM output.
    if (ts_counter == 1) {
      CheckChebOutput<FMM_Tree_t>(tree,
                                  fn_poten_,
                                  mykernel->ker_dim[1],
                                  al2,rl2,ali,rli,
                                  std::string("Output"));
    }
    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             OUTPUT_FILE_FORMAT,
             tbslas::get_result_dir().c_str(),
             OUTPUT_FILE_PREFIX,
             "diff",
             ts_counter);

    //Write2File
    tree->Write2File(out_name_buffer,tree_data.cheb_deg);
  }
  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  int num_leaves = tbslas::CountNumLeafNodes(*tree);
  tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
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

  bool  unif = (commandline_option(argc, argv, "-unif",    NULL, false,
                                   "-unif                : Uniform point distribution."        )!=NULL);
  int      m = strtoul(commandline_option(argc, argv,    "-m",    "10", false,
                                          "-m    <int> = (10)   : Multipole order (+ve even integer)."),NULL,10);
  int   test = strtoul(commandline_option(argc, argv, "-test",     "1", false,
                                          "-test <int> = (1)    : 1) Modified Laplace, Smooth Gaussian, FreeSpace Boundary"),NULL,10);

  tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();

  NUM_TIME_STEPS    = sim_config->total_num_timestep;
  TBSLAS_DT         = sim_config->dt;
  TBSLAS_DIFF_COEFF = 0.01;
  TBSLAS_ALPHA      = (1.0)/TBSLAS_DT/TBSLAS_DIFF_COEFF;

  pvfmm::Profile::Enable(sim_config->profile);
  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }
  // =========================================================================
  // RUN
  // =========================================================================
  pvfmm::Profile::Tic("RunDiffuionSolver",&comm,true);
  RunDiffusion<double>(test,
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
