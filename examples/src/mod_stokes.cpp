#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

// PVFMM
#include <pvfmm_common.hpp>
#include <profile.hpp>
#include <fmm_cheb.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>
#include <cheb_node.hpp>

// LOCAL
#include <utils.hpp>

// TBSLAS
#include <utils/common.h>
#include <utils/reporter.h>
#include <utils/fields.h>
#include <kernels/mod_stokes.h>

int NUM_TIME_STEPS = 1;
double TBSLAS_DT;
double TBSLAS_MOD_STOKES_DIFF_COEFF;
double TBSLAS_MOD_STOKES_ALPHA;

// current simulation time
double tcurr = 25;

///////////////////////////////////////////////////////////////////////////////
// Test3: Stokes problem, Smooth Gaussian, FreeSpace Boundary
///////////////////////////////////////////////////////////////////////////////
template <class Real_t>
void fn_input_t3(const Real_t* coord, int n, Real_t* out){ //Input function
  int dof=3;
  Real_t L=125;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
      out[i*dof+0]=                                        0+2*L*exp(-L*r_2)*(c[0]-0.5);
      out[i*dof+1]= 4*L*L*(c[2]-0.5)*(5-2*L*r_2)*exp(-L*r_2)+2*L*exp(-L*r_2)*(c[1]-0.5);
      out[i*dof+2]=-4*L*L*(c[1]-0.5)*(5-2*L*r_2)*exp(-L*r_2)+2*L*exp(-L*r_2)*(c[2]-0.5);
    }
  }
}
template <class Real_t>
void fn_poten_t3(const Real_t* coord, int n, Real_t* out){ //Output potential
  int dof=3;
  Real_t L=125;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
      out[i*dof+0]= 0;
      out[i*dof+1]= 2*L*(c[2]-0.5)*exp(-L*r_2);
      out[i*dof+2]=-2*L*(c[1]-0.5)*exp(-L*r_2);
    }
  }
}

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

template <class Real_t>
void RunModStokesSolver(int test_case, size_t N, size_t M, bool unif, int mult_order,
                  int cheb_deg, int depth, bool adap, Real_t tol, MPI_Comm comm) {
  typedef pvfmm::FMM_Node<pvfmm::Cheb_Node<Real_t> > FMMNode_t;
  typedef pvfmm::FMM_Cheb<FMMNode_t> FMM_Mat_t;
  typedef pvfmm::FMM_Tree<FMM_Mat_t> FMM_Tree_t;
  tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
  // Find out my identity in the default communicator
  int myrank, p;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm,&p);
  // Find out number of OMP thereads.
  int omp_p=omp_get_max_threads();

  void (*fn_input_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_poten_)(const Real_t* , int , Real_t*)=NULL;
  void (*fn_grad_ )(const Real_t* , int , Real_t*)=NULL;
  const pvfmm::Kernel<Real_t>* mykernel=NULL;
  const pvfmm::Kernel<Real_t>* mykernel_grad=NULL;;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<double> modified_stokes_kernel_d =
    pvfmm::BuildKernel<double, tbslas::modified_stokes_vel>
    (tbslas::GetModfiedStokesKernelName<double>(TBSLAS_MOD_STOKES_ALPHA, TBSLAS_MOD_STOKES_DIFF_COEFF), 3, std::pair<int,int>(3,3));

      // pvfmm::BuildKernel<double, tbslas::modified_laplace_poten>
      // (tbslas::GetModfiedLaplaceKernelName<double>(TBSLAS_ALPHA), 3, std::pair<int,int>(1,1));

  switch (test_case) {
    case 1:
      fn_input_ = fn_input_t3<Real_t>;
      fn_poten_ = fn_poten_t3<Real_t>;
      mykernel  = &modified_stokes_kernel_d;
      bndry = pvfmm::FreeSpace;
      break;
    default:
      fn_input_=NULL;
      fn_poten_=NULL;
      fn_grad_ =NULL;
      break;
  }
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
  FMM_Tree_t* tree = new FMM_Tree_t(comm);
  tree->Initialize(&tree_data);
  tree->RefineTree();
  tree->InitFMM_Tree(false,bndry);

  // ======================================================================
  // SETUP FMM
  // ======================================================================
  FMM_Mat_t* fmm_mat=NULL;
  {
    fmm_mat = new FMM_Mat_t;
    fmm_mat->Initialize(mult_order,
                        tree_data.cheb_deg,
                        comm,
                        mykernel);
  }

  // ======================================================================
  // Write2File
  // ======================================================================
  if (sim_config->vtk_save_rate) {
    tree->Write2File(tbslas::GetVTKFileName(0, sim_config->vtk_filename_variable).c_str(),
                     tree_data.cheb_deg);
  }
  double in_al2,in_rl2,in_ali,in_rli;
  CheckChebOutput<FMM_Tree_t>(tree,
                              fn_input_,
                              mykernel->ker_dim[1],
                              in_al2, in_rl2, in_ali, in_rli,
                              std::string("Input"));

  // ======================================================================
  // DIFFUSUION SOLVER
  // ======================================================================
  for (int timestep = 1; timestep < NUM_TIME_STEPS+1; timestep++) {
    pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);

    pvfmm::Profile::Tic("FMM",&comm,true);
    tree->InitFMM_Tree(false,bndry);
    tree->SetupFMM(fmm_mat);
    tree->RunFMM();
    tree->Copy_FMMOutput(); //Copy FMM output to tree Data.
    pvfmm::Profile::Toc();

    pvfmm::Profile::Toc();        // solve

    tcurr += TBSLAS_DT;

    if (sim_config->vtk_save_rate) {
      if (timestep % sim_config->vtk_save_rate == 0)
        tree->Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(),
                         tree_data.cheb_deg);
    }
  }

  // =========================================================================
  // REPORT RESULTS
  // =========================================================================
  //Find error in FMM output.
  double al2,rl2,ali,rli;
  CheckChebOutput<FMM_Tree_t>(tree,
                              fn_poten_,
                              mykernel->ker_dim[1],
                              al2,rl2,ali,rli,
                              std::string("Output"));
  int num_leaves = tbslas::CountNumLeafNodes(*tree);
  typedef tbslas::Reporter<Real_t> Rep;
  if(!myrank) {
    Rep::AddData("TOL", sim_config->tree_tolerance);
    Rep::AddData("NOCT", num_leaves, tbslas::REP_INT);

    Rep::AddData("DT", sim_config->dt);
    Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);

    Rep::AddData("DIFF",  TBSLAS_MOD_STOKES_DIFF_COEFF);
    Rep::AddData("ALPHA", TBSLAS_MOD_STOKES_ALPHA);

    Rep::AddData("InAL2", in_al2);
    Rep::AddData("OutAL2", al2);

    Rep::AddData("InRL2", in_rl2);
    Rep::AddData("OutRL2", rl2);

    Rep::AddData("InALINF", in_ali);
    Rep::AddData("OutALINF", ali);

    Rep::AddData("InRLINF", in_rli);
    Rep::AddData("OutRLINF", rli);

    Rep::Report();
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
  sim_config->vtk_filename_variable   = "conc";

  NUM_TIME_STEPS    = sim_config->total_num_timestep;
  TBSLAS_DT         = sim_config->dt;
  TBSLAS_MOD_STOKES_DIFF_COEFF = 1;//0.0001;
  TBSLAS_MOD_STOKES_ALPHA      = 1;//(1.0)/TBSLAS_DT/TBSLAS_DIFF_COEFF;

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
  pvfmm::Profile::Tic("RunModStokesSolver",&comm,true);
  RunModStokesSolver<double>(test,
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
