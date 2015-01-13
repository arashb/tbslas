#include <mpi.h>
#include <pvfmm_common.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <string>
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

int NUM_TIME_STEPS = 1;
const double TBSLAS_DT = 0.01;
const double TBSLAS_DIFF_COEFF = 0.01;
const double TBSLAS_ALPHA = (1.0)/TBSLAS_DT/TBSLAS_DIFF_COEFF;

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
  // Real_t a=-160;
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
  // Real_t a=-160;
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
void fmm_test(int test_case, size_t N, size_t M, bool unif, int mult_order,
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

  switch (test_case) {
    case 1:
      fn_input_ = fn_input_t1<Real_t>;
      fn_poten_ = fn_poten_t1<Real_t>;
      mykernel  = &tbslas::ModifiedLaplaceKernel<Real_t>::potn_ker();
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

  // **********************************************************************
  // TREE DATA
  // **********************************************************************
  //Various parameters.
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

  // **********************************************************************
  // SETUP TREE
  // **********************************************************************
  pvfmm::Profile::Tic("TreeSetup",&comm,true,1);
  {
    FMM_Tree_t* tree=new FMM_Tree_t(comm);
    tree->Initialize(&tree_data);
    tree->InitFMM_Tree(adap,bndry); //Adaptive refinement.

    pt_coord.clear();
    FMMNode_t* node=static_cast<FMMNode_t*>(tree->PreorderFirst());
    while (node!=NULL) {
      if (node->IsLeaf() && !node->IsGhost()) {
        Real_t* c=node->Coord();
        Real_t s=pow(0.5,node->Depth()+1);
        pt_coord.push_back(c[0]+s);
        pt_coord.push_back(c[1]+s);
        pt_coord.push_back(c[2]+s);
      }
      node=static_cast<FMMNode_t*>(tree->PreorderNxt(node));
    }
    delete tree;
    tree_data.pt_coord=pt_coord;
    tree_data.max_pts=1; // Points per octant.
  }
  pvfmm::Profile::Toc();

  //Create Tree and initialize with input data.
  FMM_Tree_t* tree=new FMM_Tree_t(comm);
  tree->Initialize(&tree_data);

  // **********************************************************************
  // SETUP FMM
  // **********************************************************************
  //Initialize FMM_Mat.
  FMM_Mat_t* fmm_mat=NULL;
  {
    fmm_mat=new FMM_Mat_t;
    fmm_mat->Initialize(mult_order,
                        tree_data.cheb_deg,
                        comm,
                        mykernel);
  }

  tree->InitFMM_Tree(false,bndry);
  //Find error in Chebyshev approximation.
  CheckChebOutput<FMM_Tree_t>
      (tree, (typename TestFn<Real_t>::Fn_t) fn_input_,
       mykernel->ker_dim[0], std::string("Input"));

  char out_name_buffer[300];
  snprintf(out_name_buffer,
           sizeof(out_name_buffer),
           OUTPUT_FILE_FORMAT,
           tbslas::get_result_dir().c_str(),
           OUTPUT_FILE_PREFIX,
           "diff",
           0);

  //Write2File
  tree->Write2File(out_name_buffer,tree_data.cheb_deg);

  for (int ts_counter = 1; ts_counter < NUM_TIME_STEPS+1; ts_counter++) {
    // Setup FMM
    tree->SetupFMM(fmm_mat);
    tree->RunFMM();
    tree->Copy_FMMOutput(); //Copy FMM output to tree Data.

    //Find error in FMM output.
    if (ts_counter == 1)
    CheckChebOutput<FMM_Tree_t>(tree,
                                (typename TestFn<Real_t>::Fn_t) fn_poten_,
                                mykernel->ker_dim[1], std::string("Output"));

    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             OUTPUT_FILE_FORMAT,
             tbslas::get_result_dir().c_str(),
             OUTPUT_FILE_PREFIX,
             "diff",
             ts_counter);

    //Write2File
    tree->Write2File(out_name_buffer,tree_data.cheb_deg);

    //Initialize FMM Tree for the next time step
    tree->InitFMM_Tree(false,bndry);
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
  int p;
  MPI_Comm_size(comm,&p);
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  if (p>8) { // Remove slow processors.
    MPI_Comm comm_=MPI_COMM_WORLD;
    size_t N=2048;
    pvfmm::Matrix<double> A(N,N);
    pvfmm::Matrix<double> B(N,N);
    pvfmm::Matrix<double> C(N,N);
    for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {
      A[i][j]=i+j;
      B[i][j]=i-j;
    }
    C=A*B;
    double t=-omp_get_wtime();
    C=A*B;
    t+=omp_get_wtime();

    double tt;
    int myrank, np;
    MPI_Comm_size(comm_,&np);
    MPI_Comm_rank(comm_,&myrank);
    MPI_Allreduce(&t, &tt, 1, pvfmm::par::Mpi_datatype<double>::value(), MPI_SUM, comm_);
    tt=tt/np;

    int clr=(t<tt*1.5?0:1);
    MPI_Comm_split(comm_, clr, myrank, &comm );
    if (clr) {
      MPI_Finalize();
      return 0;
    }
  }

  // Read command line options.
  commandline_option_start(argc, argv);
  omp_set_num_threads( atoi(commandline_option(argc, argv,  "-omp",     "1", false, "-omp  <int> = (1)    : Number of OpenMP threads."          )));
  size_t   N=(size_t)strtod(commandline_option(argc, argv,    "-N",     "1",  true, "-N    <int>          : Number of point sources."           ),NULL);
  size_t   M=(size_t)strtod(commandline_option(argc, argv,    "-M",     "1", false, "-M    <int>          : Number of points per octant."       ),NULL);
  bool  unif=              (commandline_option(argc, argv, "-unif",    NULL, false, "-unif                : Uniform point distribution."        )!=NULL);
  int      m=       strtoul(commandline_option(argc, argv,    "-m",    "10", false, "-m    <int> = (10)   : Multipole order (+ve even integer)."),NULL,10);
  int      q=       strtoul(commandline_option(argc, argv,    "-q",    "14", false, "-q    <int> = (14)   : Chebyshev order (+ve integer)."     ),NULL,10);
  int      d=       strtoul(commandline_option(argc, argv,    "-d",    "15", false, "-d    <int> = (15)   : Maximum tree depth."                ),NULL,10);
  double tol=        strtod(commandline_option(argc, argv,  "-tol",  "1e-5", false, "-tol <real> = (1e-5) : Tolerance for adaptive refinement." ),NULL);
  bool  adap=              (commandline_option(argc, argv, "-adap",    NULL, false, "-adap                : Adaptive tree refinement."          )!=NULL);
  int   test=       strtoul(commandline_option(argc, argv, "-test",     "1", false, "-test <int> = (1)    : 1) Modified Laplace, Smooth Gaussian, FreeSpace Boundary"),NULL,10);
  int     tn =         strtoul(commandline_option(argc, argv,    "-tn",    "1", false, "-tn   <int> = (1)   : Number of time steps."     ),NULL,10);

  NUM_TIME_STEPS = tn;
  commandline_option_end(argc, argv);
  pvfmm::Profile::Enable(true);

  // =========================================================================
  // PRINT METADATA
  // =========================================================================
  if (!myrank) {
    MetaData_t::Print();
  }

  // Run FMM with above options.
  pvfmm::Profile::Tic("FMM_Test",&comm,true);
  fmm_test<double>(test, N,M,unif, m,q, d, adap,tol, comm);
  pvfmm::Profile::Toc();

  //Output Profiling results.
  pvfmm::Profile::print(&comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
