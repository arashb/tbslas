#include <mpi.h>
#include <pvfmm_common.hpp>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <stdio.h>

#include <mpi_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>

const int DATA_DOF=3;

template <class Real_t>
void fn(const Real_t* coord, int n, Real_t* out){ //Output potential
  Real_t L=125;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-0.5)*(c[0]-0.5)+(c[1]-0.5)*(c[1]-0.5)+(c[2]-0.5)*(c[2]-0.5);
      out[i*DATA_DOF+0]= 0;
      out[i*DATA_DOF+1]= (c[2]-0.5)*exp(-L*r_2);
      out[i*DATA_DOF+2]=-(c[1]-0.5)*exp(-L*r_2);
    }
  }
}

template <class Real_t>
void tree_test(size_t N, size_t M, int cheb_deg, int depth, bool adap, Real_t tol, MPI_Comm comm){
  typedef pvfmm::Cheb_Node<Real_t> Node_t;
  typedef pvfmm::MPI_Tree<Node_t> Tree_t;

  //Various parameters.
  typename Node_t::NodeData tree_data;
  tree_data.dim=COORD_DIM;
  tree_data.max_depth=depth;
  tree_data.cheb_deg=cheb_deg;

  //Set input function pointer
  tree_data.input_fn=fn<Real_t>;
  tree_data.data_dof=DATA_DOF;
  tree_data.tol=tol;

  //Set source coordinates.
  std::vector<Real_t> pt_coord;
  pt_coord=point_distrib<Real_t>(UnifGrid,N,comm);
  tree_data.max_pts=M; // Points per octant.
  tree_data.pt_coord=pt_coord;

  //Create Tree and initialize with input data.
  Tree_t tree(comm);
  tree.Initialize(&tree_data);
  tree.RefineTree();
  //tree.Balance21(pvfmm::FreeSpace);
  //tree.RedistNodes();

  //Write2File
  tree.Write2File("result/output",4);

  //Find error in Chebyshev approximation.
  CheckChebOutput<Tree_t>(&tree, (typename TestFn<Real_t>::Fn_t) &fn<Real_t>, DATA_DOF, "Input");

  //Check Tree.
  #ifndef NDEBUG
  pvfmm::Profile::Tic("CheckTree",&comm,true,1);
  tree.CheckTree();
  pvfmm::Profile::Toc();
  #endif
}

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // Read command line options.
  commandline_option_start(argc, argv);
  omp_set_num_threads( atoi(commandline_option(argc, argv,  "-omp",     "1", false, "-omp  <int> = (1)    : Number of OpenMP threads."          )));
  size_t   N=(size_t)strtod(commandline_option(argc, argv,    "-N",     "1",  true, "-N    <int>          : Number of point sources."           ),NULL);
  size_t   M=(size_t)strtod(commandline_option(argc, argv,    "-M",     "1", false, "-M    <int>          : Number of points per octant."       ),NULL);
  int      q=       strtoul(commandline_option(argc, argv,    "-q",    "14", false, "-q    <int> = (14)   : Chebyshev order (+ve integer)."     ),NULL,10);
  int      d=       strtoul(commandline_option(argc, argv,    "-d",    "15", false, "-d    <int> = (15)   : Maximum tree depth."                ),NULL,10);
  double tol=        strtod(commandline_option(argc, argv,  "-tol",  "1e-5", false, "-tol <real> = (1e-5) : Tolerance for adaptive refinement." ),NULL);
  bool  adap=              (commandline_option(argc, argv, "-adap",    NULL, false, "-adap                : Adaptive tree refinement."          )!=NULL);
  commandline_option_end(argc, argv);

  tree_test<double>(N, M, q, d, adap, tol, comm);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}

