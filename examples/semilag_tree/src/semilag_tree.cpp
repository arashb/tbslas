#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>

#include <pvfmm_common.hpp>
#include <mpi_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>
#include <vector.hpp>
#include <cheb_utils.hpp>

#include <utils/common.h>
// Enable profiling
#define __TBSLAS_PROFILE__ 5
#include <utils/profile.h>
#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // Read command line options.
  commandline_option_start(argc, argv);
  omp_set_num_threads( atoi(commandline_option(argc, argv,  "-omp",     "1", false, "-omp  <int> = (1)    : Number of OpenMP threads."          )));
  size_t   N = (size_t)strtod(commandline_option(argc, argv,    "-N",     "1",  true, "-N    <int>          : Number of point sources."           ),NULL);
  size_t   M = (size_t)strtod(commandline_option(argc, argv,    "-M",     "1", false, "-M    <int>          : Number of points per octant."       ),NULL);
  int      q =         strtoul(commandline_option(argc, argv,    "-q",    "14", false, "-q    <int> = (14)   : Chebyshev order (+ve integer)."     ),NULL,10);
  int      d =         strtoul(commandline_option(argc, argv,    "-d",    "15", false, "-d    <int> = (15)   : Maximum tree depth."                ),NULL,10);
  double tol =         strtod(commandline_option(argc, argv,  "-tol",  "1e-5", false, "-tol <real> = (1e-5) : Tolerance for adaptive refinement." ),NULL);
  bool  adap =               (commandline_option(argc, argv, "-adap",    NULL, false, "-adap                : Adaptive tree refinement."          )!=NULL);
  int     tn =         strtoul(commandline_option(argc, argv,    "-tn",    "1", false, "-tn    <int> = (1)   : Number of time steps."     ),NULL,10);
  commandline_option_end(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    tbslas::Profile<double>::Enable(true, &comm);

    // INIT THE TREES
    Tree_t tvel_curr(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                   tbslas::get_vorticity_field<double,3>,
                                   3,
                                   tvel_curr);
    tvel_curr.ConstructLET(pvfmm::FreeSpace);
    char out_name_buffer[300];

    snprintf(out_name_buffer, sizeof(out_name_buffer),
             "%s/sltree_vel_%d_", tbslas::get_result_dir().c_str(), 0);
    tvel_curr.Write2File(out_name_buffer, q);

    double vel_max_value;
    int vel_max_depth;
    tbslas::GetMaxTreeValues<Tree_t>
        (tvel_curr, vel_max_value, vel_max_depth);

    if (!myrank)
      printf("%d: VEL MAX VALUE: %f VEL MAX DEPTH:%d\n",
             myrank,
             vel_max_value,
             vel_max_depth);

    Tree_t* tconc_curr = new Tree_t(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                   tbslas::get_gaussian_field<double,3>,
                                   1,
                                   *tconc_curr);

    snprintf(out_name_buffer, sizeof(out_name_buffer),
             "%s/sltree_val_%d_", tbslas::get_result_dir().c_str(), 0);
    tconc_curr->Write2File(out_name_buffer, q);

    double conc_max_value;
    int conc_max_depth;
    tbslas::GetMaxTreeValues<Tree_t>
        (*tconc_curr, conc_max_value, conc_max_depth);

    if (!myrank)
      printf("%d:CON MAX VALUE: %f CON MAX DEPTH:%d\n",
             myrank,
             conc_max_value,
             conc_max_depth);

    // clone a tree
    Tree_t* tconc_next = new Tree_t(comm);
    tbslas::CloneTree<Tree_t>(*tconc_curr, *tconc_next, 1);

    // set the input_fn to NULL -> need for adaptive refinement
    std::vector<Node_t*>  ncurr_list = tconc_curr->GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = NULL;
    }

    std::vector<Node_t*>  nnext_list = tconc_next->GetNodeList();
    for(int i = 0; i < nnext_list.size(); i++) {
      nnext_list[i]->input_fn = NULL;
    }

    // simulation parameters
    double cfl      = 1;
    double dx_min   = pow(0.5, conc_max_depth);
    double dt       = (cfl * dx_min)/vel_max_value;
    int num_rk_step = 1;

    // TIME STEPPING
    for (int tstep = 1; tstep < tn+1; tstep++) {
      if(!myrank) {
        printf("============================\n");
        printf("dt: %f tstep: %d\n", dt, tstep);
        printf("============================\n");
      }
      tbslas::SolveSemilagTree<Tree_t>(tvel_curr,
                                       *tconc_curr,
                                       *tconc_next,
                                       tstep,
                                       dt,
                                       num_rk_step);

      // refine the tree according to the computed values
      tbslas::Profile<double>::Tic("RefineTree",false,5);
      tconc_next->RefineTree();
      tbslas::Profile<double>::Toc();

      snprintf(out_name_buffer, sizeof(out_name_buffer),
               "%s/sltree_val_%d_", tbslas::get_result_dir().c_str(), tstep);
      tconc_next->Write2File(out_name_buffer,q);

      // prepare the next step tree
      tbslas::SyncTreeRefinement(*tconc_next, *tconc_curr);

      tbslas::swap_pointers(&tconc_curr, &tconc_next);
    }

    //Output Profiling results.
    tbslas::Profile<double>::print(&comm);

    // CLEAN UP MEM.
    delete tconc_curr;
    delete tconc_next;
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
