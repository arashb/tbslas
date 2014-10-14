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

#include <tree/tree_common.h>
#include <tree/node_field_functor.h>
#include <tree/semilag_advect_tree.h>
#include <tree/tree_utils.h>

#include <semilag/utils.h>

int main (int argc, char **argv) {
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

  {
    tbslas::Tree_t<double> tvel_curr(comm);
    tbslas::semilag_construct_tree<double>(N, M, q, d, adap, tol, comm,
                                           slas::get_vorticity_field<double,3>,
                                           3,
                                           tvel_curr);
    tvel_curr.Write2File("result/output_vel_00_", 4);

    tbslas::Tree_t<double> tconc_curr(comm);
    tbslas::semilag_construct_tree<double>(N, M, q, d, adap, tol, comm,
                                           slas::get_gaussian_field<double,3>,
                                           1,
                                           tconc_curr);
    // clone a tree
    tbslas::Tree_t<double> tconc_next(comm);
    // TODO: clone the next tree from the current tree (NOT FROM ANALYTICAL FUNCTION)
    // Tree_t tconc_next(tconc_curr);
    tbslas::semilag_construct_tree<double>(N, M, q, d, adap, tol, comm,
                                           slas::get_gaussian_field<double,3>,
                                           1,
                                           tconc_next);

    tconc_curr.Write2File("result/output_00_",4);

    tbslas::semilag_advect_tree<double>(tvel_curr, tconc_curr, tconc_next);
    tconc_next.Write2File("result/output_01_",4);
    //Find error in Chebyshev approximation.
    // CheckChebOutput<Tree_t>(&tree,
    //                         (typename TestFn<double>::Fn_t) &slas::get_gaussian_field<double,3>,
    //                         DATA_DOF,
    //                         "Input");
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
