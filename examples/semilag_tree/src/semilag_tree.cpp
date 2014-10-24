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
#include <tree/advect_tree_semilag.h>
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
    // INIT THE TREES
    tbslas::Tree_t<double> tvel_curr(comm);
    tbslas::construct_tree<double>(N, M, q, d, adap, tol, comm,
                                   tbslas::get_vorticity_field<double,3>,
                                   3,
                                   tvel_curr);
    tvel_curr.ConstructLET(pvfmm::FreeSpace);
    tvel_curr.Write2File("result/output_vel_00_", q);

    tbslas::Tree_t<double>* tconc_curr = new tbslas::Tree_t<double>(comm);
    tbslas::construct_tree<double>(N, M, q, d, adap, tol, comm,
                                   tbslas::get_gaussian_field<double,3>,
                                   1,
                                   *tconc_curr);
    // tbslas::clone_tree<double,
    //                    tbslas::Node_t<double>,
    //                    tbslas::Tree_t<double> > (tvel_curr, *tconc_curr, 1);
    // tbslas::init_tree<double,
    //                   tbslas::Node_t<double>,
    //                   tbslas::Tree_t<double> > (*tconc_curr,
    //                                             tbslas::get_gaussian_field<double,3>,
    //                                             1);
    char out_name_buffer[50];
    snprintf(out_name_buffer, sizeof(out_name_buffer), "result/output_%d_", 0);
    tconc_curr->Write2File(out_name_buffer, q);

    // clone a tree
    tbslas::Tree_t<double>* tconc_next = new tbslas::Tree_t<double>(comm);
    tbslas::clone_tree<double,
                       tbslas::Node_t<double>,
                       tbslas::Tree_t<double> >(*tconc_curr, *tconc_next, 1);

    // simulation info
    int tstep       = 1;
    double dt       = 0.5;
    int num_rk_step = 1;
    int tn          = 1;

    // TIME STEPPING
    for (int tstep = 1; tstep < tn+1; tstep++) {
      printf("====================\n");
      printf("TIME STEP: %d\n", tstep);

      tconc_curr->ConstructLET(pvfmm::FreeSpace);
      tbslas::advect_tree_semilag<double,
                                  tbslas::Node_t<double>,
                                  tbslas::Tree_t<double> >(tvel_curr,
                                                           *tconc_curr,
                                                           *tconc_next,
                                                           tstep,
                                                           dt,
                                                           num_rk_step);

      snprintf(out_name_buffer, sizeof(out_name_buffer), "result/output_%d_", tstep);
      tconc_next->Write2File(out_name_buffer,q);

      tbslas::swap_trees_pointers(&tconc_curr, &tconc_next);
    }
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
