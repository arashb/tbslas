/**
 * \file utils.hpp
 * \author Dhairya Malhotra, dhairya.malhotra@gmail.com
 * \date 1-1-2011
 */

#ifndef _UTILS_
#define _UTILS_

#include <vector>
#include <string>

#include <mpi.h>
#include <cheb_utils.hpp>
#include <fmm_tree.hpp>

#include <utils/common.h>
#include <tree/utils_tree.h>

template <class FMM_Mat_t>
void CheckFMMOutput(pvfmm::FMM_Tree<FMM_Mat_t>* mytree, const pvfmm::Kernel<typename FMM_Mat_t::Real_t>* mykernel);

template <class FMMTree_t,
          class Functor_t>
void CheckChebOutput(FMMTree_t* mytree,
                     Functor_t fn_poten,
                     int fn_dof,
                     typename FMMTree_t::Real_t& al2,
                     typename FMMTree_t::Real_t& rl2,
                     typename FMMTree_t::Real_t& ali,
                     typename FMMTree_t::Real_t& rli,
                     std::string t_name="");

enum DistribType{
  UnifGrid,
  RandUnif,
  RandGaus,
  RandElps,
  RandSphr
};

void commandline_option_start(int argc, char** argv, const char* help_text=NULL);

const char* commandline_option(int argc, char** argv, const char* opt, const char* def_val, bool required, const char* err_msg);

void commandline_option_end(int argc, char** argv);

void parse_command_line_options(int argc, char** argv);

#include <utils.txx>

#endif
