/**
 * \file utils.hpp
 * \author Dhairya Malhotra, dhairya.malhotra@gmail.com
 * \date 1-1-2011
 */

#ifndef _UTILS_
#define _UTILS_

#include <vector>
#include <mpi.h>
#include <cheb_utils.hpp>
#include <fmm_tree.hpp>

#include <utils/common.h>

extern double tcurr;
template<typename real_t, int sdim>
void
get_gaussian_field_cylinder_atT(const real_t* points_pos,
                                int num_points,
                                real_t* out) {
  real_t xc      = 0.6;
  real_t yc      = 0.5;
  real_t r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
  xc = 0.5+r*cos(tcurr);
  yc = 0.5+r*sin(tcurr);
  const real_t theta   = 0.0;
  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t A       = 1.0;

  tbslas::get_gaussian_field_cylinder<real_t, sdim>(points_pos,
                                                    num_points,
                                                    out,
                                                    xc,
                                                    yc,
                                                    theta,
                                                    sigma_x,
                                                    sigma_y,
                                                    A);
};

template <class FMM_Mat_t>
void CheckFMMOutput(pvfmm::FMM_Tree<FMM_Mat_t>* mytree, const pvfmm::Kernel<typename FMM_Mat_t::Real_t>* mykernel);

template <class Real_t>
struct TestFn{
  typedef void (*Fn_t)(Real_t* c, int n, Real_t* out);
};

template <class FMMTree_t>
void CheckChebOutput(FMMTree_t* mytree, typename TestFn<typename FMMTree_t::Real_t>::Fn_t fn_poten, int fn_dof, std::string t_name="");

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

#include <utils.txx>

#endif
