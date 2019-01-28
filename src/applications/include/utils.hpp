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
#include <tree/tree_utils.h>

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

template <class Real_t>
std::vector<Real_t> point_distrib(DistribType dist_type, size_t N, MPI_Comm comm){
  int np, myrank;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &myrank);
  static size_t seed=myrank+1; seed+=np;
  srand48(seed);

  std::vector<Real_t> coord;
  switch(dist_type){
  case UnifGrid:
    {
      size_t NN=(size_t)round(pow((double)N,1.0/3.0));
      size_t N_total=NN*NN*NN;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;
      for(size_t i=start;i<end;i++){
        coord.push_back(((Real_t)((i/  1    )%NN)+0.5)/NN);
        coord.push_back(((Real_t)((i/ NN    )%NN)+0.5)/NN);
        coord.push_back(((Real_t)((i/(NN*NN))%NN)+0.5)/NN);
      }
    }
    break;
  case RandUnif:
    {
      size_t N_total=N;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;
      size_t N_local=end-start;
      coord.resize(N_local*3);

      for(size_t i=0;i<N_local*3;i++)
        coord[i]=((Real_t)drand48());
    }
    break;
  case RandGaus:
    {
      Real_t e=2.7182818284590452;
      Real_t log_e=log(e);
      size_t N_total=N;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;

      for(size_t i=start*3;i<end*3;i++){
        Real_t y=-1;
        while(y<=0 || y>=1){
          Real_t r1=sqrt(-2*log(drand48())/log_e)*cos(2*M_PI*drand48());
          Real_t r2=pow(0.5,i*10/N_total);
          y=0.5+r1*r2;
        }
        coord.push_back(y);
      }
    }
    break;
  case RandElps:
    {
      size_t N_total=N;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;
      size_t N_local=end-start;
      coord.resize(N_local*3);

      const Real_t r=0.45;
      const Real_t center[3]={0.5,0.5,0.5};
      for(size_t i=0;i<N_local;i++){
        Real_t* y=&coord[i*3];
        Real_t phi=2*M_PI*drand48();
        Real_t theta=M_PI*drand48();
        y[0]=center[0]+0.25*r*sin(theta)*cos(phi);
        y[1]=center[1]+0.25*r*sin(theta)*sin(phi);
        y[2]=center[2]+r*cos(theta);
      }
    }
    break;
  case RandSphr:
    {
      size_t N_total=N;
      size_t start= myrank   *N_total/np;
      size_t end  =(myrank+1)*N_total/np;
      size_t N_local=end-start;
      coord.resize(N_local*3);

      const Real_t center[3]={0.5,0.5,0.5};
      for(size_t i=0;i<N_local;i++){
        Real_t* y=&coord[i*3];
        Real_t r=1;
        while(r>0.5 || r==0){
          y[0]=drand48(); y[1]=drand48(); y[2]=drand48();
          r=sqrt((y[0]-center[0])*(y[0]-center[0])
                +(y[1]-center[1])*(y[1]-center[1])
                +(y[2]-center[2])*(y[2]-center[2]));
          y[0]=center[0]+0.45*(y[0]-center[0])/r;
          y[1]=center[1]+0.45*(y[1]-center[1])/r;
          y[2]=center[2]+0.45*(y[2]-center[2])/r;
        }
      }
    }
    break;
  default:
    break;
  }
  return coord;
}

#include <utils.txx>

#endif
