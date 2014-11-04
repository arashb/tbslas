/**
 * \file profile.hpp
 * \author Dhairya Malhotra, dhairya.malhotra@gmail.com
 * \date 2-11-2011
 * \brief This file contains definition of the class Profile.
 */

#ifndef SRC_PROFILE_H_
#define SRC_PROFILE_H_

#include <mpi.h>
#include <string>
#include <vector>
#include <stack>

//#include <pvfmm_common.hpp>

#ifndef __PROFILE__
#define __PROFILE__ -1
#endif

namespace tbslas {

template<typename T>
class Profile{
 public:

  static void Add_FLOP(long long inc);

  static void Add_MEM(long long inc);

  static void Enable(bool state) {enable_state=state;}

  static void Enable(bool state, MPI_Comm* mpi_comm_) {
    enable_state=state;
    mpi_comm = mpi_comm_;
  }

  static void Tic(const char* name_, const MPI_Comm* comm_=NULL,bool sync_=false, int level=0);

  static void Tic(const char* name_, bool sync_=false, int level=0) {
    tbslas::Profile<T>::Tic(name_, mpi_comm, sync_, level);
  };

  static void Toc();

  static void print(const MPI_Comm* comm_=NULL);

  static void reset();

 private:
  static MPI_Comm* mpi_comm;
  static long long FLOP;
  static long long MEM;
  static bool enable_state;
  static std::stack<bool> sync;
  static std::stack<std::string> name;
  static std::stack<MPI_Comm*> comm;
  static std::vector<long long> max_mem;

  static unsigned int enable_depth;
  static std::stack<int> verb_level;

  static std::vector<bool> e_log;
  static std::vector<bool> s_log;
  static std::vector<std::string> n_log;
  static std::vector<double> t_log;
  static std::vector<long long> f_log;
  static std::vector<long long> m_log;
  static std::vector<long long> max_m_log;
};

}  // namespace tbslas

#include "profile.inc"

#endif  // SRC_PROFILE_H_
