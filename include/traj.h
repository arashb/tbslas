// Copyright (C) 2014 by Arash Bakhtiari

#ifndef INCLUDE_TRAJ_H_
#define INCLUDE_TRAJ_H_

#include <vec_field.h>
#include <vector>

namespace tbslas {

template<typename real_t, int sdim, int vdim>
std::vector<real_t>
integrate_rk2(VecField<real_t, sdim, vdim>& field,
              std::vector<real_t>& init_val,
              real_t dt);

template<typename real_t, int sdim, int vdim>
void
traj_rk2(VecField<real_t, sdim, vdim>& field,
         std::vector<real_t>& xinit,
         real_t tinit,
         real_t tfinal,
         const int num_rk_step,
         std::vector<real_t>& xsol);

}  // namespace tbslas

#include "traj.inc"

#endif  // INCLUDE_TRAJ_H_
