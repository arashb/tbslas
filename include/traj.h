// Copyright (C) 2014 by Arash Bakhtiari

#ifndef INCLUDE_TRAJ_H_
#define INCLUDE_TRAJ_H_

#include <vector>
#include "vec_field.h"
#include "cubic_interp_policy.h"

namespace tbslas {

template<typename real_t,
         typename InterpPolicy,
         int sdim, int vdim>
void
integrate_rk2(const VecField<real_t, sdim, vdim>& field,
              const std::vector<real_t>& xinit,
              const real_t dt,
              const InterpPolicy& interp_policy,
              std::vector<real_t>& xsol);

template<typename real_t,
         typename InterpPolicy,
         int sdim, int vdim>
void
traj_rk2(const VecField<real_t, sdim, vdim>& field,
         const std::vector<real_t>& xinit,
         const real_t tinit,
         const real_t tfinal,
         const int num_rk_step,
         const InterpPolicy& interp_policy,
         std::vector<real_t>& xsol);

}  // namespace tbslas

#include "traj.inc"

#endif  // INCLUDE_TRAJ_H_
