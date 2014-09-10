// Copyright (C) 2014 by Arash Bakhtiari

#ifndef TRAJ_H
#define TRAJ_H

#include <VecField.hpp>
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
         const int n,
         std::vector<real_t>& xsol);


}  // namespace tbslas

#include "Traj.ipp"

#endif /* TRAJ_H */
