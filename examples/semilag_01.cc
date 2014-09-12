// Copyright (C) 2014 by Arash Bakhtiari

#include <utils.h>
#include <vec_field.h>
#include <traj.h>
#include <vector>

typedef tbslas::VecField<double,3,3> VFieldD33;
typedef tbslas::VecField<double,3,1> VFieldD31;
typedef std::vector<double> VecD;

int main(int argc, char *argv[]) {
  VFieldD33 vec_field_velocity;
  VFieldD31 vec_field_concentration;
  size_t dN;

  dN = 65;
  size_t tN = dN*dN*dN;
  VecD pnts_pos           = tbslas::gen_reg_grid_points<double>(dN);
  VecD pnts_vls_vorticity = tbslas::vorticity_field(pnts_pos);
  VecD pnts_vls_gaussian  = tbslas::gaussian_field(pnts_pos);

  vec_field_velocity.init(pnts_pos, pnts_vls_vorticity);
  vec_field_concentration.init(pnts_pos, pnts_vls_gaussian);

  vec_field_velocity.write2file("semilag_vel_");
  vec_field_concentration.write2file("semilag_con_");

  return 0;
}
