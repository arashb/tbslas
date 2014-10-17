// Copyright (C) 2014 by Arash Bakhtiari
#include <vector>
#include <cmath>

#include <gtest/gtest.h>

#include "semilag/utils.h"
#include "semilag/cubic_interp_policy.h"
#include "semilag/vec_field.h"
#include "semilag/semilag.h"

typedef tbslas::VecField<double,3,3> VFieldD33;
typedef tbslas::VecField<double,3,1> VFieldD31;
typedef std::vector<double> VecD              ;

#define PI 3.14159

class SemilagTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    domain_res      = 16;
    num_points_dim  = domain_res+1;
    cfl             = 1;
    dx              = 1.0/domain_res;
    dt              = cfl*dx;

    tf              = 2*PI;
    tn              = tf/dt;

    // tn              = 1;
    // tf              = tn*dt;

    num_rk_step     = 1;

    sdim            = 3;
    num_points      = (size_t)std::pow(num_points_dim, sdim);
    num_points_size = num_points*sdim;

  }

  virtual void TearDown() {
  }

  size_t domain_res      ;
  size_t num_points_dim  ;
  int cfl                ;
  double dx              ;
  double dt              ;
  double tf              ;
  int tn                 ;
  // int tn                 ;
  // double tf              ;
  int num_rk_step        ;
  int sdim               ;
  size_t num_points      ;
  size_t num_points_size ;

};

TEST_F(SemilagTest, SteadyVorticityWithFunctor) {
  VecD points_val(num_points);  // one dimensional concentration values

  // init position of grid points
  VecD points_pos(num_points_size);
  tbslas::get_reg_grid_points<double,3>(num_points_dim, points_pos.data());

  // init velocity values
  VecD points_val_vorticity(num_points_size);
  tbslas::get_vorticity_field<double,3>(points_pos.data(), num_points,
                                        points_val_vorticity.data());

  VecD points_val_gaussian(num_points);
  tbslas::get_gaussian_field<double,3>(points_pos.data(), num_points,
                                       points_val_gaussian.data());

  VFieldD33 vel_field;
  VFieldD31 con_field;

  vel_field.init(points_pos, points_val_vorticity);
  con_field.init(points_pos, points_val_gaussian);

  for(int timestep = 1; timestep <= tn; timestep++) {
    tbslas::get_vorticity_field<double,3>(points_pos.data(),
                                          num_points,
                                          points_val_vorticity.data()
                                          );
    vel_field.push_back_values(points_val_vorticity,
                               timestep*dt);
  }

  for(int timestep = 1; timestep <= tn; timestep++) {
    std::cout << "****************************************" << std::endl;
    std::cout << "-> timestep: " << timestep
              << " time: " << dt*timestep << std::endl;

    tbslas::semilag_rk2(tbslas::get_vorticity_field<double,3>,
                        tbslas::get_gaussian_field<double,3>,
                        points_pos,
                        sdim,
                        timestep,
                        dt,
                        num_rk_step,
                        points_val
                        );

    // vel_field.push_back_values(points_val_vorticity, dt*timestep);
    con_field.push_back_values(points_val, dt*timestep);
  }

  vel_field.save("semilag_functor_vel_");
  con_field.save("semilag_functor_con_");

  // VFieldD33 vel_field_expected;
  // VFieldD31 con_field_expected;

  // for (int i = 0; i < num_points; i++) {
  //   std::cout << "expected: " << expected_val[i] << " computed: " << points_val[1] << std::endl;
  //   ASSERT_NEAR(expected_val[i], points_val[i], 0.001);
  // }
}

TEST_F(SemilagTest, SteadyVorticityWithVecField) {

}
