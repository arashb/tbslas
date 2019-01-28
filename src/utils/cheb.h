#include <cheb_utils.hpp>

#ifndef _TBSLAS_CHEB_H_
#define _TBSLAS_CHEB_H_

#ifndef __USE_SPARSE_GRID__
#define __USE_SPARSE_GRID__ 0
#endif

namespace tbslas {

template <class Real>
static void GramSchmidt(pvfmm::Matrix<Real> &M0, pvfmm::Vector<int> &flag,
                        long pivot, std::vector<Real> &q) {
  Real eps = 1.0;
  while (eps + (Real)1.0 > 1.0) eps *= 0.5;
  eps = sqrt(eps);

  Real norm = 0;
  {                                         // set q, norm=|q|_2
    for (long i = 0; i < M0.Dim(1); i++) {  // set q
      q[i] = M0[pivot][i];
      norm += q[i] * q[i];
    }
    norm = sqrt(norm);
    Real norm_inv = 1.0 / norm;
    for (long i = 0; i < M0.Dim(1); i++) {  // normalize q
      q[i] *= norm_inv;
    }
  }
  flag[pivot] = 1;
  if (norm < eps) return;
#pragma omp parallel for schedule(static)
  for (long i = 0; i < M0.Dim(0); i++) {  // compute M0=M0-M*q*q^T
    Real dot = 0;
    for (long j = 0; j < M0.Dim(1); j++) dot += M0[i][j] * q[j];
    for (long j = 0; j < M0.Dim(1); j++) M0[i][j] -= dot * q[j];
  }
}

template <class Real>
std::vector<Real> new_nodes(int deg, int dim) {
  assert(deg < 20);
  assert(dim < 4);
  static std::vector<Real> y[20][4];

  if (!y[deg][dim].size()) {
#pragma omp critical(NEW_NODES)
    if (!y[deg][dim].size()) {
      unsigned int d = deg + 1;

      std::vector<Real> x(d);
      Real scal = 1.0 / pvfmm::cos<Real>(0.5 * pvfmm::const_pi<Real>() / d);
      for (int i = 0; i < d; i++)
        x[i] =
            -pvfmm::cos<Real>((i + (Real)0.5) * pvfmm::const_pi<Real>() / d) *
                scal * 0.5 +
            0.5;

      unsigned int n1 = pvfmm::pow<unsigned int>(d, dim);
      std::vector<Real> y_(n1 * dim);
      for (int i = 0; i < dim; i++) {
        unsigned int n2 = pvfmm::pow<unsigned int>(d, i);
        for (int j = 0; j < n1; j++) {
          y_[j * dim + i] = x[(j / n2) % d];
        }
      }

      {  // build sparse mesh
#if __USE_SPARSE_GRID__
        Real eps = 1.0;
        while (eps + (Real)1.0 > 1.0) eps *= 0.5;
        eps = sqrt(eps);

        assert(dim == 3);
        std::vector<Real> coord = y_;
        pvfmm::Matrix<Real> Mcoord(coord.size() / 3, 3);
        for (long i = 0; i < coord.size(); i++)
          Mcoord[0][i] = coord[i] * 2.0 - 1.0;

        long Ncoeff = (d * (d + 1) * (d + 2)) / 6;
        std::vector<Real> buff(d * (d + 3 * 2));
        pvfmm::Matrix<Real> M0(Mcoord.Dim(0), Ncoeff);
        for (long i = 0; i < Mcoord.Dim(0); i++) {
          Real *coord = &Mcoord[i][0];
          pvfmm::cheb_eval(d - 1, coord, &M0[i][0], &buff[0]);
        }

        pvfmm::Vector<int> flag(Mcoord.Dim(0));
        flag.SetZero();
        std::vector<Real> q(Ncoeff);
        std::vector<Real> norm(Mcoord.Dim(0));

        for (long i = 0; i < Ncoeff; i++) {
          long pivot = 0;
          {  // set pivot
#pragma omp parallel for schedule(static)
            for (long j = 0; j < Mcoord.Dim(0); j++) {
              Real dot = 0;
              for (long k = 0; k < Ncoeff; k++) dot += M0[j][k] * M0[j][k];
              norm[j] = sqrt(dot);
            }
            for (long j = 0; j < Mcoord.Dim(0); j++) {
              if (norm[j] > norm[pivot]) pivot = j;
            }
            if (norm[pivot] < eps) continue;
          }

          {  // Add all symmetric points
            long j0 = (pivot) % d;
            long j1 = (pivot / d) % d;
            long j2 = (pivot / d / d) % d;

            for (int k0 = 0; k0 < 2; k0++)
              for (int k1 = 0; k1 < 2; k1++)
                for (int k2 = 0; k2 < 2; k2++) {
                  long j0_ = (k0 ? d - j0 - 1 : j0);
                  long j1_ = (k1 ? d - j1 - 1 : j1);
                  long j2_ = (k2 ? d - j2 - 1 : j2);

                  long pivot_ = 0;
                  pivot_ = j0_ + d * (j1_ + d * j2_);
                  GramSchmidt(M0, flag, pivot_, q);
                  pivot_ = j1_ + d * (j2_ + d * j0_);
                  GramSchmidt(M0, flag, pivot_, q);
                  pivot_ = j2_ + d * (j0_ + d * j1_);
                  GramSchmidt(M0, flag, pivot_, q);
                  pivot_ = j2_ + d * (j1_ + d * j0_);
                  GramSchmidt(M0, flag, pivot_, q);
                  pivot_ = j1_ + d * (j0_ + d * j2_);
                  GramSchmidt(M0, flag, pivot_, q);
                  pivot_ = j0_ + d * (j2_ + d * j1_);
                  GramSchmidt(M0, flag, pivot_, q);
                }
          }
        }

        for (long i = 0; i < M0.Dim(0); i++) {  // Add boundary points
          long j0 = i % d;
          long j1 = (i / d) % d;
          long j2 = (i / d / d) % d;
          if (j0 == 0 || j1 == 0 || j2 == 0 || j0 == d - 1 || j1 == d - 1 ||
              j2 == d - 1) {
            GramSchmidt(M0, flag, i, q);
          }
        }

        y_.clear();
        for (long i = 0; i < flag.Dim(); i++) {
          if (flag[i]) {
            y_.push_back(coord[i * 3 + 0]);
            y_.push_back(coord[i * 3 + 1]);
            y_.push_back(coord[i * 3 + 2]);
          }
        }
#endif
      }

      y[deg][dim].swap(y_);
    }
  }

  return y[deg][dim];
}

template <class Real>
void GetPt2CoeffMatrix(int cheb_deg, pvfmm::Matrix<Real> &M) {
  long Ncoeff = (cheb_deg + 1) * (cheb_deg + 2) * (cheb_deg + 3) / 6;
  if (M.Dim(1) != Ncoeff) {
    std::vector<Real> coord = new_nodes<Real>(cheb_deg, 3);
    M.ReInit(coord.size() / 3, Ncoeff);

    int cheb_deg_ = (cheb_deg + 0);
    long Ncoeff_ = (cheb_deg_ + 1) * (cheb_deg_ + 2) * (cheb_deg_ + 3) / 6;
    pvfmm::Matrix<Real> M_(coord.size() / 3, Ncoeff_);
    std::vector<Real> buff((cheb_deg_ + 1) * (cheb_deg_ + 1 + 3 * 2));
    for (long i = 0; i < M_.Dim(0); i++) {
      for (int j = i * 3; j < (i + 1) * 3; j++) coord[j] = coord[j] * 2.0 - 1.0;
      pvfmm::cheb_eval(cheb_deg_, &coord[i * 3], &M_[i][0], &buff[0]);
    }
    M_ = M_.pinv().Transpose();

    long indx = 0, indx_ = 0;
    for (int i0 = 0; i0 <= cheb_deg_; i0++) {
      for (int i1 = 0; i0 + i1 <= cheb_deg_; i1++) {
        for (int i2 = 0; i0 + i1 + i2 <= cheb_deg_; i2++) {
          if (i0 + i1 + i2 <= cheb_deg) {
            for (long j = 0; j < M.Dim(0); j++) M[j][indx] = M_[j][indx_];
            indx++;
          }
          indx_++;
        }
      }
    }
  }
}

}  // namespace tbslas

#endif  //_TBSLAS_CHEB_H_
