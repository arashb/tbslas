#include <cheb_utils.hpp>

#ifndef _TBSLAS_CHEB_H_
#define _TBSLAS_CHEB_H_

#ifndef __USE_SPARSE_GRID__
#define __USE_SPARSE_GRID__ 1
#endif

namespace tbslas{

template <class Real>
std::vector<Real> new_nodes(int deg, int dim){

  assert(deg<20);
  assert(dim<4);
  static std::vector<Real> y[20][4];

  if(!y[deg][dim].size()){
    #pragma omp critical(NEW_NODES)
    if(!y[deg][dim].size()){
      unsigned int d=deg+1;
      std::vector<Real> x(d);
      Real scal=1.0/pvfmm::cos<Real>(0.5*pvfmm::const_pi<Real>()/d);
      for(int i=0;i<d;i++) x[i]=-pvfmm::cos<Real>((i+(Real)0.5)*pvfmm::const_pi<Real>()/d)*scal*0.5+0.5;

      unsigned int n1=pvfmm::pow<unsigned int>(d,dim);
      std::vector<Real> y_(n1*dim);
      for(int i=0;i<dim;i++){
        unsigned int n2=pvfmm::pow<unsigned int>(d,i);
        for(int j=0;j<n1;j++){
          y_[j*dim+i]=x[(j/n2)%d];
        }
      }

      { // build sparse mesh
        #if __USE_SPARSE_GRID__
        assert(dim==3);
        long Ncoeff=((deg+1)*(deg+2)*(deg+3))/6;
        std::vector<Real> buff((deg+1)*(deg+1+3*2));
        std::vector<Real> coord=y_;
        pvfmm::Matrix<Real> Mcoord(coord.size()/3,3);
        for(long i=0;i<coord.size();i++) Mcoord[0][i]=coord[i]*2.0-1.0;

        pvfmm::Matrix<Real> M0(Mcoord.Dim(0),Ncoeff);
        for(long i=0;i<Mcoord.Dim(0);i++){
          Real* coord=&Mcoord[i][0];
          pvfmm::cheb_eval(deg, coord, &M0[i][0], &buff[0]);
        }

        pvfmm::Vector<int> flag(Mcoord.Dim(0));
        flag.SetZero();
        for(long i=0;i<Ncoeff;i++){
          long max_indx=0;
          std::vector<Real> norm;
          for(long j=0;j<Mcoord.Dim(0);j++){
            pvfmm::Matrix<Real> q(1, Ncoeff, &M0[j][0]);
            pvfmm::Matrix<Real> norm2=q*q.Transpose();
            norm.push_back(sqrt(norm2[0][0]));
            if(norm.back()>norm[max_indx]) max_indx=j;
          }
          flag[max_indx]=1;
          pvfmm::Matrix<Real> q(1,Ncoeff,&M0[max_indx][0]);
          for(long j=0;j<Ncoeff;j++) q[0][j]*=1.0/norm[max_indx];
          M0=M0-(M0*q.Transpose())*q;
        }

        y_.clear();
        for(long i=0;i<flag.Dim();i++){
          if(flag[i]){
            y_.push_back(coord[i*3+0]);
            y_.push_back(coord[i*3+1]);
            y_.push_back(coord[i*3+2]);
          }
        }
        #endif
      }

      y[deg][dim].swap(y_);
    }
  }

  return y[deg][dim];
}

template <class T>
static pvfmm::Matrix<T>
NewPt2ChebPt_matrix(int cheb_deg){ // There are many other ways to construct this matrix
  int d=cheb_deg+1;
  pvfmm::Matrix<T> M(d,d);

#if 1 // Using Lagrange interpolation

  std::vector<T> x=pvfmm::cheb_nodes<T>(cheb_deg,1);
  std::vector<T> y=tbslas::new_nodes<T>(cheb_deg,1);
  for(int i=0;i<d;i++){
    for(int j=0;j<d;j++){
      T M_=1;
      for(int k=0;k<d;k++){
        if(i!=k) M_*=(x[j]-y[k])/(y[i]-y[k]);
      }
      M[i][j]=M_;
    }
  }

#else // More stable but only works for scaled Chebyshev nodes

  pvfmm::Matrix<T> M_p2c;
  { // Compute M_p2c
    std::vector<T> x(d);
    for(int i=0;i<d;i++) x[i]=-pvfmm::cos<T>((i+(T)0.5)*pvfmm::const_pi<T>()/d);

    std::vector<T> p(d*d);
    pvfmm::cheb_poly(cheb_deg,&x[0],d,&p[0]);
    for(int i=d;i<d*d;i++) p[i]=p[i]*2.0;
    for(int i=0;i<d*d;i++) p[i]=p[i]/d;
    pvfmm::Matrix<T> M(d,d,&p[0],false);
    M_p2c=M.Transpose();
  }

  pvfmm::Matrix<T> M_c2p;
  { // Compute M_c2p
    std::vector<T> x(d);
    T scal=pvfmm::cos<T>(0.5*pvfmm::const_pi<T>()/d); //// This is important
    for(int i=0;i<d;i++) x[i]=-pvfmm::cos<T>((i+(T)0.5)*pvfmm::const_pi<T>()/d)*scal;

    std::vector<T> p(d*d);
    pvfmm::cheb_poly(cheb_deg,&x[0],d,&p[0]);
    M_c2p.ReInit(d,d,&p[0]);
  }
  M=M_p2c*M_c2p;

#endif

  return M;
}

template <class T>
void NewPt2ChebPt(T* fn_v, int cheb_deg, int dof, pvfmm::mem::MemoryManager* mem_mgr=NULL) {
  int d=cheb_deg+1;

  // Precompute
  pvfmm::Matrix<T>* Mp=NULL;
  static std::vector<pvfmm::Matrix<T> > precomp;
  #pragma omp critical (CHEB_APPROX)
  {
    if(precomp.size()<=(size_t)d){
      precomp .resize(d+1);
    }
    if(precomp [d].Dim(0)==0 && precomp [d].Dim(1)==0){
      precomp[d]=NewPt2ChebPt_matrix<T>(cheb_deg);
    }
    Mp=&precomp[d];
  }

  // Create work buffers
  size_t buff_size=dof*d*d*d;
  T* buff=pvfmm::mem::aligned_new<T>(2*buff_size,mem_mgr);
  T* buff1=buff+buff_size*0;
  T* buff2=buff+buff_size*1;

  { // Apply Mp along x-dimension
    pvfmm::Matrix<T> Mi(dof*d*d,d,fn_v,false);
    pvfmm::Matrix<T> Mo(dof*d*d,d,buff2,false);
    pvfmm::Matrix<T>::GEMM(Mo, Mi, *Mp);

    pvfmm::Matrix<T> Mo_t(d,dof*d*d,buff1,false);
    for(size_t i=0;i<Mo.Dim(0);i++)
    for(size_t j=0;j<Mo.Dim(1);j++){
      Mo_t[j][i]=Mo[i][j];
    }
  }
  { // Apply Mp along y-dimension
    pvfmm::Matrix<T> Mi(d*dof*d,d,buff1,false);
    pvfmm::Matrix<T> Mo(d*dof*d,d,buff2,false);
    pvfmm::Matrix<T>::GEMM(Mo, Mi, *Mp);

    pvfmm::Matrix<T> Mo_t(d,d*dof*d,buff1,false);
    for(size_t i=0;i<Mo.Dim(0);i++)
    for(size_t j=0;j<Mo.Dim(1);j++){
      Mo_t[j][i]=Mo[i][j];
    }
  }
  { // Apply Mp along z-dimension
    pvfmm::Matrix<T> Mi(d*d*dof,d,buff1,false);
    pvfmm::Matrix<T> Mo(d*d*dof,d,buff2,false);
    pvfmm::Matrix<T>::GEMM(Mo, Mi, *Mp);

    pvfmm::Matrix<T> Mo_t(d,d*d*dof,buff1,false);
    for(size_t i=0;i<Mo.Dim(0);i++)
    for(size_t j=0;j<Mo.Dim(1);j++){
      Mo_t[j][i]=Mo[i][j];
    }
  }
  { // Set output
    pvfmm::Matrix<T> Mo(d*d*d,dof,buff1,false);
    pvfmm::Matrix<T> Mo_t(dof,d*d*d,fn_v,false);
    for(size_t i=0;i<Mo.Dim(0);i++)
    for(size_t j=0;j<Mo.Dim(1);j++){
      Mo_t[j][i]=Mo[i][j];
    }
  }

  // Free memory
  pvfmm::mem::aligned_delete<T>(buff,mem_mgr);
}

}//end namespace

#endif //_TBSLAS_CHEB_H_
