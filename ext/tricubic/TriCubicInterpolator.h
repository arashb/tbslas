// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#ifndef LIKELY_TRI_CUBIC_INTERPOLATOR
#define LIKELY_TRI_CUBIC_INTERPOLATOR

#include <vector>
//#include "boost/smart_ptr.hpp"

namespace likely {
class TriCubicInterpolator {
    // Performs tri-cubic interpolation within a 3D periodic grid.
    // Based on http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.89.7835
 public:
    /* typedef boost::shared_array<double> DataCube; */
    typedef std::vector<double>& DataCube;
    // Initializes an interpolator using the specified datacube of length n1*n2*n3 where
    // data is ordered first along the n1 axis [0,0,0], [1,0,0], ..., [n1-1,0,0], [0,1,0], ...
    // If n2 and n3 are both omitted, then n1=n2=n3 is assumed. Data is assumed to be
    // equally spaced and periodic along each axis, with the coordinate origin (0,0,0) at
    // grid index [0,0,0].
    TriCubicInterpolator(DataCube data, double spacing, int n1, int n2 = 0, int n3 = 0);
    virtual ~TriCubicInterpolator();
    // Returns the interpolated data value for the specified x,y,z point. If the point lies
    // outside the box [0,n1*spacing) x [0,n2*spacing) x [0,n3*spacing), it will be folded
    // back assuming periodicity along each axis.
    double operator()(double x, double y, double z) const;
    // Returns the grid parameters.
    double getSpacing() const;
    int getN1() const;
    int getN2() const;
    int getN3() const;
 private:
    // Returns the unrolled 1D index corresponding to [i1,i2,i3] after mapping to each ik into [0,nk).
    // Assumes that i1 increases fastest in the 1D array.
    int _index(int i1, int i2, int i3) const;
    DataCube _data;
    double _spacing;
    int _n1, _n2, _n3;
    mutable int _i1, _i2, _i3;
    mutable double _coefs[64];
    mutable bool _initialized;
    static int _C[64][64];
  }; // TriCubicInterpolator

  inline double TriCubicInterpolator::getSpacing() const { return _spacing; }
  inline int TriCubicInterpolator::getN1() const { return _n1; }
  inline int TriCubicInterpolator::getN2() const { return _n2; }
  inline int TriCubicInterpolator::getN3() const { return _n3; }

  inline int TriCubicInterpolator::_index(int i1, int i2, int i3) const {
    if((i1 %= _n1) < 0) i1 += _n1;
    if((i2 %= _n2) < 0) i2 += _n2;
    if((i3 %= _n3) < 0) i3 += _n3;
    return i1 + _n1*(i2 + _n2*i3);
  }

}  // namespace likely

#endif // LIKELY_TRI_CUBIC_INTERPOLATOR
