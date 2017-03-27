#ifndef METRIC_H
#define METRIC_H

#include "Matrix.h"
#include "Vector.h"

template<typename TPrecision>
class Metric{
    
  public:
    virtual ~Metric(){};
    virtual TPrecision distance(FortranLinalg::Vector<TPrecision> &x1,
        FortranLinalg::Vector<TPrecision> &x2) = 0;
    virtual TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int ix,
        FortranLinalg::Matrix<TPrecision> &Y, int iy) = 0;
    virtual TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int i1,
        FortranLinalg::Vector<TPrecision> &x2) = 0;
};


#endif
