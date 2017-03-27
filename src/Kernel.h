#ifndef KERNEL_H
#define KERNEL_H

#include "Vector.h"
#include "Matrix.h"


template <typename TPrecision, typename TKernelParam>
class Kernel{

  public:
    virtual ~Kernel(){};

    virtual TPrecision f(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2) = 0;
    virtual TPrecision f(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Matrix<TPrecision> &X2, int i2) = 0;
    virtual TPrecision f(FortranLinalg::Matrix<TPrecision> &X1, int i1, FortranLinalg::Matrix<TPrecision> &X2, int i2) = 0;


    virtual void grad(FortranLinalg::Vector<TPrecision> &x,
        FortranLinalg::Vector<TPrecision> &x2, FortranLinalg::Vector<TPrecision> &g) = 0;
    virtual TPrecision gradf(FortranLinalg::Vector<TPrecision> &x,
        FortranLinalg::Vector<TPrecision> &x2, FortranLinalg::Vector<TPrecision> &g) = 0;

    virtual TKernelParam gradKernelParam(FortranLinalg::Vector<TPrecision> &x, FortranLinalg::Vector<TPrecision> &x2) = 0;
    
    virtual void setKernelParam(TKernelParam param) = 0;

    virtual TKernelParam getKernelParam() = 0;

};

#endif
