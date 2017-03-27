#ifndef KERNELDENSITY_H
#define KERNELDENSITY_H


#include "DenseVector.h"
#include "DenseMatrix.h"
#include "Kernel.h"
#include "Linalg.h"

template<typename TPrecision>
class KernelDensity{
      
  public:
    KernelDensity(FortranLinalg::DenseMatrix<TPrecision> &data, Kernel<TPrecision, TPrecision> &k)
                    :X(data), kernel(k){
    };

    //returns unnormalized density
    double p(int j, int leaveout = -1 ){
      TPrecision wsum = 0;
      for(int i=0; i < X.N(); i++){
        if(leaveout == i) continue;
        wsum += kernel.f(X, j, X, i);
      }
      return wsum;
    };

    //retunrs unnormalized density
    double p(FortranLinalg::DenseMatrix<TPrecision> &T, int index, bool leaveout = false){
      using namespace FortranLinalg;
      TPrecision wsum = 0;
      for(unsigned int i=0; i < X.N(); i++){
        bool use = true;
        if(leaveout){
          use = ! Linalg<TPrecision>::IsColumnEqual(T, index, X, i);
        }
        if(use){
          wsum += kernel.f(T, index, X, i);
        }
      }
      return wsum;
    };


    FortranLinalg::DenseVector<TPrecision> p(FortranLinalg::DenseMatrix<TPrecision> e){
      using namespace FortranLinalg;
      DenseVector<TPrecision> res(e.N());
      DenseVector<TPrecision> tmp(e.M());
      for(int i=0; i< e.N(); i++){
        Linalg<TPrecision>::ExtractColumn(e, i, tmp);
        res(i) = p(tmp);
      }
      return res;
    };


    double p(FortranLinalg::DenseVector<TPrecision> &x, int leaveout = -1){
      TPrecision wsum = 0;
      for(unsigned int i=0; i < X.N(); i++){
        if(leaveout != i){
          wsum += kernel.f(x, X, i);
        }
      }
      return wsum;
    };

    void setData(FortranLinalg::DenseMatrix<TPrecision> &data){
      X = data;
    };

  private:
    FortranLinalg::DenseMatrix<TPrecision> X;
    Kernel<TPrecision, TPrecision> &kernel;
};



#endif
