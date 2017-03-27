#ifndef SQUAREDEUCLIDEANMETRIC_H
#define SQUAREDEUCLIDEANMETRIC_H

#include "Metric.h"
#include <math.h>


template<typename TPrecision>
class SquaredEuclideanMetric : public Metric<TPrecision>{
  public:
    virtual ~SquaredEuclideanMetric(){};

    TPrecision distance(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = 0;
      TPrecision tmp =0;
      for(unsigned int i=0; i<x1.N(); i++){
        tmp = x1(i) - x2(i);
        result += tmp*tmp;
      }
      return result;
    };

    TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int i1, FortranLinalg::Matrix<TPrecision> &Y, int i2){
      TPrecision result = 0;
      TPrecision tmp =0;
      for(unsigned int i=0; i<X.M(); i++){
        tmp = X(i, i1) - Y(i, i2);
        result += tmp*tmp;
      }
      return result;
    };

    TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int i1, FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = 0;
      TPrecision tmp = 0;
      for(unsigned int i=0; i<x2.N(); i++){
        tmp = X(i, i1) - x2(i);
        result += tmp*tmp;
      }
      return result;
    };  

};
  

#endif
