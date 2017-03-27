#ifndef EUCLIDEANMETRIC_H
#define EUCLIDEANMETRIC_H

#include "Metric.h"
#include <math.h>


template<typename TPrecision>
class EuclideanMetric : public Metric<TPrecision>{
  public:
    virtual ~EuclideanMetric(){};

    TPrecision distance(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = distanceSquared(x1, x2);
      return sqrt(result); 
    };

    TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int i1, FortranLinalg::Matrix<TPrecision> &Y, int i2){
      TPrecision result = distanceSquared(X, i1, Y, i2);
      return sqrt(result); 
    };
    
    TPrecision distance(FortranLinalg::Matrix<TPrecision> &X, int i1, FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = distanceSquared(X, i1, x2);
      return sqrt(result); 
    };


    TPrecision distanceSquared(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = 0;
      TPrecision tmp =0;
      for(unsigned int i=0; i<x1.N(); i++){
        tmp = x1(i) - x2(i);
        result += tmp*tmp;
      }
      return result;
    };

    TPrecision distanceSquared(FortranLinalg::Matrix<TPrecision> &X, int i1, 
                               FortranLinalg::Matrix<TPrecision> &Y, int i2){
      TPrecision result = 0;
      TPrecision tmp =0;
      for(unsigned int i=0; i<X.M(); i++){
        tmp = X(i, i1) - Y(i, i2);
        result += tmp*tmp;
      }
      return result;
    };

    TPrecision distanceSquared(FortranLinalg::Matrix<TPrecision> &X, int i1,
        FortranLinalg::Vector<TPrecision> &x2){
      TPrecision result = 0;
      TPrecision tmp =0;
      for(unsigned int i=0; i< X.M(); i++){
        tmp = X(i, i1) - x2(i);
        result += tmp*tmp;
      }
      return result;
    };  

};
  

#endif
