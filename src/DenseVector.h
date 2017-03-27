#ifndef DENSEVECTOR_H
#define DENSEVECTOR_H

#include <cstddef>
#include "Vector.h"

namespace FortranLinalg{

//Simple Matrix storage to abstract row and columnwise ordering
template <typename TPrecision>
class DenseVector : public Vector<TPrecision>{


  public:
    DenseVector(){
      n = 0;
      a = NULL;
    };

    DenseVector(unsigned int nrows, TPrecision *data=NULL){
      n = nrows;
      a = data;
      if(a == NULL){
        a = new TPrecision[n];
      } 
    };

    virtual ~DenseVector(){
    };

    virtual TPrecision &operator()(unsigned int i){
     //if(i>=n) throw "Out of bounds";
     return a[i];
    };
    
    void setDataPointer(TPrecision *data){
      a = data;
    };


    unsigned int N(){
      return n;
    };

    TPrecision *data(){
      return a;
    };

    void deallocate(){
      if(a != NULL){
        delete[] a; 
        a = NULL; 
      }
    };

    void shorten(unsigned int newn){
      if(newn < n){
        n = newn;
      }
    };

  protected:
    //Access to data array
    TPrecision *a;
    unsigned int n;


};

}


#endif
