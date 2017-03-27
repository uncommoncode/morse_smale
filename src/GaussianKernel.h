#ifndef GAUSSIANKERNEL_H
#define GAUSSIANKERNEL_H

#include "Kernel.h"
#include "EuclideanMetric.h"
#include "Linalg.h"

#include <cmath>

template <typename TPrecision>
class GaussianKernel : public Kernel<TPrecision, TPrecision>{

  private:
    EuclideanMetric<TPrecision> metric;
    TPrecision var;
    TPrecision ng;
    TPrecision nh;
    TPrecision c;
    unsigned int d; 
    FortranLinalg::DenseVector<TPrecision> diff;

  public:
  

    GaussianKernel(unsigned int dim=1):diff(dim){
      d=dim;
    };
   
  
    GaussianKernel(TPrecision sigma, int dim=1):diff(dim){
      d = dim;
      setKernelParam(sigma);
    };

    virtual ~GaussianKernel(){
      diff.deallocate();
    };


    TPrecision f(TPrecision dsquared){
      return c * exp( - dsquared / var);
    };
    

    TPrecision f(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      return f( metric.distanceSquared(x1, x2) );
    };



  
    TPrecision f(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Matrix<TPrecision> &X2, int i2){
      return f( metric.distanceSquared(X2, i2, x1 ) );
    };
  
  
    TPrecision f(FortranLinalg::Matrix<TPrecision> &X1, int i1, FortranLinalg::Matrix<TPrecision> &X2, int i2){
      return f( metric.distanceSquared(X1, i1, X2, i2 ) );
    };


    void grad(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2, FortranLinalg::Vector<TPrecision> &g){
      gradf(x1, x2, g); 
    };



    //Hessian
    FortranLinalg::DenseMatrix<TPrecision> hessian(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      using namespace FortranLinalg;
      DenseMatrix<TPrecision> H(d, d);
      hessian(x1, x2, H);
      return H;

    };    
    
    //Hessian
    void hessian(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision>
        &x2, FortranLinalg::DenseMatrix<TPrecision> H){
      using namespace FortranLinalg;

      Linalg<TPrecision>::Subtract(x1, x2, diff);

      TPrecision val = 0;
      for(unsigned int i=0; i<diff.N(); i++){
        val += diff(i)*diff(i);
      }       
      val = exp( -val / var)*nh;

      Linalg<TPrecision>::OuterProduct(diff, diff, H);
      for(unsigned int i=0; i<H.N(); i++){
        H(i, i) -= 1;
      }
      Linalg<TPrecision>::Scale(H, val, H);
      
    };
    
    //Hessian
    void hessian(FortranLinalg::Matrix<TPrecision> &X1, unsigned int i1,  
                 FortranLinalg::Matrix<TPrecision> &X2, unsigned int i2,
                 FortranLinalg::DenseMatrix<TPrecision> &H){
      using namespace FortranLinalg;

      Linalg<TPrecision>::Subtract(X1, i1, X2, i2, diff);

      TPrecision val = 0;
      for(unsigned int i=0; i<diff.N(); i++){
        val += diff(i)*diff(i);
      } 
      val = exp( -val / var)*nh;

      Linalg<TPrecision>::OuterProduct(diff, diff, H);
      for(unsigned int i=0; i<H.N(); i++){
        H(i, i) -= 1;
      }
      Linalg<TPrecision>::Scale(H, val, H);
    };

    //grad and function value
    TPrecision gradf(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2,
      FortranLinalg::Vector<TPrecision> &g){
      using namespace FortranLinalg;
  
      Linalg<TPrecision>::Subtract(x1, x2, g);
      TPrecision val = 0;
      for(unsigned int i=0; i<g.N(); i++){
        val += g(i)*g(i);
      } 
      val = exp( -val / var);
      for(unsigned int i=0; i<g.N(); i++){
        g(i) *= ng * val;
      }

      return val*c;
    };


    TPrecision gradf(FortranLinalg::Matrix<TPrecision> &x1, int i1, FortranLinalg::Matrix<TPrecision> &x2, int
        i2, FortranLinalg::Vector<TPrecision> &g){
      using namespace FortranLinalg;
  
      Linalg<TPrecision>::Subtract(x1, i1, x2, i2, g);
      TPrecision val = 0;
      for(unsigned int i=0; i<g.N(); i++){
        val += g(i)*g(i);
      } 
      val = exp( -val / var);
      for(unsigned int i=0; i<g.N(); i++){
        g(i) *= ng * val;
      }

      return val*c;
    };    
    
    TPrecision gradf(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Matrix<TPrecision> &x2, int
        i2, FortranLinalg::Vector<TPrecision> &g){
      using namespace FortranLinalg;
  
      Linalg<TPrecision>::Subtract(x1, x2, i2, g);
      TPrecision val = 0;
      for(unsigned int i=0; i<g.N(); i++){
        val += g(i)*g(i);
      } 
      val = exp( -val / var);
      for(unsigned int i=0; i<g.N(); i++){
        g(i) *= ng * val;
      }

      return val*c;
    };




    TPrecision gradKernelParam(FortranLinalg::Vector<TPrecision> &x1, FortranLinalg::Vector<TPrecision> &x2){
      return -2.0/(var*sqrt(var)) * f(x1, x2);         
    };
    
  
    void setKernelParam(TPrecision sigma){
      var = 2.0*sigma*sigma;
      c = 1.0;//1.0/ pow(2.0*M_PI*sigma*sigma, d/2.0) ;
      ng = -2.0/var*c;
      nh = -ng/(sigma*sigma);
    };

  
    
    TPrecision getKernelParam(){
      return sqrt(var/2.0);
    };

    GaussianKernel<TPrecision> &operator=(const GaussianKernel<TPrecision> &rhs){
      using namespace FortranLinalg;
      if(this == &rhs){
        return *this;
      }
      this->var = rhs.var;
      this->ng = rhs.ng;
      this->nh = rhs.nh;
      this->c = rhs.c;
      this->d = rhs.d;
      diff.deallocate();
      diff = DenseVector<TPrecision>(d);
      
      return *this;
    };


};
#endif
