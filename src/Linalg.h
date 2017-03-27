//Templated C-Style vector operations

#ifndef LINALG_H
#define LINALG_H

#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <limits>
#include <list>

#include "Matrix.h"
#include "DenseMatrix.h"
#include "Vector.h"
#include "DenseVector.h"

#include "LapackDefs.h"

#ifndef FL_INT
#define FL_INT int
#endif



namespace FortranLinalg{

template<typename TPrecision>
class Linalg{

//Checks bounds for matrix multuiplications, etc.
//make it slightly slower
#define LINALG_CHECK

#define isDoubleTPrecision() sizeof(TPrecision) == sizeof(double)

  public:

  
  //--------- Methods for DenseMatrices/Vectors using blas - fast -----------//

  
  //least squares with svd
  static DenseMatrix<TPrecision> LeastSquares(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, double *sse = NULL){


    FL_INT m = a.M();
    FL_INT n = a.N();
    FL_INT lda = m;

    FL_INT nrhs = b.N();
    FL_INT ldb = std::max(m, n);
    FL_INT info = 0;
    TPrecision workTmp = 0 ;
    FL_INT query = -1;
    bool deallocate = false; 
    DenseMatrix<TPrecision> tmp = b; 
    if(ldb > (FL_INT) b.M()){
      deallocate = true;
      tmp = DenseMatrix<TPrecision>(ldb, nrhs);
      for(unsigned FL_INT i=0; i<b.M(); i++){
        for(unsigned FL_INT j=0; j<b.N(); j++){
          tmp(i, j)= b(i, j);
        }
      }
    }

    TPrecision *s = new TPrecision[std::min(m, n)];
    TPrecision rcond = -1; //use machine precision as condition number
    FL_INT rank = -1;
    FL_INT *iwork = new FL_INT[100*std::min(m,n)];
    
    if(isDoubleTPrecision()){
      lapack::dgelsd_(&m, &n, &nrhs, (double*)a.data(), &lda,
          (double*)tmp.data(), &ldb,  (double*)s, (double*)&rcond, &rank,
          (double*)&workTmp, &query, iwork, &info);

      FL_INT lwork = workTmp;
      double *work =new double[lwork];
      lapack::dgelsd_(&m, &n, &nrhs, (double*)a.data(), &lda,
          (double*)tmp.data(), &ldb,  (double*)s, (double*)&rcond, &rank,
          (double*)work, &lwork, iwork, &info);
      delete[] work;
    }
    else{
      lapack::sgelsd_(&m, &n, &nrhs, (float*)a.data(), &lda,
          (float*)tmp.data(), &ldb,  (float*)s, (float*)&rcond, &rank,
          (float*)&workTmp, &query, iwork, &info);

      FL_INT lwork = workTmp;
      float *work =new float[lwork];
 
      lapack::sgelsd_(&m, &n, &nrhs, (float*)a.data(), &lda,
          (float*)tmp.data(), &ldb,  (float*)s, (float*)&rcond, &rank,
          (float*)work, &lwork, iwork, &info);
      delete[] work;
    }

    FL_INT nrows = a.N();
    DenseMatrix<TPrecision> out(nrows, nrhs);
    for(unsigned FL_INT i=0; i<out.M(); i++){
      for(unsigned FL_INT j=0; j<out.N(); j++){
        out(i, j) = tmp(i, j);
      }
    }

    if(sse != NULL){
      for(unsigned FL_INT i=0; i<out.N(); i++){
        sse[i] = 0;
        for(unsigned FL_INT j=out.M(); j < tmp.M(); j++){
          double blah = tmp(j, i);
          sse[i] += blah *  blah;
        }
      }
    }
    
    if(deallocate){
      tmp.deallocate();
    }


    delete[] s;
    delete[] iwork;

    return out;

  };



  //--- linear equations

  static DenseMatrix<TPrecision> Solve(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b){
        DenseMatrix<TPrecision> x(b.M(), b.N());
        Copy(b, x);
        Solve2(a, x);
        return x;
  };

  static bool Solve2(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b){

    FL_INT n = a.M();
    FL_INT nrhs = b.N();
    FL_INT lda = n;
    FL_INT ldb = b.M();
    FL_INT info = 0;
    FL_INT *ipiv = new FL_INT[n];
    

    if(isDoubleTPrecision()){
      lapack::dgesv_(&n, &nrhs, (double*)a.data(), &lda, ipiv, (double*)b.data(), &ldb, &info);
    }
    else{ 
      lapack::sgesv_(&n, &nrhs, (float*)a.data(), &lda, ipiv, (float*)b.data(), &ldb, &info);
    }

    delete[] ipiv;

    return info == 0;

  };

  //--- SPD linear equations

  static DenseMatrix<TPrecision> SolveSPD(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b){
        DenseMatrix<TPrecision> x(b.M(), b.N());
        Copy(b, x);
        SolveSPD2(a, x);
        return x;
  };
  
  static DenseMatrix<TPrecision> SolveSPD(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, TPrecision &rcond, FL_INT &info){
        DenseMatrix<TPrecision> x(b.M(), b.N());
        Copy(b, x);
        SolveSPD2(a, x, rcond, info);
        return x;
  };

  
  static bool SolveSPD2(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b){

    char u='U';
    FL_INT n = a.M();
    FL_INT nrhs = b.N();
    FL_INT lda = n;
    FL_INT ldb = b.M();
    FL_INT info = 0;
    

    if(isDoubleTPrecision()){
      lapack::dposv_(&u, &n, &nrhs, (double*)a.data(), &lda, (double*)b.data(), &ldb, &info);
    }
    else{ 
      lapack::sposv_(&u, &n, &nrhs, (float*)a.data(), &lda, (float*)b.data(), &ldb, &info);
    }

    return info == 0;

  };


  static void SolveSPD2(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b,
      TPrecision &rcond, FL_INT &info){

    char fact ='E';
    char u='U';
    FL_INT n = a.M();
    FL_INT nrhs = b.N();
    FL_INT lda = n;
    TPrecision *af = new TPrecision[n*n];
    FL_INT ldaf = n;
    char equed ='0';
    TPrecision *s = new TPrecision[n];
    FL_INT ldb = b.M();
    TPrecision *x = new TPrecision[n * nrhs];
    FL_INT ldx = n;
    TPrecision *ferr = new TPrecision[nrhs];
    TPrecision *berr = new TPrecision[nrhs];
    TPrecision *work = new TPrecision[3*n];
    FL_INT *iwork = new FL_INT[n];


    if(isDoubleTPrecision()){
      lapack::dposvx_(&fact, &u, &n, &nrhs, (double*)a.data(), &lda,
      (double*) af, &ldaf, &equed, (double*) s, (double*) b.data(), &ldb, 
      (double*) x, &ldx, (double*)&rcond, (double*)ferr, (double*)berr, (double*)work,
      iwork, &info );   
    }
    else{ 
      lapack::sposvx_(&fact, &u, &n, &nrhs, (float*)a.data(), &lda,
      (float*) af, &ldaf, &equed, (float*) s, (float*) b.data(), &ldb, 
      (float*) x, &ldx, (float*)&rcond, (float*)ferr, (float*)berr, (float*)work,
      iwork, &info );   
    }

    delete[] iwork;
    delete[] x;
    delete[] s;
    delete[] work;
    delete[] berr;
    delete[] ferr;
    delete[] af;

  };

  //--- Condition numbers  ---//


  //--- Determinants ---//
  
  static TPrecision DetSPD(DenseMatrix<TPrecision> &a){
    DenseMatrix<TPrecision> ch = Cholesky(a);
    TPrecision d = DetSPDCholesky(ch);
    ch.deallocate();
    return d;
  };
  
  static TPrecision DetSPDCholesky(DenseMatrix<TPrecision> &ch){
     TPrecision d = 1;
    for(unsigned FL_INT i=0; i< ch.N(); i++){
      d *= ch(i, i);
    }
    return d*d;
  };


  static TPrecision Det(DenseMatrix<TPrecision> &a){
    DenseMatrix<TPrecision> lu = LU(a);
    TPrecision d = DetLU(lu);
    lu.deallocate();
    return d;
  };
  
  static TPrecision DetLU(DenseMatrix<TPrecision> &lu){
    TPrecision d = 1;
    for(unsigned int i=0; i< lu.N(); i++){
      d *= lu(i, i);
    }
    return d;
  };

  //---- LU factorization ---//
  static DenseMatrix<TPrecision> LU(DenseMatrix<TPrecision> &a, FL_INT *ipiv = NULL){
    FL_INT m = a.M();
    FL_INT n = a.N();
    DenseMatrix<TPrecision> lu = Copy(a);
    FL_INT lda = m;

    bool clear = false;
    if(ipiv == NULL){
     clear = true;
     ipiv = new FL_INT[std::min(n, m)];
    }
    FL_INT info = 0;
    
    if(isDoubleTPrecision()){
      lapack::dgetrf_(&m, &n, (double*) lu.data(), &lda, ipiv, &info);
    }
    else{
      lapack::sgetrf_(&m, &n, (float*) lu.data(), &lda, ipiv, &info);
    }

    if(clear){
      delete[] ipiv;
    }

    if(info !=0 ){
      //std::cerr << "LU error: " << info << std::endl;
    }
    return lu;
  };



  //---- QRfactorization ---//
  static DenseMatrix<TPrecision> QR(DenseMatrix<TPrecision> &a){
    DenseMatrix<TPrecision> q = Copy(a);
    Linalg<TPrecision>::QR_inplace(q);
    return q;    
  };

  
  static void QR_inplace(DenseMatrix<TPrecision> &q){
    FL_INT n = q.N();
    FL_INT m = q.M();
    FL_INT lda = m;
    FL_INT info = 0;
    TPrecision *work = new TPrecision[1];
    FL_INT lwork = -1;
    TPrecision *tau = new TPrecision[n];

    //workspace query
    if(isDoubleTPrecision()){
       lapack::dgeqrf_(&m, &n, (double*) q.data(), &lda, (double*) tau, (double*) work, &lwork, &info);
    }
    else{ 
       lapack::sgeqrf_(&m, &n, (float*) q.data(), &lda, (float*) tau, (float*) work, &lwork, &info);
    }
    if(info !=0 ){
      //std::cerr << "QRF workspace error: " << info << std::endl;
    }

    lwork = work[0];
    delete[] work;
    work = new TPrecision[lwork];

    FL_INT k = std::min(m, n);
    //qr
    if( isDoubleTPrecision() ){
      lapack::dgeqrf_(&m, &n, (double*) q.data(), &lda, (double*) tau, (double*) work, &lwork, &info);
      if(info !=0 ){
        //std::cerr << "QRF error: " << info << std::endl;
      }

      n = std::min(m, n);
      lapack::dorgqr_(&m, &n, &k, (double*) q.data(), &lda, (double*) tau, (double*) work, &lwork, &info);
      if(info !=0 ){
        //std::cerr << "QR error: " << info << std::endl;
      }
    }
    else{
      lapack::sgeqrf_(&m, &n, (float*) q.data(), &lda, (float*) tau, (float*) work, &lwork, &info);
      if(info !=0 ){
        //std::cerr << "QRF error: " << info << std::endl;
      }

      n = std::min(m, n);
      lapack::sorgqr_(&m, &n, &k, (float*) q.data(), &lda, (float*) tau, (float*) work, &lwork, &info);
      if(info !=0 ){
        //std::cerr << "QR error: " << info << std::endl;
      }
    }


    delete[] work;
    delete[] tau;
  };



  //---- Choelsky factorization ---//
  static DenseMatrix<TPrecision> Cholesky(DenseMatrix<TPrecision> &a, char uplo = 'U'){
    FL_INT n = a.N();
    DenseMatrix<TPrecision> ch = Copy(a);
    FL_INT lda = n;
    FL_INT info = 0;
    
    if(isDoubleTPrecision()){
      lapack::dpotrf_(&uplo, &n, (double*) ch.data(), &lda, &info);
    }
    else{ 
      lapack::spotrf_(&uplo, &n, (float*) ch.data(), &lda, &info);
    }

    if(info !=0 ){
      //std::cerr << "Cholesky error: " << info << std::endl;
    }

    return ch;
  };



  //---- Inverse ---//
  static DenseMatrix<TPrecision> Inverse(DenseMatrix<TPrecision> &a){
    FL_INT *ipiv = new FL_INT[a.N()];
    DenseMatrix<TPrecision> inv = LU(a, ipiv);
    InverseLU(inv, ipiv);

    delete[] ipiv;

    return inv;
  };

  //---- Inverse, inpout is lu factorizatiopn ---//
  static void InverseLU(DenseMatrix<TPrecision> &inv, FL_INT *ipiv){
    FL_INT n =inv.N();
    FL_INT lda = n;
    TPrecision *work = new TPrecision[2];
    FL_INT info = 0;
    FL_INT lwork = -1;

    if(isDoubleTPrecision()){
      lapack::dgetri_(&n, (double*)inv.data(), &lda, ipiv,(double*) work, &lwork, &info);
    }
    else{
      lapack::sgetri_(&n, (float*)inv.data(), &lda, ipiv, (float*) work, &lwork, &info);
    }

    lwork = work[0];
    delete[] work;

    work = new TPrecision[lwork];

    if(isDoubleTPrecision()){
      lapack::dgetri_(&n, (double*)inv.data(), &lda, ipiv, (double*)work, &lwork, &info);
    }
    else{
      lapack::sgetri_(&n, (float*)inv.data(), &lda, ipiv, (float*)work, &lwork, &info);
    }


    delete[] work;
 
    if(info !=0 ){
      //std::cerr << "Inverse error: " << info << std::endl;
    }

  };

  //--- Symmetric positive definite inverse
  static DenseMatrix<TPrecision> InverseSPD(DenseMatrix<TPrecision> &a){
    char u ='U';
    DenseMatrix<TPrecision> inv = Cholesky(a, u);
    InverseCholesky(inv, u);
    return inv;
  }; 


  //--- Symmetric positive definite inverse
  static void InverseSPD(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &inv){
    char u ='U';
    Linalg<TPrecision>::Copy(a, inv);
    Cholesky(inv, u);
    InverseCholesky(inv, u);
  };   
  
  //--- Symmetric positive definite inverse, input is the choelsky decomposition
  static void InverseCholesky(DenseMatrix<TPrecision> &inv, char u = 'U'){

    FL_INT n =inv.N();
    FL_INT lda = n;
    FL_INT info = 0;
    if(isDoubleTPrecision()){
      lapack::dpotri_(&u, &n, (double*)inv.data(), &lda, &info);
    }
    else{
      lapack::spotri_(&u, &n, (float*)inv.data(), &lda, &info);
    }

    for(FL_INT i=0; i<n;i++){
      for(FL_INT j=i+1; j<n; j++){
        inv(j, i) = inv(i, j);
      }
    }
 
    if(info !=0 ){
      //std::cerr << "InverseCholesky error: " << info << std::endl;
    }

  };
   

  //--- Matrix Matrix multiplication methods 
  
  //Matrix matrix multiply without output allocation
  static void Multiply(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b,
      DenseMatrix<TPrecision> &c, bool transposeA = false, bool transposeB = false, 
      TPrecision alpha = 1 ){

    TPrecision beta = 0;

    char transa;
    char transb;
    FL_INT m;
    FL_INT n;
    FL_INT k;
    FL_INT lda;
    FL_INT ldb;

    if(transposeA){
      m = a.N();
      k = a.M();
      transa = 'T';

#ifdef LINALG_CHECK
      if(c.M() != a.N()) throw "Invalid matrix mult";
#endif

    }
    else{
      m = a.M();
      k = a.N();
      transa = 'N';

#ifdef LINALG_CHECK
      if(c.M() != a.M()) throw "Invalid matrix mult";
#endif

    }
    lda = a.M();

    if(transposeB){
      n = b.M();
      transb = 'T';

#ifdef LINALG_CHECK
      if(c.N() != b.M()) throw "Invalid matrix mult";
#endif

    }
    else{
      n = b.N();
      transb = 'N';

#ifdef LINALG_CHECK
      if(c.N() != b.N()) throw "Invalid matrix mult";
#endif
    }
    ldb = b.M();
    

    if(isDoubleTPrecision()){
      lapack::dgemm_(&transa, &transb, &m, &n, &k, (double*)&alpha, (double*)a.data(),
          &lda, (double*)b.data(), &ldb, (double*)&beta, (double*)c.data(), &m);
    }
    else{
      lapack::sgemm_(&transa, &transb, &m, &n, &k, (float*)&alpha, (float*)a.data(), &m,
                   (float*)b.data(), &ldb, (float*)&beta, (float*)c.data(), &m);
    }
  }; 
   




  //Matrix matrix multiply with output allocation
  static DenseMatrix<TPrecision> Multiply(DenseMatrix<TPrecision> &a, DenseMatrix<TPrecision> &b,
      bool transposeA = false, bool transposeB = false, TPrecision alpha = 1 ){

    FL_INT m;
    FL_INT n;

    if(transposeA){
      m = a.N();
    }
    else{
      m = a.M();
    }
    if(transposeB){
      n = b.M();
    }
    else{
      n = b.N();
    }

    DenseMatrix<TPrecision> c(m, n); 
    Linalg<TPrecision>::Multiply(a, b, c, transposeA, transposeB, alpha);
    return c;
  }; 
   


  
  //--- Matrix Vector multiplication methods  



  //Matrix Vector multiply with output allocation
  static DenseVector<TPrecision> Multiply(DenseMatrix<TPrecision> &a,
      DenseVector<TPrecision> &v, bool transpose = false, TPrecision alpha = 1){

    FL_INT lc;

    if(transpose){
      lc = a.N();
    }
    else{
      lc = a.M();
    }

    DenseVector<TPrecision> c(lc);
    Multiply(a, v, c, transpose, alpha);
    return c; 

  };

  //Matrix Vector multiply without output allocation
  static void Multiply(DenseMatrix<TPrecision> &a,
      DenseVector<TPrecision> &v, DenseVector<TPrecision> &out, 
      bool transpose = false, TPrecision alpha = 1){
  
    TPrecision beta = 0;

    char transa;
    FL_INT ma = a.M();
    FL_INT na = a.N();
    FL_INT lda = a.M();
    if(transpose){
      transa = 'T';
    }
    else{
      transa = 'N';
    }

    FL_INT inc = 1;
    if(isDoubleTPrecision()){
      lapack::dgemv_(&transa,  &ma, &na, (double*)&alpha, (double*)a.data(),
          &lda, (double*)v.data(), &inc, (double*)&beta, (double*)out.data(), &inc);
    }
    else{
      lapack::sgemv_(&transa,  &ma, &na, (float*)&alpha, (float*)a.data(),
          &lda, (float*)v.data(), &inc, (float*)&beta, (float*)out.data(), &inc);
    }
  }; 
 
 

  //Matrix vector (matrix column) multiply with output allocation 
  static DenseVector<TPrecision> MultiplyColumn(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, FL_INT index, bool transpose = false, TPrecision alpha = 1){

    FL_INT lc;

    if(transpose){
      lc = a.M();
    }
    else{
      lc = a.N();
    }

    DenseVector<TPrecision> c(lc);
    MultiplyColumn(a, b, index, c, transpose, alpha);
    return c; 
  }; 
  
  //Matrix vector (matrix column) multiply without output allocation 
  static void MultiplyColumn(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, FL_INT index, DenseVector<TPrecision> &out, 
      bool transpose = false, TPrecision alpha = 1){

    TPrecision beta = 0;

    char transa;
    if(transpose){
      transa = 'T';
    }
    else{
      transa = 'N';
    }
    FL_INT ma = a.M();
    FL_INT na = a.N();
    FL_INT lda = a.M();

    TPrecision *v = b.data();
    FL_INT vinc = 1;
    v = &v[index*b.M()];

    FL_INT inc = 1;

    if(isDoubleTPrecision()){
      lapack::dgemv_(&transa,  &ma, &na, (double*)&alpha, (double*)a.data(),
          &lda, (double*)v, &vinc, (double*)&beta, (double*)out.data(), &inc);
    }
    else{
      lapack::sgemv_(&transa,  &ma, &na, (float*)&alpha, (float*)a.data(),
          &lda, (float*)v, &vinc, (float*)&beta, (float*)out.data(), &inc);
    }
  }; 



  //Matrix vector (matrix row) multiply with output allocation 
  static DenseVector<TPrecision> MultiplyRow(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, FL_INT index, bool transpose = false, TPrecision alpha = 1){

    FL_INT lc;

      if(transpose){
        lc = a.M();
      }
      else{
        lc = a.N();
      }

    DenseVector<TPrecision> c(lc);
    MultiplyRow(a, b, index, c, transpose, alpha);
    return c; 
  }; 
  
  //Matrix vector (matrix row) multiply without output allocation 
  static void MultiplyRow(DenseMatrix<TPrecision> &a,
      DenseMatrix<TPrecision> &b, FL_INT index, DenseVector<TPrecision> &out, 
      bool transpose = false, TPrecision alpha = 1){

    TPrecision beta = 0;

    char transa;
    FL_INT ma = a.M();
    FL_INT na = a.N();
    if(transpose){
      transa = 'T';
    }
    else{
      transa = 'N';
    }
    FL_INT lda = a.M();

    
    TPrecision *v = b.data();
    FL_INT vinc = 1;
    if(b.isRowMajor()){
      v = &v[index*b.N];
    }
    else{
      v = &v[index];
      vinc = b.M();
    }


    FL_INT inc = 1;

    if(isDoubleTPrecision()){
      lapack::dgemv_(&transa,  &ma, &na, (double*)&alpha, (double*)a.data(),
          &lda, (double*)v, &vinc, (double*)&beta, (double*)out.data(), &inc);
    }
    else{
      lapack::sgemv_(&transa,  &ma, &na, (float*)&alpha, (float*)a.data(),
          &lda, (float*)v, &vinc, (float*)&beta, (float*)out.data(), &inc);
    }


  };






  //--- Dot prioduct methods


   static TPrecision Dot(DenseVector<TPrecision> &x, DenseVector<TPrecision> &y){
     
     FL_INT n = x.N();
     FL_INT inc = 1;
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double *)x.data(), &inc, (double*)y.data(), &inc);
     }
     else{
       res = (TPrecision) lapack::sdot_( &n, (float*)x.data(), &inc, (float*)y.data(), &inc);
     }
     return res;
   }; 

   

   static TPrecision DotColumn(DenseMatrix<TPrecision> &x, FL_INT index,
      DenseVector<TPrecision> &y){
    
      TPrecision *v = x.data();
      FL_INT vinc = 1;
      v = &v[index*x.M()];

     FL_INT n = y.N();
     FL_INT inc = 1;
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double *)v, &vinc, (double*)y.data(), &inc);
     }
     else{
       res = (TPrecision) lapack::sdot_(&n, (float*)v, &vinc, (float*)y.data(), &inc);
     }
     return res;
   }; 


   static TPrecision DotRow(DenseMatrix<TPrecision> &x, FL_INT index,
      DenseVector<TPrecision> &y){
    
      TPrecision *v = x.data();
      FL_INT vinc = 1;
      if(x.isRowMajor()){
        v = &v[index * x.N()];
      }
      else{
        v = &v[index];
        vinc = x.M();
      }

     FL_INT n = y.N();
     FL_INT inc = 1;
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double*)v, &vinc, (double*)y.data(), &inc);
     }
     else{
       res = (TPrecision) lapack::sdot_(&n, (float*)v, &vinc, (float*)y.data(), &inc);
     }
     return res;
   }; 


   static TPrecision DotRowRow(DenseMatrix<TPrecision> &x, FL_INT xindex,
      DenseMatrix<TPrecision> &y, FL_INT yindex ){
    
      TPrecision *v = x.data();
      FL_INT vinc = 1;
      if(x.isRowMajor()){
        v = &v[xindex * x.N()];
      }
      else{
        v = &v[xindex];
        vinc = x.M();
      }


      TPrecision *w = y.data();
      FL_INT winc = 1;
      if(y.isRowMajor()){
        w = &w[yindex * y.N()];
      }
      else{
        w = &w[yindex];
        winc = y.M();
      }


     FL_INT n = x.M();
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double*)v, &vinc, (double*)w, &winc);
     }
     else{
       res = (TPrecision) lapack::sdot_(&n, (float*)v, &vinc,(float*)w, &winc);
     }
     return res;
   }; 

   
   static TPrecision DotColumnRow(DenseMatrix<TPrecision> &x, FL_INT xindex,
      DenseMatrix<TPrecision> &y, FL_INT yindex ){
    return DotRowColumn(y, yindex, x, xindex);
   };


   static TPrecision DotRowColumn(DenseMatrix<TPrecision> &x, FL_INT xindex,
      DenseMatrix<TPrecision> &y, FL_INT yindex ){
    
      TPrecision *v = x.data();
      FL_INT vinc = 1;
      if(x.isRowMajor()){
        v = &v[xindex * x.N()];
      }
      else{
        v = &v[xindex];
        vinc = x.M();
      }


      TPrecision *w = y.data();
      FL_INT winc = 1;
      if(y.isRowMajor()){
        w = &v[yindex];
        winc = y.N();
      }
      else{
        w = &w[yindex*y.M()];
      }

     FL_INT n = x.M();
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double*)v, &vinc, (double*)w, &winc);
     }
     else{
       res = (TPrecision) lapack::sdot_(&n, (float*)v, &vinc,(float*)w, &winc);
     }
     return res;
   }; 

   
   static TPrecision DotColumnColumn(DenseMatrix<TPrecision> &x, FL_INT xindex,
      DenseMatrix<TPrecision> &y, FL_INT yindex ){
    
      TPrecision *v = x.data();
      FL_INT vinc = 1;

      v = &v[xindex*x.M()];

      TPrecision *w = y.data();
      FL_INT winc = 1;
     // if(y.isRowMajor()){
      //  w = &v[yindex];
      //  winc = y.N();
      //}
      //else{
        w = &w[yindex*y.M()];
      //}

     FL_INT n = x.M();
    
     TPrecision res = 0; 
     if( isDoubleTPrecision() ){
       res = (TPrecision) lapack::ddot_(&n, (double*)v, &vinc, (double*)w, &winc);
     }
     else{
       res = (TPrecision) lapack::sdot_(&n, (float*)v, &vinc,(float*)w, &winc);
     }
     return res;
   }; 









  //--------------- Methods for general matrices/vectors - slow ---------------// 




  //TODO
  /*static DenseMatrix<TPrecision> Multiply(Matrix<TPrecision> &a, Matrix<TPrecision> &b,
      bool transposeA = false, bool transposeB = false, TPrecision alpha = 0){

    DenseMatrix<TPrecision> c(a.M(), b.N(), false);
    
    return c; 
  };*/



  static DenseVector<TPrecision> ExtractColumn(Matrix<TPrecision> &a, FL_INT index){
    DenseVector<TPrecision> v(a.M());
    for(unsigned FL_INT i=0; i<a.M(); i++){
      v(i) = a(i, index);
    }
    return v;
  };
  
  static DenseMatrix<TPrecision> ExtractColumns(Matrix<TPrecision> &a, unsigned FL_INT start, unsigned FL_INT end){
    DenseMatrix<TPrecision> m(a.M(), end-start);
    for(unsigned FL_INT i=0; i<end-start; i++){
      for(unsigned FL_INT j=0; j< a.M(); j++){
        m(j,i) = a(j, start+i);
      }
    }
    return m;
  };


  static void ExtractColumn(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v){
    for(unsigned FL_INT i=0; i<std::min(a.M(),v.N()); i++){
      v(i) = a(i, index);
    }
  };
  

  static DenseVector<TPrecision> ExtractRow(Matrix<TPrecision> &a, FL_INT index){
    DenseVector<TPrecision> v(a.N());
    for(unsigned FL_INT i=0; i<a.N(); i++){
      v(i) = a(index, i);
    }
    return v;
  };


  static void ExtractRow(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v){
    for(unsigned FL_INT i=0; i<a.N(); i++){
      v(i) = a(index, i);
    }
  };
  



  static void SetColumn(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v){
    for(unsigned FL_INT i=0; i < std::min(a.M(), v.N()); i++){
      a(i, index) = v(i);
    }
  };
  
  
  static void SetColumns(Matrix<TPrecision> &T, FL_INT s1, FL_INT e1, DenseMatrix<TPrecision> &F, FL_INT s2 = 0){
    for(unsigned FL_INT i=0; i< e1-s1; i++){
      for(FL_INT j=0; j<T.M(); j++){
        T(j, s1+i) = F(j, s2+i);
      }
    }
  };


  static void SetRow(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v){
    for(unsigned FL_INT i=0; i<std::min(a.N(), v.N()); i++){
      a(index, i) = v(i);
    }
  };

  
  static void SetColumn(Matrix<TPrecision> &a, FL_INT index, TPrecision v){
    for(unsigned FL_INT i=0; i < a.M(); i++){
      a(i, index) = v;
    }
  };

  static void SetRow(Matrix<TPrecision> &a, FL_INT index, TPrecision v){
    for(unsigned FL_INT i=0; i<a.N(); i++){
      a(index, i) = v;
    }
  }; 



  static void SetColumn(Matrix<TPrecision> &a, FL_INT aindex,
      DenseMatrix<TPrecision> &b, FL_INT bindex){
    for(unsigned FL_INT i=0; i<std::min(a.M(), b.M()); i++){
      a(i, aindex) = b(i, bindex);
    }
  };



  static void SetRow(Matrix<TPrecision> &a, FL_INT aindex,
      DenseMatrix<TPrecision> &b, FL_INT bindex){
    for(unsigned FL_INT i=0; i<std::min(a.N(), b.N()); i++){
      a(aindex, i) = b(bindex, i);
    }
  };


  static void SetRowFromColumn(Matrix<TPrecision> &a, FL_INT aindex,
      DenseMatrix<TPrecision> &b, FL_INT bindex){
    for(unsigned FL_INT i=0; i<std::min(a.N(), b.M()); i++){
      a(aindex, i) = b(i, bindex);
    }
  };


  static DenseVector<TPrecision> SumColumns(Matrix<TPrecision> &a){
    DenseVector<TPrecision> sum(a.M());
    SumColumns(a, sum);
    return sum;
  };

  


  static void SumColumns(Matrix<TPrecision> &a, Vector<TPrecision> &v){
      for(unsigned FL_INT i=0; i < a.M(); i++){
        v(i) = 0;
        for(unsigned FL_INT j=0; j < a.N(); j++){
          v(i) += a(i, j);
        }
      }
  };

  
  static double SumColumn(Matrix<TPrecision> &a, FL_INT index){
    double sum = 0;
    for(unsigned FL_INT i=0; i < a.M(); i++){
      sum += a(i, index);
    }
    return sum;
  };
  
   static double Product(Vector<TPrecision> &a){
    double p = 1;
    for(unsigned FL_INT i=0; i < a.N(); i++){
      p *= a(i);
    }
    return p;
  };


  static double Sum(Vector<TPrecision> &a){
    double sum = 0;
    for(unsigned FL_INT i=0; i < a.N(); i++){
      sum += a(i);
    }
    return sum;
  };
   
  static double Sum(Matrix<TPrecision> &a){
    double sum = 0;
    for(unsigned FL_INT i=0; i < a.M(); i++){
      for(unsigned FL_INT j=0; j < a.N(); j++){
        sum += a(i,j);
      }
    }
    return sum;
  };
    


  static DenseVector<TPrecision> SumRows(Matrix<TPrecision> &a){
    DenseVector<TPrecision> sum(a.N());
    SumRows(a, sum);
    return sum;
  };



  static void SumRows(Matrix<TPrecision> &a, Vector<TPrecision> &v){
      for(unsigned FL_INT j=0; j < a.N(); j++){
        v(j) = 0;
        for(unsigned FL_INT i=0; i < a.M(); i++){
          v(j) += a(i, j);
        }
      }
  };




  static void SubtractColumnwise(Matrix<TPrecision> &a, Vector<TPrecision> &v,
      Matrix<TPrecision> &out){
      for(unsigned FL_INT i=0; i < a.M(); i++){
        for(unsigned FL_INT j=0; j < a.N(); j++){
          out(i, j) = a(i, j) - v(i);
        }
      }
  };
  


  static void AddColumnwise(Matrix<TPrecision> &a, Vector<TPrecision> &v,
      Matrix<TPrecision> &out){
      for(unsigned FL_INT i=0; i < a.M(); i++){
        for(unsigned FL_INT j=0; j < a.N(); j++){
          out(i, j) = a(i, j) + v(i);
        }
      }
  };
  


  static void SubtractColumn(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v,
      Matrix<TPrecision> &out)  {
      for(unsigned FL_INT i=0; i < a.M(); i++){
        out(i, index) = a(i, index) - v(i);
      }
  };


  static void AddColumn(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v,
      Matrix<TPrecision> &out){
      for(unsigned FL_INT i=0; i < a.M(); i++){
        out(i, index) = a(i, index) + v(i);
      }
  };

  static void SubtractRowwise(Matrix<TPrecision> &a, Vector<TPrecision> &v,
      Matrix<TPrecision> &out){
      for(unsigned FL_INT i=0; i < a.M(); i++){
        for(unsigned FL_INT j=0; j < a.N(); j++){
          out(i, j) = a(i, j) - v(j);
        }
      }
  };
  


  static void SubtractRow(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &v,
      Matrix<TPrecision> &out){
      for(unsigned FL_INT j=0; j < a.N(); j++){
        out(index, j) = a(index, j) - v(j);
      }
  };



  //result = a + b
  static void Add(Vector<TPrecision> &a, Vector<TPrecision> &b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) + b(i);
    }    
  };

  //result = a + B(:, index)
  static void Add(Vector<TPrecision> &a, Matrix<TPrecision> &B, FL_INT index, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) + B(i, index);
    }    
  };


  //B(:, index) = a + B(:, index)
  static void Add(Matrix<TPrecision> &B, FL_INT index, Vector<TPrecision> &a){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      B(i, index) = a(i) + B(i, index);
    }    
  };
  
  //C = A + B
  static void Add(Matrix<TPrecision> &A, Matrix<TPrecision> &B, Matrix<TPrecision> &C){
    for(unsigned FL_INT i = 0; i < A.N(); i++){
      for(unsigned FL_INT j = 0; j< A.M(); j++){
        C(j, i) = A(j, i) + B(j, i);
      }
    }    
  };    
  
  
  //C = A - B
  static void Subtract(Matrix<TPrecision> &A, Matrix<TPrecision> &B, Matrix<TPrecision> &C){
    for(unsigned FL_INT i = 0; i < A.N(); i++){
      for(unsigned FL_INT j = 0; j< A.M(); j++){
        C(j, i) = A(j, i) - B(j, i);
      }
    }    
  };     
  
  static void Scale(Matrix<TPrecision> &A, TPrecision s, Matrix<TPrecision> &C){
    for(unsigned FL_INT i = 0; i < A.N(); i++){
      for(unsigned FL_INT j = 0; j< A.M(); j++){
        C(j, i) = s*A(j, i);
      }
    }    
  };
  
  //result = a + s*b
  static void AddScale(Vector<TPrecision> &a, TPrecision s, Vector<TPrecision> &b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) + s * b(i);
    }    
  };
 


  static void AddScale(Vector<TPrecision> &a, TPrecision s, Matrix<TPrecision> &b,
      FL_INT index, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) + s * b(i, index);
    }    
  };

  static void AddScale(Matrix<TPrecision> &a, FL_INT i1, TPrecision s, Vector<TPrecision> &b,
      Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      result(i) = a(i,i1) + s * b(i);
    }    
  };

  static void AddScale(Matrix<TPrecision> &a, FL_INT i1, TPrecision s, Matrix<TPrecision> &b, FL_INT i2,
      Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      result(i) = a(i,i1) + s * b(i, i2);
    }    
  };


  
  static void AddScale(Matrix<TPrecision> &a, TPrecision s, Matrix<TPrecision> &b, Matrix<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      for(unsigned FL_INT j = 0; j < a.N(); j++){
        result(i, j) = a(i, j) + s * b(i, j);
      }
    }    
  };

  static void SubtractScale(Matrix<TPrecision> &a, TPrecision s, Matrix<TPrecision> &b, Matrix<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      for(unsigned FL_INT j = 0; j < a.N(); j++){
        result(i, j) = a(i, j) - s * b(i, j);
      }
    }    
  };



  
  static void ColumnAddScale(DenseMatrix<TPrecision> &m, FL_INT index, TPrecision s, Vector<TPrecision> &b){
    for(unsigned FL_INT i = 0; i < m.M(); i++){
      m(i, index) = m(i, index) + s * b(i);
    }    
  };

  static void ColumnAddScale(DenseMatrix<TPrecision> &m, FL_INT i1, TPrecision s, DenseMatrix<TPrecision> &b, FL_INT i2){
    for(unsigned FL_INT i = 0; i < m.M(); i++){
      m(i, i1) = m(i, i1) + s * b(i, i2);
    }    
  };


  static DenseVector<TPrecision> ColumnwiseSquaredNorm(Matrix<TPrecision> &A){
    DenseVector<TPrecision> norms(A.N());
    ColumnwiseSquaredNorm(A, norms);
    return norms;
  };
	  
  static void ColumnwiseSquaredNorm(Matrix<TPrecision> &A, Vector<TPrecision> &v){	  
    for(unsigned FL_INT i=0; i< A.N(); i++){
      TPrecision norm =0;
      TPrecision tmp =0;
      for(unsigned FL_INT j=0; j<A.M(); j++){
	      tmp = A(j, i);
        norm += tmp *tmp;
      }
      v(i) = norm;
    }       
  };


  static DenseVector<TPrecision> ColumnwiseNorm(Matrix<TPrecision> &A){
    DenseVector<TPrecision> norms(A.N());
    ColumnwiseNorm(A, norms);
    return norms;
  };
	  
  static void ColumnwiseNorm(Matrix<TPrecision> &A, Vector<TPrecision> &v){
    for(unsigned FL_INT i=0; i< A.N(); i++){
      TPrecision norm =0;
      TPrecision tmp =0;
      for(unsigned FL_INT j=0; j<A.M(); j++){
	tmp = A(j, i);
        norm += tmp *tmp;
      }
      v(i) = sqrt(norm);
    }       
  };



  //result = a - b
  static void Subtract(Vector<TPrecision> &a, TPrecision b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) - b;
    }    
  };

    //result = a - b
  static void Subtract(Matrix<TPrecision> &a, TPrecision b, Matrix<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      for(unsigned FL_INT j = 0; j < a.N(); j++){
        result(i,j) = a(i, j) - b;
      }
    }    
  };

  //result = a - b
  static void Add(Vector<TPrecision> &a, TPrecision b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      result(i) = a(i) + b;
    }    
  };




  //result = a - b
  static void Subtract(Vector<TPrecision> &a, Vector<TPrecision> &b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < std::min(a.N(), b.N()); i++){
      result(i) = a(i) - b(i);
    }    
  };



  //returns a - b
  static DenseVector<TPrecision> Subtract(Vector<TPrecision> &a, Vector<TPrecision> &b){
    DenseVector<TPrecision> result(a.N());
    Subtract(a, b, result);    
    return result;   
  };



  //result = a(:, index) - b
  static void Subtract(Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &b, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      result(i) = a(i, index) - b(i);
    }    
  };
  


  //result =b - a(:, index)
  static void Subtract(Vector<TPrecision> &b, Matrix<TPrecision> &a, FL_INT index, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      result(i) = b(i) - a(i, index);
    }    
  };



  //result = a(:, index) - b(:, index)
  static void Subtract(Matrix<TPrecision> &a, FL_INT aindex, Matrix<TPrecision> &b,
      FL_INT bindex, Vector<TPrecision> &result){
    for(unsigned FL_INT i = 0; i < a.M(); i++){
      result(i) = a(i, aindex) - b(i, bindex);
    }    
  };

  static DenseMatrix<TPrecision> Copy(DenseMatrix<TPrecision> &from){
    DenseMatrix<TPrecision> to(from.M(), from.N());
    Copy(from, to);
    return to;
  };



  static void Copy(Matrix<TPrecision> &from, Matrix<TPrecision> &to){
    for(unsigned FL_INT i = 0; i < std::min(from.N(), to.N()); i++){
      for(unsigned FL_INT j= 0; j < std::min(from.M(), to.M()); j++){
        to(j, i) = from(j, i);
      }
    } 
  };  
 

  static void CopyColumn(Matrix<TPrecision> &from, unsigned FL_INT fi, Matrix<TPrecision>
      &to, unsigned FL_INT ti){
    for(unsigned FL_INT j= 0; j < std::min(from.M(), to.M()); j++){
      to(j, ti) = from(j, fi);
    }
  };



  static DenseVector<TPrecision> Copy(Vector<TPrecision> &from){
    DenseVector<TPrecision> to(from.N());
    Copy(from, to);
    return to;
  };




  static void Copy(Vector<TPrecision> &from, Vector<TPrecision> &to){
    for(unsigned FL_INT i = 0; i < std::min( from.N(), to.N() ); i++){
      to(i) = from(i);
    } 
  };


  static TPrecision SquaredLengthColumn(Matrix<TPrecision> &v, FL_INT index){
    TPrecision l = 0;
    for(unsigned FL_INT i = 0; i < v.M(); i++){
      l += v(i, index) * v(i, index);
    }  
    return l;
  };


  static DenseVector<TPrecision> ColumnLengths(Matrix<TPrecision> &X){
    DenseVector<TPrecision> l(X.N());
    for(FL_INT i=0; i<l.N(); i++){
      l(i) = LengthColumn(X, i);
    }
    return l;
  };  
  
  static void ColumnLengths(Matrix<TPrecision> &X, DenseVector<TPrecision> &l){
    for(FL_INT i=0; i<l.N(); i++){
      l(i) = LengthColumn(X, i);
    }
  };

  static TPrecision LengthColumn(Matrix<TPrecision> &v, FL_INT index){
    TPrecision l = 0;
    for(unsigned FL_INT i = 0; i < v.M(); i++){
      l += v(i, index) * v(i, index);
    }  
    return sqrt(l);
  };

  static TPrecision LengthRow(Matrix<TPrecision> &v, FL_INT index){
    TPrecision l = 0;
    for(unsigned FL_INT i = 0; i < v.N(); i++){
      l += v(index, i) * v(index, i);
    }  
    return sqrt(l);
  };

  static TPrecision Length(Vector<TPrecision> &v){
    TPrecision l = 0;
    for(unsigned FL_INT i = 0; i < v.N(); i++){
      l += v(i) * v(i);
    }  
    return sqrt(l);
  };
  
  static TPrecision SquaredLength(Vector<TPrecision> &v){
    TPrecision l = 0;
    for(unsigned FL_INT i = 0; i < v.N(); i++){
      l += v(i) * v(i);
    }  
    return l;
  };



  static TPrecision Dot(Vector<TPrecision> &a, Vector<TPrecision> &b){
    TPrecision tmp = 0;
    for(unsigned FL_INT i = 0; i < a.N(); i++){
      tmp += b(i) * a(i);
    }  
    return tmp; 
  };

  static DenseMatrix<TPrecision> OuterProduct(Vector<TPrecision> &a, Vector<TPrecision> &b){
    DenseMatrix<TPrecision> op(a.N(), b.N());
    OuterProduct(a, b, op);
    return op;
  };


  static void OuterProduct(Vector<TPrecision> &a,
    Vector<TPrecision> &b, Matrix<TPrecision> &op){
    for(unsigned FL_INT i=0; i<a.N(); i++){
      for(unsigned FL_INT j=0; j<b.N(); j++){
        op(j,i) = a(i) * b(j);
      }
    }
  };  
  
  static void AddOuterProduct(Matrix<TPrecision> &c, Vector<TPrecision> &a,
    Vector<TPrecision> &b, Matrix<TPrecision> &op){
    for(unsigned FL_INT i=0; i<a.N(); i++){
      for(unsigned FL_INT j=0; j<b.N(); j++){
        op(j,i) = c(j, i) + a(i) * b(j);
      }
    }
  };  
  
  static void SubtractOuterProduct(Matrix<TPrecision> &c, Vector<TPrecision> &a,
    Vector<TPrecision> &b, Matrix<TPrecision> &op){
    for(unsigned FL_INT i=0; i<a.N(); i++){
      for(unsigned FL_INT j=0; j<b.N(); j++){
        op(j,i) = c(j, i) - a(i) * b(j);
      }
    }
  };  
  
  static void AddOuterProduct(Matrix<TPrecision> &c, Matrix<TPrecision> &a,
      unsigned FL_INT ai, Matrix<TPrecision> &b, unsigned FL_INT bi, 
      Matrix<TPrecision> &op){

    for(unsigned FL_INT i=0; i<a.M(); i++){
      for(unsigned FL_INT j=0; j<b.M(); j++){
        op(j,i) = c(j,i) + a(i, ai) * b(j, bi);
      }
    }
  };  
  
  static void AddOuterProduct(Matrix<TPrecision> &c, Vector<TPrecision> &a, 
      Matrix<TPrecision> &b, unsigned FL_INT bi, Matrix<TPrecision> &op){

    for(unsigned FL_INT i=0; i<a.N(); i++){
      for(unsigned FL_INT j=0; j<b.M(); j++){
        op(j,i) = c(j,i) + a(i) * b(j, bi);
      }
    }
  };






  //result = v*s
  static void Scale(Vector<TPrecision> &v, TPrecision s, Vector<TPrecision> &result){
    for(unsigned FL_INT i=0; i < v.N(); i++){
      result(i) = v(i)*s;
    }
  };   
  

  static void Sqrt(Vector<TPrecision> &v, Vector<TPrecision> &result){
    for(unsigned FL_INT i=0; i < v.N(); i++){
      result(i) = sqrt(v(i));
    }
  };   
  
  //returns 0 vector for vectors with 1/length = inf
  static void Normalize(Vector<TPrecision> &v){
      TPrecision l = 1.0/Length(v);
      if( l!=l || 
          std::numeric_limits<TPrecision>::infinity() == l ){
        l = 0;
      }
      Scale(v, l, v);
  };
  
  static void NormalizeColumns(Matrix<TPrecision> &m){
    for(FL_INT i=0; i< m.N(); i++){
      NormalizeColumn(m, i);
    }
  };

  static void NormalizeColumn(Matrix<TPrecision> &V, FL_INT index){
      TPrecision l = 1.0/LengthColumn(V, index);
      if( l!=l || 
          std::numeric_limits<TPrecision>::infinity() == l ){
        l = 0;
      }
      ScaleColumn(V, index, l);
  }; 
  



  static void ScaleColumn(Matrix<TPrecision> &m, FL_INT index, TPrecision s){
    for(unsigned FL_INT i=0; i < m.M(); i++){
      m(i,index) = m(i, index) * s;
    }
  };

  static void ScaleRow(DenseMatrix<TPrecision> &m, FL_INT index, TPrecision s){
    ScaleRow(m, index, s, m);
  };
  
  static void ScaleRow(DenseMatrix<TPrecision> &m, FL_INT index, TPrecision s, DenseMatrix<TPrecision> &out){
    for(unsigned FL_INT i=0; i < m.N(); i++){
      out(index, i) = m(index, i) * s;
    }
  };



  static void Set(Vector<TPrecision> &v, TPrecision s){
    for(unsigned FL_INT i=0; i < v.N(); i++){
      v(i) = s;
    }
  };
  


  static void Set(Matrix<TPrecision> &a, TPrecision s){
    for(unsigned FL_INT i=0; i<a.M(); i++){
      for(unsigned FL_INT j=0; j<a.N(); j++){
        a(i, j) = s;
      }
    }
  };



  static void Zero(DenseMatrix<TPrecision> &m){
    memset(m.data(), 0, m.N() * m.M() * sizeof(TPrecision));
  };
  


  static void Zero(DenseVector<TPrecision> &v){
    memset(v.data(), 0, v.N() * sizeof(TPrecision));
  };
  
  static void Zero(Vector<TPrecision> &v){
    Set(v, 0);
  };  
  
  static void Zero(Matrix<TPrecision> &m){
    Set(m, 0);    
  };
  





  static void Scale(DenseMatrix<TPrecision> &m, TPrecision s,
      DenseMatrix<TPrecision> &r){
    TPrecision *a = m.data();
    TPrecision *b = r.data();
    unsigned FL_INT l = m.N() * m.M();
    for(unsigned FL_INT i=0; i<l; i++){
      b[i] = a[i]*s;
    }
  };



  static TPrecision FrobeniusNorm(Matrix<TPrecision> &m){
    return sqrt(SquaredFrobeniusNorm(m));
  };

  static TPrecision SquaredFrobeniusNorm(Matrix<TPrecision> &m){
    TPrecision tmp = 0;
    TPrecision norm = 0;
    for(unsigned FL_INT i=0; i<m.M(); i++){
      for(unsigned FL_INT j=0; j<m.N(); j++){
        tmp = m(i, j);
        norm += tmp*tmp;
      }
    }
    return norm;
  };

  
  static DenseMatrix<TPrecision> Identity(FL_INT n){
    DenseMatrix<TPrecision> eye(n, n);
    eye.zero();
    for(unsigned FL_INT i=0; i<n; i++){
      eye(i, i) = 1;
    }
    return eye;
  };

  static DenseMatrix<TPrecision> Transpose(DenseMatrix<TPrecision> A){
    DenseMatrix<TPrecision> B(A.N(), A.M());
    Transpose(A, B); 
    return B;
  };

  static void Transpose(DenseMatrix<TPrecision> A, DenseMatrix<TPrecision> B){
    for(unsigned FL_INT i=0; i<A.M(); i++){
      for(unsigned FL_INT j=0; j<A.N(); j++){
        B(j, i) = A(i, j);
      }
    }
  };



  static bool IsColumnEqual(DenseMatrix<TPrecision> A, FL_INT ai, DenseMatrix<TPrecision>
      B, FL_INT bi){
    bool equal = true;
    for(unsigned FL_INT i=0; i< A.M() && equal; i++){
      equal = A(i, ai) == B(i, bi);
    }
    return equal;
  };



  static DenseVector<TPrecision> RowMax(Matrix<TPrecision> &A){
    DenseVector<TPrecision> m = Linalg<TPrecision>::ExtractColumn(A, 0);
    for(unsigned  FL_INT i=0;i<A.N(); i++){
      for(unsigned FL_INT j=0; j<A.M(); j++){
        if(A(j, i) > m(j)){
          m(j) = A(j, i);
        }
      }
    }
    return m;
  };




  static DenseVector<TPrecision> RowMin(Matrix<TPrecision> &A){
    DenseVector<TPrecision> m = Linalg<TPrecision>::ExtractColumn(A, 0);
    for(unsigned FL_INT i=0;i<A.N(); i++){
      for(unsigned FL_INT j=0; j<A.M(); j++){
        if(A(j, i) < m(j)){
          m(j) = A(j, i);
        }
      }
    }
    return m;
  };


  static TPrecision MaxAll(Matrix<TPrecision> &A){
    TPrecision m = -std::numeric_limits<TPrecision>::max();
    for(unsigned  FL_INT i=0;i<A.N(); i++){
      for(unsigned FL_INT j=0; j<A.M(); j++){
        if(A(j, i) > m){
          m = A(j, i);
        }
      }
    }
    return m;
  };

  static TPrecision MinAll(Matrix<TPrecision> &A){
    TPrecision m = std::numeric_limits<TPrecision>::max();
    for(unsigned  FL_INT i=0;i<A.N(); i++){
      for(unsigned FL_INT j=0; j<A.M(); j++){
        if(A(j, i) < m){
          m = A(j, i);
        }
      }
    }
    return m;
  };







  static TPrecision Max(Vector<TPrecision> &A){
    TPrecision m = A(0);
    for(unsigned FL_INT i=0;i<A.N(); i++){
        if(A(i) > m){
          m = A(i);
        }
    }
    return m;
  };
  
  static TPrecision Min(Vector<TPrecision> &A){
    TPrecision m = A(0);
    for(unsigned FL_INT i=0;i<A.N(); i++){
        if(A(i) < m){
          m = A(i);
        }
    }
    return m;
  };


  static TPrecision MaxColumn(Matrix<TPrecision> &A, FL_INT index){
    TPrecision m = A(0, index);
    for(unsigned FL_INT i=1;i<A.M(); i++){
        if(A(i, index) > m){
          m = A(i, index);
        }
    }
    return m;
  };
  
  static TPrecision MinColumn(Matrix<TPrecision> &A, FL_INT index){
    TPrecision m = A(0, index);
    for(unsigned FL_INT i=1;i<A.M(); i++){
        if(A(i, index) < m){
          m = A(i, index);
        }
    }
    return m;
  };




  static void Print(DenseMatrix<TPrecision> &m){
    for(FL_INT i=0; i<m.M(); i++){
      for(FL_INT j=0; j< m.N(); j++){
        std::cout << m(i, j) << ", ";
      }
      std::cout << std::endl;
    }
  };

    static void Print(DenseVector<TPrecision> &v){
      for(FL_INT j=0; j< v.N(); j++){
        std::cout << v(j) << ", ";
      }
      std::cout << std::endl;
  };



  static DenseMatrix<TPrecision> ToMatrix(std::list< DenseVector<TPrecision> > v){
    DenseMatrix<TPrecision> res;
    if( v.empty() ){
      return res;
    }

    typename std::list< DenseVector<TPrecision> >::iterator it = v.begin();
    res = DenseMatrix<TPrecision>( (*it).N(), v.size() );
    for(FL_INT i=0; it != v.end(); ++it, ++i){
      Linalg<TPrecision>::SetColumn(res, i, *it);
    }
    return res;
  };

 static DenseVector<TPrecision> ToVector(std::list< TPrecision > v){
    DenseVector<TPrecision> res;
    if( v.empty() ){
      return res;
    }

    res = DenseVector<TPrecision>( v.size() );
    typename std::list< TPrecision>::iterator it = v.begin();
    for(FL_INT i=0; it != v.end(); ++it, ++i){
      res(i) = *it;
    }
    return res;
 };


 static DenseMatrix<TPrecision> Center(DenseMatrix<TPrecision> X){
 
        DenseMatrix<TPrecision> Xc = Linalg<TPrecision>::Copy(X);
        DenseVector<TPrecision> c = Linalg<TPrecision>::SumColumns(Xc);
        Linalg<TPrecision>::Scale(c, 1.0/Xc.N(), c);
        Linalg<TPrecision>::SubtractColumnwise(Xc, c, Xc);
        c.deallocate();
        return Xc;
 };


 static DenseVector<TPrecision> Expand(DenseVector<TPrecision> &v, FL_INT
     newLength, TPrecision fill = 0){
   DenseVector<TPrecision> a(newLength);
   Linalg<TPrecision>::Copy(v, a);
   for(FL_INT i=v.N(); i<newLength; i++){
     a(i) = fill;
   }
   return a;
 };

 static DenseMatrix<TPrecision> ExpandColumns(DenseMatrix<TPrecision> A, FL_INT
     ncol, TPrecision fill = 0){
   DenseMatrix<TPrecision> M(A.M(), ncol);
   Linalg<TPrecision>::Copy(A, M);
   for(FL_INT i=A.N(); i<M.N(); i++ ){
     for(FL_INT j=0; j<M.M(); j++){
       M(j, i) = fill;
     }
   }
   return M;
 };  
 
 static DenseMatrix<TPrecision> ExpandRows(DenseMatrix<TPrecision> A, FL_INT
     nRow, TPrecision fill = 0){
   DenseMatrix<TPrecision> M(nRow, A.N());
   Linalg<TPrecision>::Copy(A, M);
   for(FL_INT i=0; i<M.N(); i++ ){
     for(FL_INT j=A.M(); j<M.M(); j++){
       M(j, i) = fill;
     }
   }
   return M;
 };

 static DenseMatrix<TPrecision> Expand(DenseMatrix<TPrecision> A, FL_INT
     nCol, FL_INT nRow, TPrecision fill = 0){
   DenseMatrix<TPrecision> M(nRow, nCol);
   Linalg<TPrecision>::Copy(A, M);
   for(FL_INT i=A.N(); i<M.N(); i++ ){
     for(FL_INT j=A.M(); j<M.M(); j++){
       M(j, i) = fill;
     }
   }
   return M;
 };  

 
};

};

#endif
