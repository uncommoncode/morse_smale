#ifndef LAPACKDEFS_H
#define LAPACKDEFS_H



#ifndef FL_INT
#define FL_INT int
#endif

//blas & lapack routines
namespace lapack{

extern "C"{

  //---- matrix matrix multiply ----//
  void dgemm_(char *transa, char *transb, FL_INT *m, FL_INT *n, FL_INT *k,
      double *alpha, double *A, FL_INT *lda, double *B, FL_INT *ldb, double *beta,
      double *C, FL_INT *ldc);
    
  void sgemm_(char *transa, char *transb, FL_INT *m, FL_INT *n, FL_INT *k, 
      float *alpha, float *A, FL_INT *lda, float *B, FL_INT *ldb, float *beta, 
      float *C, FL_INT *ldc);
  



  //---- matrix vector multiply ----//
  void dgemv_(char *trans, FL_INT *m, FL_INT *n, double *alpha, double *A, FL_INT *lda, 
      double *x, FL_INT *incx, double *beta, double *y, FL_INT *incy);
    
  void sgemv_(char *trans, FL_INT *m, FL_INT *n, float *alpha, float *A, FL_INT *lda, 
      float *x, FL_INT *incx, float *beta, float *y, FL_INT *incy);





  //---- dot product ----/
  double ddot_(FL_INT *m, double *x, FL_INT *incx, double *y, FL_INT *incy);

  float sdot_(FL_INT *m, float *x, FL_INT *incx, float *y, FL_INT *incy);


  //---- QR decomposition ----//
  void dgeqrf_(FL_INT *m, FL_INT *n, double *a, FL_INT *lda, double *tau, double *work, FL_INT *lwork, FL_INT *info);
  void sgeqrf_(FL_INT *m, FL_INT *n, float *a, FL_INT *lda, float *tau, float *work, FL_INT *lwork, FL_INT *info);

  void dorgqr_(FL_INT *m, FL_INT *n, FL_INT *k, double *a, FL_INT *lda, double *tau, double *work, FL_INT *lwork, FL_INT *info);
  void sorgqr_(FL_INT *m, FL_INT *n, FL_INT *k, float *a, FL_INT *lda, float *tau, float *work, FL_INT *lwork, FL_INT *info);

  //---- LU decvomposition ----//
  void dgetrf_(FL_INT * m, FL_INT *n, double *a, FL_INT *lda, FL_INT *ipiv, FL_INT *info);

  void sgetrf_(FL_INT * m, FL_INT *n, float *a, FL_INT *lda, FL_INT *ipiv, FL_INT *info);
  
  //---- Cholesky decvomposition ----//
  void dpotrf_(char *u, FL_INT *n, double *a, FL_INT *lda, FL_INT *info);

  void spotrf_(char *u, FL_INT *n, float *a, FL_INT *lda, FL_INT *info);
  
  //---- inverse ---//
  void dgetri_(FL_INT *n, double *a, FL_INT *lda, FL_INT *ipiv, double *work, FL_INT
      *lwork, FL_INT *info);

  void sgetri_(FL_INT *n, float *a, FL_INT *lda, FL_INT *ipiv, float *work, FL_INT
      *lwork, FL_INT *info);


  //---- inverse symmetric positive definite matrix ----/
  void dpotri_(char *uplo, FL_INT *n, double *A, FL_INT *lda, FL_INT *info);
  
  void spotri_(char *uplo, FL_INT *n, float *A, FL_INT *lda, FL_INT *info);

  //---- Solve SPD linear equations ----//
  void dposv_(char *uplo, FL_INT *n, FL_INT *nrhs, double *A, FL_INT *lda, double *B,
      FL_INT *ldb, FL_INT *info );

  void sposv_(char *uplo, FL_INT *n, FL_INT *nrhs, float *A, FL_INT *lda, float *B,
      FL_INT *ldb, FL_INT *info );

  void dposvx_( char *fact, char *uplo, FL_INT *n, FL_INT *nrhs, double *A, FL_INT *lda,
      double *af, FL_INT *ldaf, char *equed, double *s, double *b, FL_INT *ldb, double
      *x, FL_INT *ldx, double *rcond, double *ferr, double *berr, double *work,
      FL_INT *iwork, FL_INT *info );
  
  void sposvx_( char *fact, char *uplo, FL_INT *n, FL_INT *nrhs, float *A, FL_INT *lda,
      float *af, FL_INT *ldaf, char *equed, float *s, float *b, FL_INT *ldb, float
      *x, FL_INT *ldx, float *rcond, float *ferr, float *berr, float *work,
      FL_INT *iwork, FL_INT *info );

  //---- Solve linear equations ----/
  void dgesv_( FL_INT *n, FL_INT *nrhs, double *A, FL_INT *lda, FL_INT *ipiv, double *b, FL_INT
      *ldb, FL_INT *info );

  void sgesv_( FL_INT *n, FL_INT *nrhs, float *A, FL_INT *lda, FL_INT *ipiv, float *b, FL_INT
      *ldb, FL_INT *info );

  //---- Least squares  or minimum norm with QR or LQ, A is assumed to have full rank----// 
  void dgels_(char *trans, FL_INT *m, FL_INT *n, FL_INT *nhrs, double *A, FL_INT *lda, double
      *B, FL_INT *ldb, double *work, FL_INT *lwork, FL_INT *info);

  void sgels_(char *trans, FL_INT *m, FL_INT *n, FL_INT *nhrs, float *A, FL_INT *lda, float 
      *B, FL_INT *ldb, float *work, FL_INT *lwork, FL_INT *info);

  //---- Least squares  or minimum norm with SVD// 
  void dgelsd_(FL_INT *m, FL_INT *n, FL_INT *nhrs, double *A, FL_INT *lda, double
      *B, FL_INT *ldb, double *S, double *rcond, FL_INT *rank, double *work, FL_INT
      *lwork, FL_INT *iwork, FL_INT *info);

  void sgelsd_(FL_INT *m, FL_INT *n, FL_INT *nhrs, float *A, FL_INT *lda, float
      *B, FL_INT *ldb, float *S, float *rcond, FL_INT *rank, float *work, FL_INT *lwork,
      FL_INT *iwork, FL_INT *info);


  //--- Cholesky condition number ---//
  void dpocon(char *uplo, FL_INT *n, double *a, FL_INT *lda, double *anorm, double
      *rcond, double *work, FL_INT *iwork, FL_INT *info);
  
  void spocon(char *uplo, FL_INT *n, float *a, FL_INT *lda, float *anorm, float
      *rcond, float *work, FL_INT *iwork, FL_INT *info);


   //---- SVD ---//
   void dgesvd_(char *JOBU, char *JOBVT, FL_INT* M, FL_INT *N, double *A, 
               FL_INT *LDA, double *S, double *U, FL_INT* LDU, double *VT, 
               FL_INT *LDVT, double *WORK, FL_INT *LWORK, FL_INT *INFO );
   
   void sgesvd_(char *JOBU, char *JOBVT, FL_INT* M, FL_INT *N, float *A, 
               FL_INT *LDA, float *S, float *U, FL_INT* LDU, float *VT, 
               FL_INT *LDVT, float *WORK, FL_INT *LWORK, FL_INT *INFO );
   
   
   void dgesdd_(char *JOBU, FL_INT* M, FL_INT *N, double *A, 
               FL_INT *LDA, double *S, double *U, FL_INT* LDU, double *VT, 
               FL_INT *LDVT, double *WORK, FL_INT *LWORK, FL_INT *IWORK, FL_INT *INFO );
   

   void sgesdd_(char *JOBU, FL_INT* M, FL_INT *N, float *A, 
               FL_INT *LDA, float *S, float *U, FL_INT* LDU, float *VT, 
               FL_INT *LDVT, float *WORK, FL_INT *LWORK, FL_INT *IWORK, FL_INT *INFO );





   //------ lapack symmetric eigensystem routine ----//
   void dsyevr_( char *jobz, char *range, char *uplo, FL_INT *N, double *A,
                    FL_INT *LDA, double *VL, double *VU, FL_INT *IL, FL_INT *IU, 
                    double *ABSTOL, FL_INT *M, double *W, double *Z, FL_INT *LDZ, 
                    FL_INT *ISUPPZ, double *WORK, FL_INT *LWORK, FL_INT *IWORK, 
                    FL_INT *LIWORK, FL_INT *INFO );

    void ssyevr_( char *jobz, char *range, char *uplo, FL_INT *N, float *A,
                    FL_INT *LDA, float *VL, float *VU, FL_INT *IL, FL_INT *IU, 
                    float *ABSTOL, FL_INT *M, float *W, float *Z, FL_INT *LDZ, 
                    FL_INT *ISUPPZ, float *WORK, FL_INT *LWORK, FL_INT *IWORK, 
                    FL_INT *LIWORK, FL_INT *INFO );
}; 

};

#endif
