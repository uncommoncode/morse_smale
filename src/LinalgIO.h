//templated C-Style IO functions
#ifndef LINALGIO_H
#define LINALGIO_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#include "DenseMatrix.h"
#include "DenseVector.h"

namespace FortranLinalg{

template <typename TPrecision>
class LinalgIO{


  public:

  
  //Read Vector from header file
  static DenseVector<TPrecision> readVector(const std::string &filename){
    std::ifstream hdr;
    hdr.open(filename.c_str());
    
    std::string token;
    getline(hdr, token);
    if(token.compare("DenseVector") != 0){
      throw "Not a vector header file";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);
    int size = atoi(token.c_str());

    getline(hdr, token, ' ');
    getline(hdr, token);
    int elemSize = atoi(token.c_str());
    if(elemSize != sizeof(TPrecision)){
      throw "element size not equal to size of Vector precision";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);

    DenseVector<TPrecision> vector(size);
    size_t start = filename.find_last_of("/\\");
    std::stringstream ss;
    if(start != std::string::npos){
      ss << filename.substr(0, start+1);
    }
    ss << token;


    readVector(ss.str(), vector);

    return vector;
  };



  //Reads a  N Vector from binary file
  static bool readVector(const std::string &filename, DenseVector<TPrecision> &vector){
    std::ifstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.read((char*)vector.data(), sizeof(TPrecision)*vector.N());

    if(file.fail()){
      std::cout << "Reading failed" << std::endl;
      // get length of file:
      file.seekg (0, std::ios::beg);
      file.seekg (0, std::ios::end);
      long length = file.tellg();
      file.seekg (0, std::ios::beg);
      std::cout << length << std::endl;
      return false;  
    }  
    file.close();
    return true;
  };



  static void writeVector(const std::string &filename, DenseVector<TPrecision> &vector,
      bool writeHeader = true){
    std::ofstream file;
    file.open(filename.c_str(), std::ios::binary);
    char *data = (char*)vector.data();
    file.write(data, vector.N()*sizeof(TPrecision));
    file.close();

    if(writeHeader){
      size_t start = filename.find_last_of("/\\");
      std::string localFile;
      if(start != std::string::npos){
        localFile = filename.substr(start+1);
      }
      else{
        localFile = filename;
      };
      std::stringstream ss;
      ss << filename << ".hdr";

      std::ofstream hdr;
      hdr.open(ss.str().c_str());
      hdr << "DenseVector" << std::endl;
      hdr << "Size: " << vector.N() << std::endl;
      hdr << "ElementSize: " << sizeof(TPrecision) << std::endl;
      hdr << "DataFile: " << localFile << std::endl;
      hdr.close();
    }
  };



  //Read matrix from header file
  static DenseMatrix<TPrecision> readMatrix(const std::string &filename){
    std::ifstream hdr;
    hdr.open(filename.c_str());
    
    std::string token;
    getline(hdr, token);
    if(token.compare("DenseMatrix") != 0){
      throw "Not a matrix header file";
    }

    getline(hdr, token, ' ');
    getline(hdr, token, ' ');
    int m = atoi(token.c_str());
    getline(hdr, token, ' ');
    getline(hdr, token);
    int n = atoi(token.c_str());

    getline(hdr, token, ' ');
    getline(hdr, token);
    int elemSize = atoi(token.c_str());
    if(elemSize != sizeof(TPrecision)){
      throw "Element size not equal to size of matrix precision";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);
    bool rowMajor = atoi(token.c_str()) != 0;

    getline(hdr, token, ' ');
    getline(hdr, token);

    DenseMatrix<TPrecision> matrix(m, n);

    size_t start = filename.find_last_of("/\\");
    std::stringstream ss;
    if(start != std::string::npos ){
      ss << filename.substr(0, start+1);
    }
    ss << token;

    readMatrix(ss.str(), matrix, rowMajor);

    return matrix;
  }



  //Reads a M x N Matrix from binary file
  static bool readMatrix(const std::string &filename, DenseMatrix<TPrecision> &matrix,
      bool rowMajor = false){
    std::ifstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.read((char*)matrix.data(), sizeof(TPrecision)*matrix.N()*matrix.M());

    if(file.fail()){
      std::cout << "Reading binary matrix file failed" << std::endl;
      std::cout << matrix.M() << " x " << matrix.N() << std::endl;
      std::cout << sizeof(TPrecision) << std::endl;
      // get length of file
      file.seekg (0, std::ios::beg);
      file.seekg (0, std::ios::end);
      long length = file.tellg();
      file.seekg (0, std::ios::beg);
      std::cout << length << std::endl;
      return false;     
    }  
    
    file.close();

    //convert to column major if necessary
    if(rowMajor){
      DenseMatrix<TPrecision> a( matrix.N(), matrix.M());
      for(unsigned int i=0; i<a.M(); i++){
        for(unsigned int j=0; j<a.N(); j++){
          a(i, j) = matrix(j, i);
        }
      }
      matrix.deallocate();
      matrix = a;
    }

    return true;
  };




  static void writeMatrix(const std::string &filename, DenseMatrix<TPrecision> &matrix,
      bool writeHeader = true){
    std::ofstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.write((char *)matrix.data(), matrix.N()*matrix.M()*sizeof(TPrecision));
    file.close();
    
    if(writeHeader){
      size_t start = filename.find_last_of("/\\");
      std::string localFile;
      if(start ==std::string::npos){
        localFile = filename;
      }
      else{
        localFile = filename.substr(start+1);
      }
      
      std::stringstream ss;
      ss << filename << ".hdr";

      std::ofstream hdr;
      hdr.open(ss.str().c_str());

      hdr << "DenseMatrix" << std::endl;
      hdr << "Size: " << matrix.M() << " x " << matrix.N() << std::endl;
      hdr << "ElementSize: " << sizeof(TPrecision) << std::endl;
      hdr << "RowMajor: " << false << std::endl;
      hdr << "DataFile: " << localFile << std::endl;
      hdr.close();
    }
  };


};

}

#endif
