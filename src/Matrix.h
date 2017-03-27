#ifndef MATRIX_H
#define MATRIX_H

namespace FortranLinalg{

//Simple Matrix storage to abstract row and columnwise ordering
template <typename TPrecision>
class Matrix{

  public:
    virtual ~Matrix(){};
    virtual TPrecision &operator()(unsigned int i, unsigned int j) = 0;
    virtual unsigned int M() = 0;    
    virtual unsigned int N() = 0;    
    virtual void deallocate() = 0;
};

}

#endif
