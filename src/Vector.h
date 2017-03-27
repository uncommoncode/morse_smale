#ifndef VECTOR_H
#define VECTOR_H

namespace FortranLinalg{

template <typename TPrecision>
class Vector{
  public:    
    
    virtual ~Vector(){
    };

    virtual TPrecision &operator()(unsigned int i) = 0;
    virtual unsigned int N() = 0;

    virtual void deallocate() = 0;
};

}

#endif
