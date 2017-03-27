#ifndef MINHEAP_H
#define MINHEAP_H

#include "Heap.h"

#include <math.h>

template <typename T>
class MinHeap : public Heap<T>{
  
  public:
    typedef Heap<T> Super;

    MinHeap(int n):Heap<T>(n){};

    MinHeap(T *k, int n):Heap<T>(n){
      Super::buildHeap(k, n);
    };

    virtual ~MinHeap(){};
   
  protected: 
    void changeElement(int i){
      while(i > 0 && Super::bigger( Super::parent(i), i) ){
        Super::exchange(i, Super::parent(i));
        i = Super::parent(i);  
      }    
    };
    
    void heapify(int i){
      int l = Super::left(i);
      int r = Super::right(i);

      int smallest = i;
      if( l < Super::getNumberOfElements() && Super::smaller(l, smallest) ){
        smallest = l;
      }

      
      if( r < Super::getNumberOfElements() && Super::smaller( r,smallest) ){
        smallest = r;
      }
     
      if(smallest != i){
        Super::exchange(i, smallest);
        heapify(smallest);
      }
    };

};

#endif
