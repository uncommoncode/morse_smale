#ifndef MAXHEAP_H
#define MAXHEAP_H

#include "Heap.h"

#include <math.h>

template <typename T>
class MaxHeap : public Heap<T>{
  
  public:
    typedef Heap<T> Super;

    MaxHeap(int n):Heap<T>(n){};

    MaxHeap(T *k, int n):Heap<T>(k, n){
      Super::buildHeap();
    };
    
      
    void changeElement(int i, T &newElem){
      Super::elems[i] = newElem;
      while(i > 0 && Super::elems[Super::parent(i)] < Super::elems[i]){
        Super::exchange(i, Super::parent(i));
        i = Super::parent(i);  
      }    
    };
    
    void heapify(int i){
      int l = Super::left(i);
      int r = Super::right(i);

      int largest;
      if( l < Super::nElements && Super::elems[l] > Super::elems[i]){
        largest = l;
      }
      else{
        largest = i;
      }

      
      if( r < Super::nElements && Super::elems[r] > Super::elems[largest]){
        largest = r;
      }

      if(largest != i){
        Super::exchnage(i, largest);
        heapify(largest);
      }
    };

};

#endif
