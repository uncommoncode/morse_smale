#ifndef HEAP_H
#define HEAP_H

#include <math.h>

template <typename T>
class Heap{
  
  public:
    Heap(int n){
      elems = new T[n];
      orig2heap = new int[n];
      heap2orig = new int[n];
      for(int i=0; i < n; i++){
        orig2heap[i] = i;
        heap2orig[i] = i;
      }
      capacity = n;
      nElements = 0;
    };

    virtual ~Heap(){
      delete[] elems;
      delete[] heap2orig;
      delete[] orig2heap;
    }

    void init(int nElems, double val){
      nElements = nElems;
       for(int i=0; i < nElements; i++){
        elems[i] = val;
      }
    }
    
    void init(int nElems, double val, int root, double rootval){
      nElements = nElems;
       for(int i=0; i < nElements; i++){
        elems[i] = val;
      }
      elems[root] = rootval;
      exchange(0, root); 
    }

    void buildHeap(T *a, int n){
      for(int i = 0; i< n; i++){
        elems[i] = a[i];
      }
      nElements = n;
      for(int i= nElements/2; i>=0; i--){
        heapify(i);
      }
    };

    int parent(int i){
      return (i-1)/2;
    };

    int left(int i){
      return 2*i + 1;
    };

    int right(int i){
      return 2*i + 2;
    };
    
     
    bool insert(const T &elem){
      if(nElements == capacity){
        return false;
      }


      elems[nElements] = elem;
      orig2heap[nElements] = nElements;
      heap2orig[nElements] = nElements;
      nElements++;
      changeElement(nElements-1);
       

      return true;
    };
   
    int getRootIndex(){
      return heap2orig[0];
    }; 
    
    int extractRootIndex(){
      int root= heap2orig[0];
      exchange(0, nElements-1);
      nElements--;
      heapify(0);
      return root;

    }; 

    T getRoot(){
      return elems[heap2orig[0]];
    };

    T extractRoot(){
      T root = elems[heap2orig[0]];
      exchange(0, nElements-1);
      nElements--;
      heapify(0);
      return root;
    };

    void changeOrigElement(int i, const T &elem){
      elems[i] = elem;
      changeElement(orig2heap[i]);
    };

    void changeHeapElement(int i, const T &elem){
      elems[heap2orig[i]] = elem;
      changeElement(i);
    };

     

    int getCapacity(){
      return capacity;
    };

    int getNumberOfElements(){
      return nElements;
    };

    
    T &operator[](int i){
      return elems[heap2orig[i]];
    };


  protected:

    bool smaller(int i, int j){
      return elems[heap2orig[i]] < elems[heap2orig[j]];
    };
    
    bool bigger(int i, int j){
      return elems[heap2orig[i]] > elems[heap2orig[j]];
    };


    void exchange(int i, int j){
      int tmp = heap2orig[i];
      heap2orig[i] = heap2orig[j];
      heap2orig[j] = tmp;
      
      orig2heap[heap2orig[i]] = i;
      orig2heap[heap2orig[j]] = j;
    };

    virtual void heapify(int i) = 0;
    virtual void changeElement(int i) = 0;

  private:
    int *heap2orig;
    int *orig2heap;

    T *elems;
    int capacity;
    int nElements;



};

#endif
