//Morse-Smale complex computation as described in:
// Gerber S, Bremer PT, Pascucci V, Whitaker R (2010). 
// “Visual Exploration of High Dimensional Scalar Functions.” 
// IEEE Transactions on Visualization and Computer Graphics, 16(6), 1271–1280.

//author: Samuel Gerber

#ifndef NNMSCOMPLEX_H
#define NNMSCOMPLEX_H

#include "Geometry.h"
#include "EuclideanMetric.h"
#include "SquaredEuclideanMetric.h"
#include "DenseMatrix.h"
#include "DenseVector.h"
#include "Linalg.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>



template<typename TPrecision>
class NNMSComplex{

  private:
    typedef std::map< std::pair<int, int>, int> map_pi_i;
    typedef map_pi_i::iterator map_pi_i_it;

    typedef std::multimap< int, std::pair<int, int> > mmap_i_pi;
    typedef mmap_i_pi::iterator mmap_i_pi_it;
    
    typedef typename std::map< std::pair<int, int>, TPrecision > map_pi_f;
    typedef typename map_pi_f::iterator map_pi_f_it;
    
    typedef typename std::map< TPrecision, std::pair<int, int>, std::less<TPrecision>  > map_f_pi;
    typedef typename map_f_pi::iterator map_f_pi_it;


    //Steepest ascending KNNG(0,) and descending KNNG(1, ) neighbors for each point    
    FortranLinalg::DenseMatrix<int> KNNG;

    //Data points
    FortranLinalg::DenseMatrix<TPrecision> X;
    FortranLinalg::DenseVector<TPrecision> y;
      
    FortranLinalg::DenseMatrix<int> KNN;
    FortranLinalg::DenseMatrix<TPrecision> KNND;
    

    //extrema ID for ach point --- max extrema(0, ) and min extrema(1, ) 
    FortranLinalg::DenseMatrix<int> extrema;

    //map of crystals as <max, min> -> crystal ID    
    map_pi_i crystals;
    //map after merging for reduced set of crystals
    map_pi_i pcrystals;
    
    
    //inital persistencies for each crystal
    map_f_pi persistence;

    //extrema ID to index into X
    FortranLinalg::DenseVector<int> extremaIndex;
    //extrema after merging of crystals e.g extrema(i) -> merge(extrema(i))
    FortranLinalg::DenseVector<int> merge;
    //number of maxima, first nMax entries in extremaIndex are maxima
    int nMax;


    EuclideanMetric<TPrecision> l2;


    int ascending(int index){
      return KNNG(0, index);
    };

    int descending(int index){
      return KNNG(1, index);
    };

    int followChain(int i){
      while(merge(i) != i){
        i = merge(i);
      }
      return i;
    };


    //compute crystals based on merge chain
    void mergeCrystals(){
      //initalize simplified crystals to original crystal assignments
      int nCrystals = 0;
      //reassign crystals based on merge chain
      pcrystals.clear();
      for(map_pi_i_it it = crystals.begin(); it != crystals.end(); ++it){
        std::pair<int, int> p = (*it).first;
        
        
        //follow merge chains for min and max
        p.first = merge(p.first);
        p.second = merge(p.second);
        
        //check if we created a new crystal otherwise assign to existing
        //crystal
        map_pi_i_it ito = pcrystals.find(p);
        if(ito == pcrystals.end()){
          pcrystals[p] = nCrystals;
          nCrystals++;
        }
      }
   };

  public:


    NNMSComplex(FortranLinalg::DenseMatrix<TPrecision> &Xin, FortranLinalg::DenseVector<TPrecision> &yin,
                FortranLinalg::DenseMatrix<int> &KNNin, FortranLinalg::DenseMatrix<TPrecision> &KNNDin, 
                bool smooth = false, double sigma2=0) : X(Xin), y(yin), KNN(KNNin), KNND(KNNDin){
      runMS(smooth, sigma2);

    };

 
    NNMSComplex(FortranLinalg::DenseMatrix<TPrecision> &Xin, FortranLinalg::DenseVector<TPrecision> &yin, int
        knn, bool smooth = false, double eps=0.01, double sigma2=0) : X(Xin), y(yin){
      if(knn > (int) X.N()){
        knn = X.N();
      }
      KNN = FortranLinalg::DenseMatrix<int>(knn, X.N());
      KNND = FortranLinalg::DenseMatrix<TPrecision>(knn, X.N());

      //Compute nearest nieghbors
      //Geometry<TPrecision>::computeANN(X, KNN, KNND, eps);
      SquaredEuclideanMetric<TPrecision> dist;
      Geometry<TPrecision>::computeKNN(X, KNN, KNND, dist);
      
      runMS(smooth, sigma2);      
      KNN.deallocate();
      KNND.deallocate();
    };







    //Compute the MS crystals for the given persistence level. Neighboring
    //extrema with a absolute difference between saddle and lower exterma
    //smaller than pLevel, are recursively joined into a single extrema.
    void mergePersistence(TPrecision pLevel){
      //compute merge chain 
      for(unsigned int i=0; i<merge.N(); i++){
        merge(i) = i;
      } 
      

  
      for(map_f_pi_it it = persistence.begin(); it != persistence.end() && (*it).first < pLevel; ++it){
        std::pair<int, int> p = (*it).second;        
        p.first = followChain(p.first);
        p.second = followChain(p.second);
        if(p.first < nMax){
          if( y(extremaIndex(p.first)) > y(extremaIndex(p.second)) ){ 
            std::swap(p.second, p.first); 
          }
        }
        else{
          if( y(extremaIndex(p.first)) < y(extremaIndex(p.second)) ){ 
            std::swap(p.second, p.first);
          }
        }
        merge(p.first) = p.second;
      }
      for(unsigned int i=0; i<merge.N(); i++){
        merge(i) = followChain(i);
      }


      //compute crystals based on merge chain
      mergeCrystals();
      
    };

   


    //Get partioning accordinng to the crystals of the MS-complex for the
    //currently set persistence level
    FortranLinalg::DenseVector<int> getPartitions(){
      FortranLinalg::DenseVector<int> crys(X.N());
      getPartitions(crys);
      return crys;
    };




    void getPartitions(FortranLinalg::DenseVector<int> &crys){
      for(unsigned int i=0; i<X.N(); i++){
        std::pair<int, int> p( merge(extrema(0, i)), merge(extrema(1, i)) );
        crys(i) = pcrystals[p];
      }
    };


    int getNCrystals(){
      return pcrystals.size();
    };


    int getNAllExtrema(){
      return extremaIndex.N();
    }; 

    //return extrema indicies (first row is max, secon is min) for each crystal
    FortranLinalg::DenseMatrix<int> getCrystals(){
       FortranLinalg::DenseMatrix<int> e(2, pcrystals.size());
       getCrystals(e);
       return e;
    };


    void getCrystals(FortranLinalg::DenseMatrix<int> ce){
     for(map_pi_i_it it = pcrystals.begin(); it != pcrystals.end(); ++it){
        std::pair<int, int> p = (*it).first;
        ce(0, (*it).second) = extremaIndex(p.first);
        ce(1, (*it).second) = extremaIndex(p.second);
      }
    };



    void getMax(FortranLinalg::DenseVector<int> vmaxs){
     for(map_pi_i_it it = pcrystals.begin(); it != pcrystals.end(); ++it){
        std::pair<int, int> p = (*it).first;
        vmaxs((*it).second) = extremaIndex(p.first);
      }
    };



    void getMin(FortranLinalg::DenseVector<int> vmins){
     for(map_pi_i_it it = pcrystals.begin(); it != pcrystals.end(); ++it){
        std::pair<int, int> p = (*it).first;
        vmins((*it).second) = extremaIndex(p.second);
      }
    };


    //get persistencies
    FortranLinalg::DenseVector<TPrecision> getPersistence(){
      FortranLinalg::DenseVector<TPrecision> pers(persistence.size()+1);
      getPersistence(pers);
      return pers;
    };


    void getPersistence(FortranLinalg::DenseVector<TPrecision> pers){
      int index = 0;
      for(map_f_pi_it it = persistence.begin(); it != persistence.end(); ++it, ++index){
        pers(index) = (*it).first;
      }
      pers(index) = std::numeric_limits<TPrecision>::max();
    };



    void cleanup(){
      extrema.deallocate();
      merge.deallocate();
      extremaIndex.deallocate();
      KNNG.deallocate();

    };





private:



    void runMS(bool smooth, double sigma2){


      int knn = KNN.M();

      FortranLinalg::DenseVector<TPrecision> ys;
      if(smooth){
        ys = FortranLinalg::DenseVector<TPrecision>(y.N());
        for(unsigned int i=0; i< ys.N(); i++){
          ys(i) = 0;
          double wsum = 0;
          for(int k=0; k<knn; k++){
            double w = exp( -KNND(k, i) / sigma2 );
            ys(i) += w*y(KNN(k, i));
            wsum += w;
          }
          ys(i) /= wsum;//*(knn+1)/2;
        }
        //y.deallocate();
        y = ys;
      }
      else{
        ys = y;
      }


      KNNG = FortranLinalg::DenseMatrix<int>(2, X.N());
      FortranLinalg::Linalg<int>::Set(KNNG, -1);
      FortranLinalg::DenseMatrix<TPrecision> G = FortranLinalg::DenseMatrix<TPrecision>(2, X.N());
      FortranLinalg::Linalg<TPrecision>::Zero(G);


      //compute steepest asc/descending neighbors
      for(unsigned int i=0; i<X.N(); i++){
        for(unsigned int k=1; k<KNN.M(); k++){
          int j = KNN(k, i);
          double d = sqrt(KNND(k, i));
          double g = ys(j) - ys(i);
          if(d == 0 ){
            /*if(g > 0){
              g = std::numeric_limits<double>::max(); 
            }
            else{
              g = std::numeric_limits<double>::min(); 
            }*/
            g = 0;
          }
          else{
            g = g / d; 
          }
          
          if(G(0, i) < g){
            G(0, i) = g;
            KNNG(0, i) = j;
          }
          else if(G(1, i) > g){
            G(1, i) = g;
            KNNG(1, i) = j;
          }          
          if(G(0, j) < -g){
            G(0, j) = -g;
            KNNG(0, j) = i;
          }
          else if(G(1, j) > -g){
            G(1, j) = -g;
            KNNG(1, j) = i;
          }
        }
      }
      //if(smooth){
      // ys.deallocate();
      //}



      //compute for each point its minimum and maximum based on
      //steepest ascent/descent
      extrema = FortranLinalg::DenseMatrix<int>(2, X.N()); 
      FortranLinalg::Linalg<int>::Set(extrema, -1);

      std::list<int> extremaL;
      std::list<int> path;
      int nExt = 0;
      nMax = 0;
      std::vector<int> nPnts;
      for(int e=0; e<2; e++){
        for(unsigned int i=0; i<extrema.N(); i++){
          
/*
          int current = i;
          int prev = current;
          while(current != -1) {
            prev = current;
            if(e==0){
              current = ascending(prev);
            }
            else{
              current = descending(prev);
            }
          }
          int ind = 0;
          for(std::list<int>::iterator it = extremaL.begin();  it!=extremaL.end(); ++it){
            if((*it) == prev){
              break; 
            }
            ++ind;
          }      
          extrema(e, i) = ind;
          if(ind == extremaL.size()){
            extremaL.push_back(prev);  
            nExt++;
            if(e==0){
              nMax++;
            }
          }

*/

          if(extrema(e, i) == -1){
            path.clear();
            int prev = i;
            while(prev != -1 && extrema(e, prev) == -1){
              path.push_back(prev);
              if(e==0){
                prev = ascending(prev);
              }
              else{
                prev = descending(prev);
              }
            }
            int ext = -1;
            if(prev == -1){
              int extIndex = path.back();
              extremaL.push_back(extIndex);
              ext = nExt;
              nExt++;
              if(e==0){
                nMax++;
              }
            }
            else{
              ext = extrema(e, prev);
            }
            for(std::list<int>::iterator it = path.begin(); it!=path.end(); ++it){
              extrema(e, *it) = ext;
            }   
          }

        }
      }


      //setup crystals for zero peristence level
      int crystalID = 0;
      for(unsigned int i=0; i<extrema.N(); i++){
        std::pair<int, int> id(extrema(0, i), extrema(1, i));
        if(crystals.find(id) == crystals.end()){
          crystals[id] = crystalID;
          ++crystalID;
        }
      }


      //Persistence

      //initalize:
      //-persistence levels: difference between saddle point and extrema of
      //neighboring crystals 
      //-merge indices: merging to extrema
      extremaIndex = FortranLinalg::DenseVector<int>(nExt);
      merge = FortranLinalg::DenseVector<int>(nExt);
      int index = 0;
      for(std::list<int>::iterator it = extremaL.begin(); it != extremaL.end(); ++it, ++index){
        extremaIndex(index) = *it;
      }

      //inital persistencies
      //store as pairs of extrema such thats p.first merges to p.second (e.g.
      //p.second is the max/min with the larger/smaller function value
      map_pi_f pinv;
      for(int e=0; e<2; e++){
        for(unsigned int i=0; i < extrema.N(); i++){
          int e1 = extrema(e, i);
          for(unsigned int k=1; k < KNN.M(); k++){
            int e2 = extrema(e, KNN(k, i));
            if(e1 != e2){
              std::pair<int, int> p;
              TPrecision pers = 0;
              if(e==0){
                if(y(extremaIndex(e1)) > y(extremaIndex(e2)) ){
                  p.first = e2;
                  p.second = e1;
                }
                else{
                  p.first = e1;
                  p.second = e2;
                }
                pers = y(extremaIndex(p.first)) -  std::min(y(i), y(KNN(k, i)));
                //pers = l2.distance(X, extremaIndex(e1), X, extremaIndex(e2) ); 
              }
              else{
                if(y(extremaIndex(e1)) < y(extremaIndex(e2)) ){
                  p.first = e2;
                  p.second = e1;
                }
                else{
                  p.first = e1;
                  p.second = e2;
                }
                pers = std::max(y(i), y(KNN(k, i))) - y(extremaIndex(p.first));
                //pers = l2.distance(X, extremaIndex(e1), X, extremaIndex(e2) ); 
              }

              map_pi_f_it it = pinv.find(p);
              if(it!=pinv.end()){
                TPrecision tmp = (*it).second;
                if(pers < tmp){
                  (*it).second = pers;
                }
              }
              else{
                pinv[p] = pers;
              }
            }
          }
        }
      }
      

      for(map_pi_f_it it = pinv.begin(); it != pinv.end(); ++it){
        persistence[(*it).second] = (*it).first;
      }


      //compute final persistencies - recursively merge smallest persistence
      //extrema and update remaining peristencies depending on the merge
      for(unsigned int i=0; i<merge.N(); i++){
        merge(i) = i;
      } 
      
      map_f_pi ptmp;
      map_pi_f pinv2;
      while(!persistence.empty()){
        map_f_pi_it it = persistence.begin(); 
        std::pair<int, int> p = (*it).second;
	
        //store old extrema merging pair and persistence
        std::pair<int, int> pold = p;
        double pers = (*it).first;

        //find new marging pair, based on possible previous merges
        //make sure that p.first is the less significant extrema as before
        p.first = followChain(p.first);
        p.second = followChain(p.second);
        if(p.first < nMax){
          if( y(extremaIndex(p.first)) > y(extremaIndex(p.second)) ){ 
            std::swap(p.second, p.first); 
          }
        }
        else{
          if( y(extremaIndex(p.first)) < y(extremaIndex(p.second)) ){ 
            std::swap(p.second, p.first);
          }
        }
        
        //remove current merge pair from list
        persistence.erase(it);

        //are the extrema already merged?
        if(p.first == p.second) continue;

        //check if there is new merge pair with increased persistence
        TPrecision diff = 0;
        if(p.first < nMax) {
          diff = y(extremaIndex(p.first)) - y(extremaIndex(pold.first)) ;
        }
        else{
          diff = y(extremaIndex(pold.first)) - y(extremaIndex(p.first)) ;
        }

        if( diff > 0  ){
          //if the persistence increased insert into the persistence list and
          //merge possible other extrema with smaller persistence values first
          double npers = pers + diff;
          persistence[npers] = p;
        }
        //otherwise merge the pair
        else{
          //check if the pair has not been previously merged
          map_pi_f_it invIt = pinv2.find(p);
          if(pinv2.end() == invIt){
            merge(p.first) = p.second;
            ptmp[pers] = p;
            pinv2[p] = pers;
          }
        }
      }
      persistence = ptmp;


      //initialize to 0 persistence
      mergePersistence(0);
  
    };


};

#endif 

