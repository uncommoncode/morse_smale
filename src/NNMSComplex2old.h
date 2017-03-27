//Compute persistence based on crytsal size (number of points) instaed of
//function value differences
#ifndef NNMSCOMPLEXR2_H
#define NNMSCOMPLEXR2_H

#include "Geometry.h"
#include "EuclideanMetric.h"
#include "DenseMatrix.h"
#include "DenseVector.h"
#include "Linalg.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <set>

template<typename TPrecision>
class NNMSComplexR2{

  private:
    typedef std::pair<int, int> pair_i;
    
    typedef std::set<int> set_i;
    typedef set_i::iterator set_i_it;
    
    typedef std::set<pair_i> set_pi;
    typedef set_pi::iterator set_pi_it;

    typedef std::map< pair_i, int> map_pi_i;
    typedef map_pi_i::iterator map_pi_i_it;

    typedef std::map< int, set_i> map_i_si;
    typedef map_i_si::iterator map_i_si_it;

    typedef std::map< pair_i, double> map_pi_d;
    typedef map_pi_d::iterator map_pi_d_it;

    typedef std::map< pair_i, set_i > map_pi_si;
    typedef map_pi_si::iterator map_pi_si_it;
    

    //Steepest ascending KNNG(0,) and descending KNNG(1, ) neighbros for each point    
    FortranLinalg::DenseMatrix<int> KNNG;

    //Data points
    FortranLinalg::DenseMatrix<TPrecision> X;
    FortranLinalg::DenseVector<TPrecision> y;
    //extrema ID for ach point --- max extrema(0, ) and min extrema(1, ) 
    FortranLinalg::DenseMatrix<int> extrema;

    //map of crystals as <max, min> -> set of points in crystal 
    map_pi_si crystals;
    map_pi_si *crystalHistory;
    int nLevels;
    //map of crystals to R2 based on linear fit or -1/n if number of points n
    //too small for R2 statistic
    map_pi_d r2;
    
    //connections of maxima (to minima) and minima (to maxima)
    map_i_si connections;

    FortranLinalg::DenseVector<TPrecision> persistence;
    
    
    //extrema ID to index into X
    FortranLinalg::DenseVector<int> extremaIndex;
    
    //number of maxima, first nMax entries in extremaIndex are maxima
    int nMax;


    EuclideanMetric<TPrecision> l2;


    int ascending(int index){
      return KNNG(0, index);
    };

    int descending(int index){
      return KNNG(1, index);
    };




   double fitLM(set_i &points){
     using namespace FortranLinalg;
      //TODO decide based on power of R2 statistics
      int n = points.size();
      if(n < X.M()+2 ){
        return  - 100000.0 / n;
      }
      else{
        DenseMatrix<TPrecision> A(n, X.M() + 1);
        DenseMatrix<TPrecision> b(n, 1);
        int index = 0;
        double mean = 0;
        for(set_i_it it = points.begin(); it != points.end(); ++it, ++index){
          A(index, 0 ) = 1;
          int xindex = *it;
          for(int i=0; i<X.M(); i++){
            A(index, i+1) = X(i, xindex);
          }
          b(index, 0) = y( xindex );
          mean += y(xindex);
        }
        mean/=n;
        double sse = 0; 
        DenseMatrix<TPrecision> coef = Linalg<TPrecision>::LeastSquares(A, b, &sse);
        double sst = 0;
	for(set_i_it it = points.begin(); it != points.end(); ++it, ++index){
          int xindex = *it;
          double tmp = y(xindex) - mean;
          sst += tmp*tmp;
        }
        A.deallocate();
        b.deallocate();
        coef.deallocate();
        //adjusted R^2 weighted by crystal size
        return ( 1-(sse/sst)*(n-1.0)/(n-X.M()-1.0) )  * n / (double) X.N();
      } 
   };



   void doMerge(pair_i m){  
      int e1 = m.first;
      int e2 = m.second;
      set_i &cons1 = connections[e1];
      set_i &cons2 = connections[e2];
      for(set_i_it it = cons1.begin(); it != cons1.end(); ++it){
        int em = *it;
        pair_i p1;
        pair_i p2;
        if(e1 < nMax){
          p1 = pair_i(e1, em);
          p2 = pair_i(e2, em);
        }
        else{
          p1 = pair_i(em, e1);
          p2 = pair_i(em, e2);
        }
        map_pi_si_it p1it = crystals.find(p1);
        //is there something to merge?
        if( p1it != crystals.end()){ 
          map_pi_si_it p2it = crystals.find(p2);
          if(p2it != crystals.end() ){         
            set_i &s1 = p1it->second;
            set_i &s2 = p2it->second;
            s2.insert(s1.begin(), s1.end());
            double r = fitLM(s2);
            r2[p2] = r;
          }
          else{
            crystals[p2] = crystals[p1];
            r2[p2] = r2[p1];
          }
          r2.erase(p1);
          crystals.erase(p1);
        }
      }
      cons2.insert( cons1.begin(), cons1.end() );

      //remove references to e1
      connections.erase(e1);
      for(map_i_si_it it = connections.begin(); it != connections.end(); ++it){
        set_i &c = it->second;
        if( c.erase(e1) == 1){
          c.insert(e2);
        }
      } 
   };


   double testMerge(pair_i m, double curR2){
      //curR2 = curR2*crystals.size(); 
      int nmerges = 0;
      int e1 = m.first;
      int e2 = m.second;
      double improve = 0;
      for(set_i_it it = connections[e1].begin(); it != connections[e1].end(); ++it){
        int em = *it;
        pair_i p1;
        pair_i p2;
        if(e1 < nMax){
          p1 = pair_i(e1, em);
          p2 = pair_i(e2, em);
        }
        else{
          p1 = pair_i(em, e1);
          p2 = pair_i(em, e2);
        }
        map_pi_si_it p1it = crystals.find(p1);
        map_pi_si_it p2it = crystals.find(p2);
        //is there something to merge?
        if( p2it != crystals.end() && p1it != crystals.end()){
          ++nmerges;
          set_i &s2 = p2it->second;
          set_i &s1 = p1it->second;
          set_i sa = s1;
          sa.insert(s2.begin(), s2.end());
          double r2a = fitLM(sa);
          double r21 = r2[p1];  
          double r22 = r2[p2];
          if(r2a < 0 ){
            improve -=   (r21 + r22);
          }
          else{
            improve += r2a - (r21 + r22);
          }
            //+ (X.M()+1)( 2*( log( s1.size() + s2.size()) ) - log( sa.size() ) );  
        }
      }
      return improve;// / (crystals.size() - nmerges);
   };

   

  public:

    NNMSComplexR2(FortranLinalg::DenseMatrix<TPrecision> &Xin, FortranLinalg::DenseVector<TPrecision> &yin, int
        knn, int nl = -1, bool smooth = false, double eps=0.1) : X(Xin), y(yin){
     using namespace FortranLinalg;
      nLevels = nl;
      if(knn > X.N()){
        knn = X.N();
      }
      DenseMatrix<int> KNN(knn, X.N());
      DenseMatrix<TPrecision> KNND(knn, X.N());

      //Compute nearest neighbors
      Geometry<TPrecision>::computeANN(X, KNN, KNND, eps);


      DenseVector<TPrecision> ys;
      if(smooth){
        ys = DenseVector<TPrecision>(y.N());
        for(int i=0; i< ys.N(); i++){
          ys(i) = 0;
          for(int k=0; k<knn; k++){
            ys(i) += y(KNN(k, i));
          }
          ys(i) /= knn;//*(knn+1)/2;
        }
      }
      else{
        ys = y;
      }

      KNNG = DenseMatrix<int>(2, X.N());
      Linalg<int>::Set(KNNG, -1);
      DenseMatrix<TPrecision> G = DenseMatrix<TPrecision>(2, X.N());
      Linalg<TPrecision>::Zero(G);

      //compute steepest asc/descending neighbors
      for(int i=0; i<X.N(); i++){
        for(int k=1; k<KNN.M(); k++){
          int j = KNN(k, i);
          TPrecision g = (ys(j) - ys(i)) / sqrt(KNND(k, i));
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
      G.deallocate();
      KNND.deallocate();
      if(smooth){
       ys.deallocate();
      }


      //compute for each point its minimum and maximum based on
      //steepest ascent/descent
      extrema = DenseMatrix<int>(2, X.N()); 
      Linalg<int>::Set(extrema, -1);

      std::list<int> extremaL;
      std::list<int> path;
      int nExt = 0;
      nMax = 0;
      std::vector<int> nPnts;
      for(int e=0; e<2; e++){
        for(int i=0; i<extrema.N(); i++){
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

      extremaIndex = DenseVector<int>(nExt);
      int index = 0;
      for(std::list<int>::iterator it = extremaL.begin(); it != extremaL.end(); ++it, ++index){
        extremaIndex(index) = *it;
      }

      //Persistence based on linear fit
      //for each crystal fit linear model - if not enough points for r2
      //statistic set to 1/n with n the number of crystals in the peak. This
      //merges smallest crystals first and then based on linear fit.
      
      //create lowest persistence point assignments
      //store extrema connections min->max and max -> min
      for(int i=0; i<extrema.N(); i++){
        int e1 = extrema(0, i);
        int e2 = extrema(1, i);
        connections[e1].insert(e2);
        connections[e2].insert(e1);
        pair_i id(e1, e2);
        crystals[id].insert(i);
      }

      /*for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it){
        pair_i p = it->first;
        set_i s = it->second;
        s.insert(extremaIndex(p.first));
        s.insert(extremaIndex(p.second));
      }*/
     
      /* 
      for(int i=0; i < extrema.N(); i++){
        int max1 = extrema(0, i);
        int min1 = extrema(1, i);
        for(int k=1; k < KNN.M(); k++){
          int max2 = extrema(0, KNN(k, i));
          int min2 = extrema(1, KNN(k, i));
          if(min1 != min2 || max1 != max2){
            connections[min1].insert(max2);
            connections[min2].insert(max1);
            connections[max1].insert(min2);
            connections[max2].insert(min1);
          }
        }
      }*/

      //look for disconnected crystals
      for(map_i_si_it it =  connections.begin(); it != connections.end(); ++it){
        set_i &s = it->second;
        if(s.size() < 2){
          int e2 = *s.begin();
          if(connections[e2].size() < 2){
            //std::cout << "Disconnected barf" << std::endl;
          }
        }
      }
      
      


      //compute inital r2 of crystals
      for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it){
        set_i points = it->second;
        double r = fitLM(points);
        r2[it->first] = r;
      }
      

     
      //check possible merges (adjacent mins and maxs) for reduction in r2
      
      //stire history from 1 through nLevels
      //if nLevels < 0 store best one
      if(nLevels < 0){
        crystalHistory = new map_pi_si[1];
      }
      else{
        nLevels = std::min((int)connections.size()-1, nLevels);
        crystalHistory = new map_pi_si[nLevels];
      }

      //Store Average r2 at each level
      persistence = DenseVector<TPrecision>(connections.size()-1);
      Linalg<TPrecision>::Set(persistence,
          std::numeric_limits<TPrecision>::infinity());

      for(;;){
        double curR2 = 0;
        for(map_pi_d_it it = r2.begin(); it != r2.end(); ++it){
          curR2 += it->second;
        }
        //curR2 /= r2.size();
        persistence(persistence.N() - connections.size() + 1) = curR2;
        
 
        //std::cout << "#Crystals: " << crystals.size() << std::endl;
        //std::cout << "#Extrema: " << connections.size() << std::endl;
        //std::cout << "R2: " << curR2 << std::endl;

        //find best possible merge in terms of R2
        double best = -std::numeric_limits<double>::infinity();
        set_pi merged;
        pair_i bestMerge;
        for(map_i_si_it it = connections.begin(); it != connections.end(); ++it){
          set_i &adj = it->second;
	        for( set_i_it i1 = adj.begin(); i1 != adj.end(); ++i1){
            set_i_it i2 = i1;
            ++i2;
            for(; i2 != adj.end(); ++i2){
              pair_i m( *i1, *i2);
              std::pair<set_pi_it, bool> insert = merged.insert(m);
              if( insert.second ){
                double tmp = testMerge(m, curR2);
                if(tmp > best){
                  best = tmp;
             	    bestMerge = m;
                }
              } 
            } 
          }
        }

        //Is there an R2 reducing merge?
        if(nLevels < 0){
           if(best <= 0 ){
             crystalHistory[0] = crystals;
             break;       
           }
        }
        else if(nLevels >= connections.size() - 1){
          crystalHistory[nLevels - connections.size() +1] = crystals;
        }
        int next = connections.size();
        if(next == 2){
          break;
        }
        doMerge(bestMerge);
        if(next == connections.size()){
          if(next > nLevels){
            nLevels = -1;
            crystalHistory[0] = crystals;
          }
          else{
            for(int n = next; n < nLevels; n++){
              //std::cout << n <<std::endl;
              crystalHistory[n] = crystals; 
            }
          }
          break;
        }
      }
  
      KNN.deallocate();
    };




    //Get partioning accordinng to the crystals of the MS-complex for the
    //currently set persistence level
    FortranLinalg::DenseVector<int> getPartitions(){
     using namespace FortranLinalg;
      DenseVector<int> crys(X.N());
      getPartitions(crys);
      return crys;
    };


    void setLevel(int level){
      if(level >= 0 && level < nLevels){
        crystals = crystalHistory[nLevels-level-1];
      }
    };


    void getPartitions(FortranLinalg::DenseVector<int> &crys){
      int crystalIndex = 0;
      for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it,  ++crystalIndex){
        set_i &points = it->second;
        for(set_i_it sit = points.begin(); sit != points.end(); ++sit){
          crys(*sit) = crystalIndex;
        }
      }
    };



    int getNCrystals(){
      return crystals.size();
    };


    int getNAllExtrema(){
      return extremaIndex.N();
    }; 

    //return extrema indicies (first row is max, secon is min) for each crystal
    FortranLinalg::DenseMatrix<int> getExtrema(){
     using namespace FortranLinalg;
       DenseMatrix<int> e(2, crystals.size());
       getExtrema(e);
       return e;
    };


    void getExtrema(FortranLinalg::DenseMatrix<int> ce){
      int crystalIndex = 0;
      for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it, ++crystalIndex){
        pair_i p = it->first;
        ce(0, crystalIndex) = extremaIndex(p.first);
        ce(1, crystalIndex) = extremaIndex(p.second);
      }
    };



    void getMax(FortranLinalg::DenseVector<int> vmaxs){
      int crystalIndex = 0;
     for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it, ++crystalIndex){
        pair_i p = it->first;
        vmaxs(crystalIndex) = extremaIndex(p.first);
      }
    };


    void getMin(FortranLinalg::DenseVector<int> vmins){
      int crystalIndex = 0;
      for(map_pi_si_it it = crystals.begin(); it != crystals.end(); ++it, ++crystalIndex){
        pair_i p = it->first;
        vmins(crystalIndex) = extremaIndex(p.second);
      }
    };
    
    //get persistencies
    FortranLinalg::DenseVector<TPrecision> getPersistence(){
     using namespace FortranLinalg;
      DenseVector<TPrecision> pers(persistence.N());
      getPersistence(pers);
      return pers;
    };

    void getPersistence(FortranLinalg::DenseVector<TPrecision> pers){
      for(int i=0; i<persistence.N(); i++){
        pers(i) = persistence(i);
      }
    };

    void cleanup(){
      delete[] crystalHistory;
      extrema.deallocate();
      extremaIndex.deallocate();
      KNNG.deallocate();
    };


};

#endif 

