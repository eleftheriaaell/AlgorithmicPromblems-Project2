#include <iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <time.h>
#include "data.hpp"

using namespace std;

class LSH{ // LSH CLASS

    private:
        int k; 
        int w;
        int tableSize;
        int current_ID;//ID FOR THE CURRENT QUERY POINT 
        double delta;
        double t_x;//t FOR X
        double t_y;//t FOR Y
        vector<double> *v_array;//ARRAY WITH K V FOR EACH Hi FUNCTION
        struct bucketList{//HASH TABLE NODE
            long ID;
            int kati=0;
            vector<double> new_curve;//CURVE AFTER FILTERING AND SNAPING
            data* p;
        };

        vector<bucketList> *LSH_hashTable;//HASH TABLE
        int gp(bucketList* node);//g(p)
        int DiscreteFrechet_gp(bucketList* node);//g(p) for DiscreteFrechet
        int ContinuousFrechet_gp(bucketList* node);//g(p) for DiscreteFrechet
    
    public:
        LSH();
        LSH(int kk, int window, int tablesize, double delta);
        ~LSH();

        int current_q_bucket;//G BUCKET FOR CURRENT QUERY
        
        /* LSH */
        void hashPush(data* p);
        int hi(data* p, int w, int i);
        void v_calculator(int d);
        double inner_product(data* p, int i);
        long modulo(int rh, long M);
        int gp_cluster(data* node);
        vector<data> kNNsearch(vector<double> q, int d, int N);
        vector<data*> range_search(vector<double> q, int d, int minR, int maxR);    

        /* Discrete Frechet LSH */
        void DiscreteFrechet_hashPush(data* p);
        void DiscreteFrechet_create_new_curve(bucketList* node);
        double DiscreteFrechet_inner_product(bucketList* node, int i);
        int DiscreteFrechet_hi(bucketList* node, int w, int i);
        int DiscreteFrechet_gp_cluster(data* node);
        vector<data> DiscreteFrechet_kNNsearch(vector<double> q, int d, int N);
        vector<data*> DiscreteFrechet_range_search(vector<double> q, int d, int minR, int maxR);  
        
        /* Continuous Frechet LSH */
        void ContinuousFrechet_hashPush(data* p);
        void filtering(bucketList* node);
        void ContinuousFrechet_create_new_curve(bucketList* node);
        double ContinuousFrechet_inner_product(bucketList* node, int i);
        int ContinuousFrechet_hi(bucketList* node, int w, int i);
        vector<data> ContinuousFrechet_kNNsearch(vector<double> q, int d, int N);
        void filtering(data* node);
        
        void print(); 


};