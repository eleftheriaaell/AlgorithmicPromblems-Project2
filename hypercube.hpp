#include <list>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <time.h>
#include <string>

#include "data.hpp"

using namespace std;

class HyperCube{

    private:
        int k;
        int w;
        int tableSize;
        vector<float> *v_array;//ARRAY WITH K V FOR EACH Hi
        struct bucketList{
            string bit_label;
            data* p;
        };

        vector<bucketList> *hyper_hashTable;
        string fi(int hi, int i);           //hashfunction fi

        vector<int> *zero_hi;//ARRAY THAT I STORE THE Hi VALUES WITH BIT 0
        vector<int> *one_hi;//ARRAY THAT I STORE THE Hi VALUES WITH BIT 1
    
    public:
        HyperCube();
        HyperCube(int d, int window);
        ~HyperCube();
        string label(bucketList node);
        void hashPush(data* p);
        int hi(data* p, int w, int i);
        void v_calculator(int d);
        float inner_product(data* p, int i);
        int hamming_distance(int a, int b);
        vector<data> kNNsearch(vector<double> q, int d, int N, int M, int probes);
        vector<data*> range_search(vector<double> q, int d, int minR, int maxR, int M, int probes);
        string label_cluster(data* node);
        void print();
        

};