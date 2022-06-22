#include "LSH.hpp"

LSH::LSH(){}

LSH::LSH(int kk, int window, int tablesize, double delt){ //CONSTRUCTOR
    k = kk;
    tableSize = tablesize;
    w = window;
    LSH_hashTable = new vector<bucketList>[tableSize];
    v_array = new vector<double>[k]; //ARRAY WITH K V
    
    delta = delt;

    unsigned seed_x = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen_x(seed_x);

    uniform_real_distribution<> dist_x(0, delta);
    t_x = dist_x(gen_x);

    unsigned seed_y = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen_y(seed_y);

    uniform_real_distribution<> dist_y(0, delta);
    t_y = dist_y(gen_y);
    
}

LSH::~LSH(){}
        
void LSH::hashPush(data* p){ //PUSH POINTS IN HASH TABLE
    bucketList node;
    node.p = p;

    LSH_hashTable[gp(&node)].push_back(node); //PUSH THE NODE IN BUCKET G
}

int LSH::hi(data* p, int w, int i){ //CALCULATE Hi

    int inner_p = inner_product(p, i); //INNER PRODUCT P*V

    int t = rand() % w; //t E[0,W]

    int ans = floor((inner_p + t) / w);//FINAL Hi

    return ans;
}

void LSH::v_calculator(int d){//CULCULATE V FOR INNER PRODUCT

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
      
    normal_distribution<double> di(0, 1); 

    for(int i = 0; i < k; i++) //CALCULATE K V WITH NORMAL DISTRIBUTION FOR EACH Hi HASH FUNCTION
        for(int j = 0; j < d; j++) 
            v_array[i].push_back(di(gen));

}

double LSH::inner_product(data* p, int i){//INNER PRODUCT V*P
    
    double product = 0.0;
 
    for(int j = 0; j < p->d; j++)  //V*P
        product = product + v_array[i][j] * p->p_data[j];
    
    return product;
}

int LSH::gp(bucketList* node){//CALCULATE G (BUCKET NUMBER) FOR EACH POINT
    
    int r; 
    long IDp = 0;

    for(int i = 0; i < k; i++){
        srand(time(NULL));
        r = rand() % 100;

        IDp = IDp + modulo((r*hi(node->p, w, i)), 4294967291); //CALCULATE ID
      
    }

    node->ID = IDp;
    return modulo(IDp, tableSize); //RETURN G
}

int LSH::gp_cluster(data* node){ // CALCULATE G (BUCKET NUMBER) FOR CLUSTERS 
    int r; 
    long IDp = 0;

    for(int i = 0; i < k; i++){
        srand(time(NULL));
        r = rand() % INT32_MAX + 1;

        IDp = IDp + modulo((r*hi(node, w, i)), 4294967291);//CALCULATE ID
      
    }
    
    return modulo(IDp, tableSize);//RETURN G
}

long LSH::modulo(int rh, long M){//MODULO
    
    if (rh > 0) //r*h POSITIVE
        return rh % M;

    else{//r*h NEGATIVE
        int rem = rh % M;
        if(rem == 0)
            return 0;
        return rem + M;
    }
   
}

vector<data> LSH::kNNsearch(vector<double> q, int d, int N){ //FIND N NEAREST NEIGHBORS

    vector<data> kNN_arr; //ARRAY FOR THE N NEAREST NEIGHBORS 
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY POINT NODE
    node_temp.p = &temp;

    current_q_bucket = gp(&node_temp);//G FOR THE CURRENT QUERY POINT 
    current_ID = node_temp.ID;
   
    int i = 0;

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EACH POINT IN THE CURRENT BUCKET

        if (v.ID == current_ID){ //IF IT HAS THE SAME ID PUSH THE POINT IN KNN ARRAY

            v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);
            kNN_arr.push_back(*v.p);
            i++;

        }
    }

    if(i >= N){// IF ALL THE N HAVE THE SAME ID
        sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

        int size = kNN_arr.size();
        for(int j = N; j < size; j++)//RETUR THE KNN ARRAY WITH ONLY N POINTS
            kNN_arr.pop_back();
      
        return kNN_arr;
    }

    kNN_arr.clear();

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EIASH POINT IN THE CURRENT BUCKET 
        v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);//CALCULATE DISTANCE 
        kNN_arr.push_back(*v.p);
    }

    sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

    int size = kNN_arr.size();
    for(int j = N; j < size; j++)
        kNN_arr.pop_back();
    
    return kNN_arr;
    
}

vector<data*> LSH::range_search(vector<double> q, int d, int minR, int maxR){ //RANGE SEARCH
    vector<data*> R_arr;//ARRAY FOR THE POINTS WITH DISTANCE >MINR AND <= MAXR
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY POINT NODE
    node_temp.p = &temp;

    for(auto v = LSH_hashTable[current_q_bucket].begin(); v != LSH_hashTable[current_q_bucket].end(); v++){//FOR EACH POINT IN THE CURRENT BUCKET 
 
        v->p->current_distance = distance(v->p->p_data, temp.p_data, temp.d); //CALCULATE DISTANCE
       
        if(v->p->current_distance > minR && v->p->current_distance <= maxR) //MINR> CURRENT DISTANCE <= MAXARR
            R_arr.push_back(v->p);
        
    }

    return R_arr;

}

void LSH::print( ){
    for(int i=0; i< 4; i++){
        cout<<"///////////////////////////////////////////////////////"<<endl;
        for(int j=0; j< LSH_hashTable[i].size(); j++)
            cout<<LSH_hashTable[i][j].kati<<" ";
        cout<<endl;
    }
}