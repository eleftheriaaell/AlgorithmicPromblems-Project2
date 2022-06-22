#include "LSH.hpp"

/* Discrete Frechet LSH --------------------------------------------------------------------------------------------------------------------*/

void LSH::DiscreteFrechet_hashPush(data* p){ //PUSH POINTS IN HASH TABLE

    bucketList node;
    node.p = p; 
    DiscreteFrechet_create_new_curve(&node);

    LSH_hashTable[DiscreteFrechet_gp(&node)].push_back(node); //PUSH THE NODE IN BUCKET G
}

void LSH::DiscreteFrechet_create_new_curve(bucketList* node){//CREATE NEW 2d CURVE AFTER GRID FOR HASHING

    int current_size = node->p->p_data.size();
    
    for(int i = 0; i < current_size; i++){

        double qx = floor(( abs((i+1) - t_x) )/delta + 0.5) * delta + t_x;//GRID FOR X
        node->new_curve.push_back(qx);

        double piy = node->p->p_data[i];
        double qy = floor(abs(piy - t_y)/delta + 0.5) * delta + t_y;//GRID FOR Y
        node->new_curve.push_back(qy);

    }

    int current_element_x = 2;
    int current_element_y = 3;

    while(current_element_x < node->new_curve.size()){//ERASE DUPLICATES

        //IF CURRENT (X,Y) IS THE SAME AS THE PREVIOUS (X,Y) ERASE THE CURRENT
        if(node->new_curve[current_element_x] == node->new_curve[current_element_x - 2] && node->new_curve[current_element_y] == node->new_curve[current_element_y - 2]){
            node->new_curve.erase(node->new_curve.begin() + current_element_x); //ERASE X
            node->new_curve.erase(node->new_curve.begin() + current_element_x); //ERASE Y
        }else{
            current_element_x += 2;
            current_element_y += 2;
        }
    }

    int curve_size = node->new_curve.size();
    int d_size = 2 * node->p->d;

    if(curve_size < d_size)//PADDING
        for(int i = curve_size; i < d_size; i++){
            node->new_curve.push_back(33554432.0);
        }

}

int LSH::DiscreteFrechet_hi(bucketList* node, int w, int i){ //CALCULATE Hi

    int inner_p = DiscreteFrechet_inner_product(node, i); //INNER PRODUCT P*V

    int t = rand() % w; //t E[0,W]

    int ans = floor((inner_p + t) / w);//FINAL Hi

    return ans;
}

double LSH::DiscreteFrechet_inner_product(bucketList* node, int i){//INNER PRODUCT V*P
    
    double product = 0.0;
    
    for(int j = 0; j < 2 * node->p->d; j++)  //V*P
        product = product + v_array[i][j] * node->new_curve[j];
    
    return product;
}

int LSH::DiscreteFrechet_gp(bucketList* node){//CALCULATE G (BUCKET NUMBER) FOR EACH POINT
    
    int r; 
    long IDp = 0;

    for(int i = 0; i < k; i++){
        srand(time(NULL));
        r = rand() % 100;

        IDp = IDp + modulo((r*DiscreteFrechet_hi(node, w, i)), 4294967291); //CALCULATE ID
        
    }

    node->ID = IDp;

    return modulo(IDp, tableSize); //RETURN G
}

int LSH::DiscreteFrechet_gp_cluster(data* node){ // CALCULATE G (BUCKET NUMBER) FOR CLUSTERS 

    int r; 
    long IDp = 0;

    bucketList node_temp;//QUERY POINT NODE
    node_temp.p = node;
    node_temp.p->d = node->p_data.size();

    DiscreteFrechet_create_new_curve(&node_temp);//CREATE THE GRID CURVE FOR HASHING

    for(int i = 0; i < k; i++){
        srand(time(NULL));
        r = rand() % INT32_MAX + 1;

        IDp = IDp + modulo((r*DiscreteFrechet_hi(&node_temp, w, i)), 4294967291);//CALCULATE ID
      
    }
    
    return modulo(IDp, tableSize);//RETURN G
}

vector<data> LSH::DiscreteFrechet_kNNsearch(vector<double> q, int d, int N){ //FIND N NEAREST NEIGHBORS

    vector<data> kNN_arr; //ARRAY FOR THE N NEAREST NEIGHBORS 
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY CURVE NODE
    node_temp.p = &temp;

    DiscreteFrechet_create_new_curve(&node_temp);//CREATE THE GRID CURVE FOR THE QUERY
   
    current_q_bucket = DiscreteFrechet_gp(&node_temp);//G FOR THE CURRENT QUERY CURVE 
    current_ID = node_temp.ID;

    int i = 0;

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EACH CURVE IN THE CURRENT BUCKET

        if (v.ID == current_ID){ //IF IT HAS THE SAME ID PUSH THE CURVE IN KNN ARRAY

            v.p->current_distance = Discrete_Frechet_Distance(v.p->p_data, temp.p_data);
            kNN_arr.push_back(*v.p);
            i++;

        }
    }

    if(i >= N){// IF ALL THE N HAVE THE SAME ID
        sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

        int size = kNN_arr.size();
        for(int j = N; j < size; j++)//RETUR THE KNN ARRAY WITH ONLY N CURVES
            kNN_arr.pop_back();

        return kNN_arr;
    }

    kNN_arr.clear();

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EACH CURVE IN THE CURRENT BUCKET 
        v.p->current_distance = Discrete_Frechet_Distance(v.p->p_data, temp.p_data);//CALCULATE DISTANCE 
        kNN_arr.push_back(*v.p);
    }

    sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

    int size = kNN_arr.size();
    for(int j = N; j < size; j++)
        kNN_arr.pop_back();
      

    return kNN_arr;
    
}

vector<data*> LSH::DiscreteFrechet_range_search(vector<double> q, int d, int minR, int maxR){ //RANGE SEARCH
    vector<data*> R_arr;//ARRAY FOR THE CURVES WITH DISTANCE >MINR AND <= MAXR
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY CURVE NODE
    node_temp.p = &temp;

    for(auto v = LSH_hashTable[current_q_bucket].begin(); v != LSH_hashTable[current_q_bucket].end(); v++){//FOR EACH CURVE IN THE CURRENT BUCKET 
 
        v->p->current_distance = Discrete_Frechet_Distance(v->p->p_data, temp.p_data); //CALCULATE DISTANCE
       
        if(v->p->current_distance > minR && v->p->current_distance <= maxR) //MINR> CURRENT DISTANCE <= MAXARR
            R_arr.push_back(v->p);
        
    }

    return R_arr;

}

/* Continuous Frechet LSH ----------------------------------------------------------------------------------------------------------------------*/

void LSH::ContinuousFrechet_hashPush(data* p){ //PUSH CURVE IN HASH TABLE

    bucketList node;
    node.p = p;
    ContinuousFrechet_create_new_curve(&node);

    LSH_hashTable[ContinuousFrechet_gp(&node)].push_back(node); //PUSH THE NODE IN BUCKET G
}

void LSH::ContinuousFrechet_create_new_curve(bucketList* node){ //CREATE NEW CURVE AFTER GRID FOR HASHING

    int current_size = node->p->p_data.size();

    for(int i = 0; i < current_size; i++)
       node->new_curve.push_back(floor( (node->p->p_data[i] + t_x)/delta) * delta);//GRID

    
    int current_element = 1;

    while(current_element < node->new_curve.size() - 1){//MINIMA AND MAXIMA

        double a = node->new_curve[current_element - 1];
        double b = node->new_curve[current_element];
        double c = node->new_curve[current_element + 1];

        if(a < c){
            if(a <= b && b <= c)
                node->new_curve.erase(node->new_curve.begin() + current_element);
            
        }else{ 
            if(c <= b && b <= a)
                node->new_curve.erase(node->new_curve.begin() + current_element);
        }

        current_element++;
        
    }

    int curve_size = node->new_curve.size();
    int d_size = node->p->d;

    if(curve_size < d_size)//PADDING
        for(int i = curve_size; i < d_size; i++){
            node->new_curve.push_back(33554432.0);
        }

}

int LSH::ContinuousFrechet_hi(bucketList* node, int w, int i){ //CALCULATE Hi

    int inner_p = ContinuousFrechet_inner_product(node, i); //INNER PRODUCT P*V

    int t = rand() % w; //t E[0,W]

    int ans = floor((inner_p + t) / w);//FINAL Hi

    return ans;
}

double LSH::ContinuousFrechet_inner_product(bucketList* node, int i){//INNER PRODUCT V*P
    
    double product = 0.0;
 
    for(int j = 0; j < node->p->d; j++)  //V*P
        product = product + v_array[i][j] * node->new_curve[j];
    
    return product;
}

int LSH::ContinuousFrechet_gp(bucketList* node){//CALCULATE G (BUCKET NUMBER) FOR EACH CURVE
    
    int r; 
    long IDp = 0;

    for(int i = 0; i < k; i++){
        srand(time(NULL));
        r = rand() % INT32_MAX + 1;

        IDp = IDp + modulo((r*ContinuousFrechet_hi(node, w, i)), 4294967291); //CALCULATE ID
        
    }

    node->ID = IDp;
    return modulo(IDp, tableSize); //RETURN G
}


vector<data> LSH::ContinuousFrechet_kNNsearch(vector<double> q, int d, int N){ //FIND N NEAREST NEIGHBORS
    
    vector<data> kNN_arr; //ARRAY FOR THE N NEAREST NEIGHBORS 
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY CURVE NODE
    node_temp.p = &temp;

    ContinuousFrechet_create_new_curve(&node_temp);//CREATE GRID CURVE FOR HASHING
    
    current_q_bucket = ContinuousFrechet_gp(&node_temp);//G FOR THE CURRENT QUERY CURVE 
    current_ID = node_temp.ID;
  
    int i = 0;

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EACH CURVE IN THE CURRENT BUCKET

        if (v.ID == current_ID){ //IF IT HAS THE SAME ID PUSH THE CURVE IN KNN ARRAY
            v.p->current_distance = ContinuousFrechet_distance(v.p->p_data, temp.p_data);
            kNN_arr.push_back(*v.p);
            i++;

        }
    }

    if(i >= N){// IF ALL THE N HAVE THE SAME ID
        sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

        int size = kNN_arr.size();
        for(int j = N; j < size; j++)//RETUR THE KNN ARRAY WITH ONLY N CURVES
            kNN_arr.pop_back();
      
        return kNN_arr;
    }

    kNN_arr.clear();

    for(auto v: LSH_hashTable[current_q_bucket]){//FOR EIASH CURVE IN THE CURRENT BUCKET 
        v.p->current_distance = ContinuousFrechet_distance(v.p->p_data, temp.p_data);//CALCULATE DISTANCE 
        kNN_arr.push_back(*v.p);
    }

    sort(kNN_arr.begin(), kNN_arr.end(), comp);//SORT N NEAREST NEIGHBORS

    int size = kNN_arr.size();
    for(int j = N; j < size; j++)
        kNN_arr.pop_back();
    
    return kNN_arr;
    
}

