#include "hypercube.hpp"

HyperCube::HyperCube(){}

HyperCube::HyperCube(int d, int window){//CONSTRUCTOR
    k = d;
    tableSize = pow(2, k);
    w = window;
    hyper_hashTable = new vector<bucketList>[tableSize];

    v_array = new vector<float>[k];//ARRAY WITH K V FOR EACH Hi FUNCTION
    zero_hi = new vector<int>[k];//ARRAY TO STORE Hi VALUES WITH 0 BIT
    one_hi = new vector<int>[k];//ARRAY TO STORE Hi VALUES WITH 1 BIT
}

HyperCube::~HyperCube(){}

string HyperCube::label(bucketList node){//LABELING FOR EACH POINT
    int h;
    string bit_string;
    bit_string.clear();

    for(int i = 0; i < k; i++){// FOR EACH Hi FUNCTION
        h = hi(node.p, w, i);
        
        if(binary_search(zero_hi[i].begin(), zero_hi[i].end(), h)){//IF THE VALUE OF Hi EXIST WITH BIT 0
            bit_string = bit_string + '0';
        }
        else if(binary_search(one_hi[i].begin(), one_hi[i].end(), h)){//IF THE VALUE IF Hi EXIST WITH BIT 1
            bit_string = bit_string + '1';
        }
        else{// IF Hi DOESNT EXIST 
            string f = fi(h, i);
            bit_string = bit_string + f;
        }        
    }

    return bit_string;//RETURN LABEL
}

string HyperCube::label_cluster(data* node){//LABELING FOR EACH POINT IN CLUSTERING
    int h;
    string bit_string;
    bit_string.clear();

    for(int i = 0; i < k; i++){//FOR Hi FUNCTION
        h = hi(node, w, i);
        
        if(binary_search(zero_hi[i].begin(), zero_hi[i].end(), h)){//IF THE VALUE OF Hi EXIST WITH BIT 0
            bit_string = bit_string + '0';
        }
        else if(binary_search(one_hi[i].begin(), one_hi[i].end(), h)){//IF THE VALUE OF Hi EXIST WITH BIT 1
            bit_string = bit_string + '1';
        }
        else{//IF Hi DOESNT EXIST
            string f = fi(h, i);
            bit_string = bit_string + f;
        }        
    }

    return bit_string;//RETURN LABEL
}
        
void HyperCube::hashPush(data* p){//PUSH POINT IN HASH TABLE
    bucketList node;
    node.p = p; 
    
    node.bit_label = label(node);//CURRENT POINT LABEL
    int bucket_label = stoi(node.bit_label, 0, 2);//CONVER BIT LABEL INTO INT
    hyper_hashTable[bucket_label].push_back(node);//PUSH POINT IN HASH TABLE
}  

string HyperCube::fi(int hi, int i){//SET RANDOM BIT FOR Hi

    int f = rand() % 2;
    if(f == 0){
        zero_hi[i].push_back(hi);
        return "0";
    }
    else{ 
        one_hi[i].push_back(hi);
        return "1";
    }

}

int HyperCube::hi(data* p, int w, int i){//CALCULATE Hi

    int inner_p = inner_product(p, i);//INNER PRODUCT P*V

    int t = rand() % w;

    int ans = floor((inner_p + t) / w);

    return ans;//RETURN Hi VALUE
}

void HyperCube::v_calculator(int d){//CALCULATE V FOR INNER PRODUCT

    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
      
    normal_distribution<float> di(0, 1); 

    for(int i = 0; i < k; i++)//CALCULATE K V WITH NORMAL DISTRIBUTION
        for(int j = 0; j < d; j++) 
            v_array[i].push_back(di(gen));

}

float HyperCube::inner_product(data* p, int i){//INNER PRODUCT V*P
    
    float product = 0.0;
 
    for(int j = 0; j < p->d; j++)//V*P
        product = product + v_array[i][j] * p->p_data[j];
    
    return product;
}

int HyperCube::hamming_distance(int a, int b){//CALCULATE THE HAMMING DISTANCE BETWEEN 2 NUMBERS
    int x = a ^ b;
    int ham_dist = 0;

    while (x > 0) {
        ham_dist += x & 1;
        x >>= 1;
    }

    return ham_dist;
}


vector<data> HyperCube::kNNsearch(vector<double> q, int d, int N, int M, int probes){//FIND N NEAREST NEGHBORS
    
    vector<data> M_max_points;//ARRAY FOR M POINTS (MAX NUMBER OF POINTS)
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY POINT NODE
    node_temp.p = &temp;

    node_temp.bit_label = label(node_temp);//LABEL FOR CURRENT QUERY NODE
    int label = stoi(node_temp.bit_label, 0, 2);
    int q_label = label;

    random_device rd;
    mt19937 g(rd());

    int M_size = M;

    int current_ham_dist = 1;
    int current_bucket = 0;

    for(int probe_cnt = 1; probe_cnt <= probes; probe_cnt++){//LOOP FOR MAX NUMBER OF BUCKET THAT I CAN CHECK THE M POINTS

        int bucket_size = hyper_hashTable[label].size();
        
        if(M_size <= bucket_size){//IF THE POINTS INSITE THE BUCKET ARE MORE THAN M
            int M_cnt = 1;

            vector<data> points;

            for(auto v: hyper_hashTable[label])
                points.push_back(*v.p);
            
            shuffle(points.begin(), points.end(), g);
            for(auto v: hyper_hashTable[label]){//TAKE M RANDOM POINTS FROM THIS BUCKET
                if(M_cnt <= M_size){
                    v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);
                    M_max_points.push_back(*v.p);
                }
                else    
                    break;
                M_cnt++;
            }
            break;
        }

        if(M_size >= bucket_size){//IF THE POINTS INSITE THE BUCKET ARE LESS THAN M

            for(auto v: hyper_hashTable[label]){//TAKE ALL THE POINTS FROM THIS BUCKET
                v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);
                M_max_points.push_back(*v.p);
            }

            M_size = M_size - bucket_size;

            for(int i = current_bucket; i < tableSize; i++){//FIND THE LABEL FOR THE NEXT BUCKET WITH HAMMING DISTANCE
                
                if(hamming_distance(q_label, i) == current_ham_dist){
                    label = i;              //NEXT BUCKET LABEL
                    current_bucket = i + 1;
                
                    if(current_bucket > tableSize){//IF I HAVE CHECKED ALL THE BUCKETS WITH THE CURRENT HAMMING DISTANCE ITS TIME TO INCREASE THE HAMMING DISTANCE
                        current_bucket = 0;
                        current_ham_dist++;
                    }

                    break;
                
                }
            }
        }

    }

    sort(M_max_points.begin(), M_max_points.end(), comp);//SORT THE FINAL M POINTS BY DISTANCE

    int size = M_max_points.size();
    for(int j = N; j < size; j++)// RETURN ARRAY WITH N POINTS
        M_max_points.pop_back();
    
    return M_max_points;


}

vector<data*> HyperCube::range_search(vector<double> q, int d, int minR, int maxR, int M, int probes){//RANGE SEARCH
   
    vector<data*> M_max_points;//ARRAY FOR M POINTS (MAX NUMBER OF POINTS)
    data temp;

    temp.p_data = q;
    temp.ID = "query";
    temp.d = d;

    bucketList node_temp;//QUERY POINT NODE
    node_temp.p = &temp;

    node_temp.bit_label = label(node_temp);//LABEL FOR CURRENT QUERY NODE
    int label = stoi(node_temp.bit_label, 0, 2);
    int q_label = label;

    random_device rd;
    mt19937 g(rd());

    int M_size = M;

    int current_ham_dist = 1;
    int current_bucket = 0;
    
    for(int probe_cnt = 1; probe_cnt <= probes; probe_cnt++){//LOOP FOR MAX NUMBER OF BUCKET THAT I CAN CHECK THE M POINTS

        int bucket_size = hyper_hashTable[label].size();
        
        if(M_size <= bucket_size){//IF THE POINTS INSITE THE BUCKET ARE MORE THAN M
            int M_cnt = 1;
            
            vector<data> points;
            for(auto v: hyper_hashTable[label])
                points.push_back(*v.p);
            
            //shuffle(points.begin(), points.end(), g);
            for(auto v: hyper_hashTable[label]){//TAKE M RANDOM POINTS FROM THIS BUCKET
                if(M_cnt <= M_size){
                    v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);
                    M_max_points.push_back(v.p);
                }
                else    
                    break;
                M_cnt++;
            }
            break;
        }

        if(M_size >= bucket_size){//IF THE POINTS INSITE THE BUCKET ARE LESS THAN M

            for(auto v: hyper_hashTable[label]){//TAKE ALL THE POINTS FROM THIS BUCKET
                v.p->current_distance = distance(v.p->p_data, temp.p_data, temp.d);
                M_max_points.push_back(v.p);
            }

            M_size = M_size - bucket_size;

            for(int i = current_bucket; i < tableSize; i++){//FIND THE LABEL FOR THE NEXT BUCKET WITH HAMMING DISTANCE
                
                if(hamming_distance(q_label, i) == current_ham_dist){
                    label = i;              //NEXT BUCKET LABEL
                    current_bucket = i + 1;
                
                    if(current_bucket > tableSize){//IF I HAVE CHECKED ALL THE BUCKETS WITH THE CURRENT HAMMING DISTANCE ITS TIME TO INCREASE THE HAMMING DISTANCE
                        current_bucket = 0;
                        current_ham_dist++;
                    }

                    break;
                
                }
            }
        }

    }

   vector<data*> inR_points;

    int size = M_max_points.size();
    for(auto v: M_max_points){
        if( v->current_distance > minR && v->current_distance <= maxR)//MINR < CURRENT DISTANCE <= MAXR
            inR_points.push_back(v);
    }

    return inR_points;//RETURN THE POINTS WIT DISTANCE BETWEEN MINR AND MAXR

}

void HyperCube::print( ){
    for(int i=0; i< tableSize; i++){
        cout<<"///////////////////////////////////////////////////////"<<endl;
        for(int j=0; j< hyper_hashTable[i].size(); j++)
            cout<<hyper_hashTable[i][j].p->ID<<" ";
        cout<<endl;
    }
    cout<<tableSize<<endl;
}