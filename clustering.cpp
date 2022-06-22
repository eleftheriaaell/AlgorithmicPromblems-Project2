#include "clustering.hpp"

Clustering::Clustering(){}

Clustering::Clustering(int K_medians){//CONSTRUCTOR
    K = K_medians;
    cluster_Table = new list<clusterList>[K];
}

Clustering::~Clustering(){
}

void Clustering::init_Push(data* p){//PLACE ALL THE POINTS TO THE FIRST CLUSTER
    clusterList node;
    
    node.p = p;
    cluster_Table[0].push_back(node);       
}

void Clustering::Di_first_centroid(int d){//INITIALIZE THE Di(DISTANCE BETWEEN POINT AND CENTROID) FOR EACH POINT

    clusterList centroid = cluster_Table[0].front();    

    double max = distance(cluster_Table->begin()->p->p_data, centroid.p->p_data, d); //MAX DISTANCE BETWEEN CENTROID AND POINT

    for(auto v = cluster_Table->begin(); v != cluster_Table->end(); v++){      
        v->Di = distance(v->p->p_data, centroid.p->p_data, d); 

        if( v->Di > max)
            max = v->Di;
    }

    current_max_dist = max;
    
}

void Clustering::k_means(int d, int n, vector<data*> dataList){//K_MEANS++
    
    Di_first_centroid(d);

    probability node1;
    node1.P = 0.0;

    P_array.push_back(node1);
    
    for(int t = 0; t < K - 1; t++){//ASSIGN ALL CENTROIDS FOR THE FIRST TIME
        double prev_P = 0.0;             
        double curr_P = 0.0; 
        
        probability node1;
        node1.P = 0.0;

        P_array.push_back(node1);

        for(int i = 0; i <= t; i++){// FOR EACH CENTROID 
            
            for(auto v: cluster_Table[i]){//CALCULATE PROBABILITY SUM(D(i))^2
                curr_P = prev_P + pow((v.Di/current_max_dist), 2);
                
                if(v.Di == 0)
                    continue;
                
                probability node;
                node.P = curr_P;
                node.r_point = v.p;

                P_array.push_back(node);//PROBABILITY ARRAY

                prev_P = curr_P;

            }
            
        }                                                          
        
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine gen(seed);

        uniform_real_distribution<> dist(0, P_array[P_array.size() - 1].P);//GET RANDOM NUMBER BETWEEN THE PROBABILITY ARRAY

        double x = dist(gen);
       
        int r = binary_search(P_array, P_array.size(), x);//FIND r POINT (THE NEXT CENTROID)    

        clusterList node;
        node.Di = 0.0;
        node.p = P_array[r].r_point;
        
        cluster_Table[t + 1].push_back(node);

        for(int i = 0; i <= t; i++){//ERASE ALL THE POINTS
            // if(cluster_Table[i].size() < 2)
            //     continue;
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());
        }
        
        table_Push(dataList, t + 2);//ASSIGN ALL THE POINTS TO THE NEAREST CENTROID
        
        P_array.clear();       

    }

}

void Clustering::table_Push(vector<data*> dataList, int t){//ASSIGN ALL THE POINTS FROM DATALIST TO THE NEAREST CENTROID
    
    for(auto point: dataList){//FOR EACH POINT
  
        double min = distance(point->p_data, cluster_Table[0].begin()->p->p_data, point->d);
        int centroid = 0;
      
        double max = distance(point->p_data, cluster_Table[0].begin()->p->p_data, point->d);
        for(int i = 1; i < t; i++){//CHECK THE DISTANCE BETWEEN THE CURRENT POINT AND THE CENTROIDS

            double current = distance(point->p_data, cluster_Table[i].begin()->p->p_data, point->d);

            if(current < min){//FIND THE NEAREST CENTROID (WITH MIN DISTANCE)
                min = current;
                centroid = i;
            }

            if(current > max)
                max = current;
        }

        if(min != 0){//IF THE POINT ITS NOT A CENTROID ASSIGN IT TO THE NEAREST CENTROID
            clusterList node;
            node.Di = min;
            node.p = point;
            cluster_Table[centroid].push_back(node);
        }

        current_max_dist = max;
    }

}

void Clustering::Lloyds(int d, vector<data*> dataList){//LLOYDS

    int iteration_flag = 0;
    data* new_centro = new data();
    while(iteration_flag < K){//STOP THE CENTROIDS FROM UPDATING IF THE DISTANCES BETWEEN THE CURRENT AND NEW CENTROIDS ITS ZERO

        iteration_flag = 0;
        for(int i = 0; i < K; i++){//FOR EACH CLUSTER
            new_centro = new_centroid(i, d);//CALCULATE THE NEW CENTROID
            
            if(distance(cluster_Table[i].begin()->p->p_data, new_centro->p_data, d) == 0)//IF THE DISTUNCE BETWEEN THE CURRENT AND THE NEW CENTROID IS ZERO
                iteration_flag++;
            
            cluster_Table[i].begin()->p = new_centro;
            
            cluster_Table[i].begin()->Di = -1;
            
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());
        
        }
        
        table_Push(dataList, K);  //ASSING ALL THE POINTS TO THE NEAREST CENTROID     
    }

}

data* Clustering::new_centroid(int i, int d){//CALCULATE NEAREST CENTROID
    
    double average_vector[d] = {0.0};
    int cluster_size = cluster_Table[i].size();

    data* new_centroid_data = new data();
    new_centroid_data->ID = "centroid";
    new_centroid_data->d = d;

    if(cluster_size == 1)//IF THE CLUSTER HAS ONLY THE CENNTROID
        return cluster_Table[i].begin()->p;

    for(auto v: cluster_Table[i]){//CALCULATE THE VECTOR OF THE NEW CENTROID WITH AVERAGE VECTOR METHOD
        for(int j = 0; j < d; j++)
            average_vector[j] = average_vector[j] + v.p->p_data[j];
    }       
  
    for(int i = 0; i < d; i++){
        average_vector[i] = average_vector[i] / cluster_size;
        new_centroid_data->p_data.push_back(average_vector[i]);
    }
    
    return new_centroid_data;//RETURN THE DATA OF THE NEW CENTROID

}

int Clustering::binary_search(vector<probability> P_array, int n, double target){//BINARY SEARCH
    
    if (target <= P_array[0].P)//IF THE TARGET SMALLER THAN MIN PROBABILITY
        return 1;
    if (target >= P_array[n - 1].P)//IF THE TARGET IS BIGER THAN MAX PROBABILITY
        return n - 1;

    int i = 0;
    int j = n;
    int mid = 0;

    while (i < j){
        mid = (i + j) / 2;
 
        if (P_array[mid].P == target)//IF TARGET ITS IN THE MIDDLE
            return mid;
 
        if (target < P_array[mid].P){//IF TARGED IS SMALLER THAN THE MIDDLE PROBABILITY
            if (mid > 0 && target > P_array[mid - 1].P)
                return mid;
 
            j = mid;
        }
 
        else{//IF TARGET IS BIGER THAN THE MIDDLE PROBABILITY
            if (mid < n - 1 && target < P_array[mid + 1].P)
                return mid + 1;
            i = mid + 1;
        }
    }

    return mid;
}

void Clustering::table_Push_reverse(vector<data*> dataList, int t){//ASSIGN THE POINTS TO CENTROIDS FOR THE REVERSE LSH/HYPERCUBE
    
    for(auto point: dataList){//FOR EACH POINT
        if(point->cluster_counter == 1){//IF THE POINT LIES IN ONE CLUSTER
            point->cluster_counter = 0;
            continue;
        }

        double min = distance(point->p_data, cluster_Table[0].begin()->p->p_data, point->d);
        int centroid = 0;

        point->cluster_counter = 0;
        for(int i = 1; i < t; i++){//CHECK THE DISTANCE BETWEEN THE CURRENT POINT AND THE CENTROIDS

            double current = distance(point->p_data, cluster_Table[i].begin()->p->p_data, point->d);
            if(current < min){//FIND THE NEAREST CENTROID (WITH MIN DISTANCE)
                min = current;
                centroid = i;
            }
        }
        
        if(min != 0){//IF THE POINT ITS NOT A CENTROID ASSIGN IT TO THE NEAREST CENTROID
            clusterList node;
            node.Di = min;
            node.p = point;
            
            cluster_Table[centroid].push_back(node);
        }

    }

}

void Clustering::reverse_range_search_LSH(LSH* objLSH, int L, int R, vector<data*> dataList){//REVERS RANGE SEARCH LSH
      
    int minR = 0;
    int maxR = R;

    int g_array[K][L];

    int iteration_flag = 0;
    data* new_centro = new data();
    
    while(iteration_flag < K/2){//STOP THE CENTROIDS FROM UPDATING IF THE DISTANCES BETWEEN THE CURRENT AND NEW CENTROIDS ITS ZERO
        
        for(int i = 0; i < K; i++)
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());
        
        iteration_flag = 0;

        for(int i = 0; i < K; i++){//FIND THE G FOR EACH HASH TABLES
            for(int j = 0; j < L; j++){
                g_array[i][j] = objLSH[j].gp_cluster(cluster_Table[i].begin()->p);
            }
        }
        
        int radii_flag = 0;

        while(radii_flag < K/2){//STOP INCREASING RANGE IF THE ALL CLUSTERS GET NO NEW POINT
            radii_flag = 0;
            

            for(int j = 0; j < K; j++){//FOR EACH CENTROID FIND THE POINTS BETWEEN MIN AND MAX RANGE
                
                vector<vector<data*>> R_arr;
                for(int i = 0; i < L; i++){
                    objLSH[i].current_q_bucket = g_array[j][i];
                    R_arr.push_back(objLSH[i].range_search(cluster_Table[j].front().p->p_data, cluster_Table[j].front().p->d, minR, maxR));
                }
                
                vector<data*> Rarr;
                for(auto k: R_arr){
                    for(auto l: k){
                        Rarr.push_back(l);
                    }
                }               
                
                if(Rarr.size() == 0){
                    radii_flag++;
                    continue;
                }
               
                sort(Rarr.begin(), Rarr.end(), compare_ID_pointer);
                
                data* prev = Rarr.front();

                for(int i = 1; i < Rarr.size(); i++){// ASSAGN THE POINTS BETWEEN MIN AND MAX RANGE IN CURRENT CLUSTER ONLY ONCE
                    
                    if(Rarr[i]->ID != prev->ID){
                        clusterList node;
                        node.Di = 0.0;
                        node.p = prev; 
                        node.p->cluster_counter++;//INCREASE THE COUNTER EACH TIME THE POINT ASSIGN TO A CLUSTER

                        cluster_Table[j].push_back(node); 
                        prev = Rarr[i];  
                    
                        continue;           
                    }
                    prev = Rarr[i];  
                
                }
            }

            minR = maxR;
            maxR *= 2;
        }

        for(int i = 0; i < K; i++){//REMOVE THE POINTS THAT LIES MORE THAN ONE CLUSTER
            cluster_Table[i].begin()->p->cluster_counter = 0;
            cluster_Table[i].remove_if(cond);
        }
       
        table_Push_reverse(dataList, K);//ASSIGN THE REMAIN POINTS TO THE NEAREST CLUSTER
       
        for(int i = 0; i < K; i++){//FOR EACH CLUSTER FIND NEW  CENTROID
            new_centro = new_centroid(i, dataList[0]->d);
           
            if(distance(cluster_Table[i].begin()->p->p_data, new_centro->p_data, dataList[0]->d) == 0)//IF THE DISTANCE BETWEEN THE CURRENT AND THE NEW CENTROID IS 0
                iteration_flag++;
                    
            cluster_Table[i].begin()->p = new_centro;
            
            cluster_Table[i].begin()->Di = -1;        
        }       
    }
    
}

bool Clustering::cond(const clusterList& value){//CONDITION FOR THE POINTS THAT LIES IN 2+ CLUSTERS
    return value.p->cluster_counter >= 2;
}

void Clustering::reverse_range_search_HYPER(HyperCube objHYPER, int R, vector<data*> dataList, int M, int probes){//REVERSE RANGE SEARCH WITH HYPERCUBE
        
    int minR = 0;
    int maxR = R;      /* ITS THE SAME WITH REVERSE RANGE SEARCH FOR LSH BUT WITH HYPERCUBE RANGE  SEARCH FUNCTION */

    int g_array[K];

    int iteration_flag = 0;
    data* new_centro = new data();
   
    while(iteration_flag < K){
        for(int i = 0; i < K; i++)
        cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());

        iteration_flag = 0;
                
        for(int i = 0; i < K; i++)
            g_array[i] = stoi(objHYPER.label_cluster(cluster_Table[i].begin()->p)); 

        int radii_flag = 0;
       
        while(radii_flag < K){
            radii_flag = 0;
             
            for(int j = 0; j < K; j++){
                vector<data*> Rarr;
                Rarr = objHYPER.range_search(cluster_Table[j].front().p->p_data, cluster_Table[j].front().p->p_data.size(), minR, maxR, M, probes);           
                
                if(Rarr.size() == 0){
                    radii_flag++;
                    continue;
                }
                
                sort(Rarr.begin(), Rarr.end(), compare_ID_pointer);
                
                data* prev = Rarr.front();
           
                for(int i = 1; i < Rarr.size(); i++){
                   
                    if(Rarr[i]->ID != prev->ID){
                        clusterList node;
                        node.Di = 0.0;
                        node.p = prev; 
                        node.p->cluster_counter++;
                        
                        cluster_Table[j].push_back(node); 
                        prev = Rarr[i];  

                        continue;           
                    }
                    prev = Rarr[i];  
                    
                }

            
            }

            minR = maxR;
            maxR *= 2;
        }
    
        for(int i = 0; i < K; i++){
            cluster_Table[i].begin()->p->cluster_counter = 0;
            cluster_Table[i].remove_if(cond);
        }
        
        table_Push_reverse(dataList, K);
        
        for(int i = 0; i < K; i++){
            new_centro = new_centroid(i, dataList[0]->d);
            
            if(distance(cluster_Table[i].begin()->p->p_data, new_centro->p_data, dataList[0]->d) == 0)
                iteration_flag++;
            
            cluster_Table[i].begin()->p = new_centro;
            
            cluster_Table[i].begin()->Di = -1;        
        }       
    }
    
}

void Clustering::silhouttes(ofstream& outfile){//SILHOUTTES

    outfile << "Silhouette: [";
    
    double s[K] = {0};//S(I) ARRAY
    double stotal = 0;
    
    for(int i = 0; i < K; i++){//FOR EACH CLUSTER
        for(auto x: cluster_Table[i]){//FOR EACH POINT IN CURRENT CLUSTER
            double a = 0, b = 0;
            
            for(auto y: cluster_Table[i])//CALCULATE AVERAGE DISTANCE BETWEEN ALL POINTS IN SAME CLUSTER
                a = a + distance(x.p->p_data, y.p->p_data, x.p->d);
            a = a / cluster_Table[i].size();

            double min;
            int centroid;
            if( i != K - 1){//FIND THE NEXT NEAREST CENTROID
                min = distance(x.p->p_data, cluster_Table[i + 1].begin()->p->p_data , x.p->d);
                centroid = i + 1;
            }
            else{
                min = distance(x.p->p_data, cluster_Table[i - 1].begin()->p->p_data , x.p->d);
                centroid = i - 1;
            }

            for(int j = 0; j < K; j++){
                if(i == j)
                    continue;

                double current = distance(x.p->p_data, cluster_Table[j].begin()->p->p_data, x.p->d);
                if(current < min){
                    min = current;
                    centroid = j;
                }
            }

            for(auto y: cluster_Table[centroid])//CALCULATE AVERAGE DISTANCE BETWEEN ALL POINTS IN SECOND BEST CLUSTER
                b = b + distance(x.p->p_data, y.p->p_data, x.p->d);
            b = b / cluster_Table[centroid].size();

            if(a > b)//SUM(S(i))
                s[i] = s[i] + ((b - a) / a);
            else
                s[i] = s[i] + ((b - a) / b);
        }

        s[i] /= cluster_Table[i].size();//AVERAGE S(i)
        stotal = stotal + s[i];

        outfile << s[i] << ", ";
    }

    outfile << stotal/K << "]" << endl << endl;//AVERAGE S[i]
}

void Clustering::output(ofstream& outfile, int d){

    for(int i = 0; i < K; i++){
        outfile << "Cluster-" << i + 1 << " {size: " << cluster_Table[i].size() << ", centroid: ";
        for(int j = 0; j < d; j++)
            outfile << cluster_Table[i].begin()->p->p_data[j] << " ";
        outfile << endl << endl;
    }

}

void Clustering::optional(ofstream& outfile, int d){
    
    for(int i = 0; i < K; i++){
        outfile << "Cluster-" << i + 1 << " {centroid: ";
        for(int j = 0; j < d; j++)
            outfile << cluster_Table[i].begin()->p->p_data[j] << " ";
        for(auto v: cluster_Table[i])
            outfile << ", " << v.p->ID;
        outfile << "}" << endl << endl;
    }

}