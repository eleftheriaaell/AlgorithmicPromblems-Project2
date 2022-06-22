#include "clustering.hpp"

vector<double> Clustering::Mean_Curve(vector<double> x, vector<double> y, int d){ //CALCULATE THE MEAN CURVE BETWEEN TWO CURVES

    int m1 = x.size();
    int m2 = y.size();

    vector<vector<double>> Dynamic_Array(m1, vector<double> (m2));//DYNAMIC ARRAY
    double distance = 0.0;
    
    Dynamic_Array[0][0] = abs(x[0] - y[0]);//C(0,0)
    distance = Dynamic_Array[0][0];
   
    for(int i = 1; i < m1; i++){//C(i,0)

        double cij = abs(x[i] - y[0]);

        if(Dynamic_Array[i - 1][0] > cij){//MAX BETWEEN C(i-1,0) AND C(i,0)
            Dynamic_Array[i][0] = Dynamic_Array[i -1][0];
            continue;
        }
        Dynamic_Array[i][0] = cij;
    }
    
    for(int j = 1; j < m2; j++){//C(0,j)

        double cij = abs(x[0] - y[j]);

        if(Dynamic_Array[0][j - 1] > cij){//MAX BETWEEN C(0,j-1) AND C(0, j)
            Dynamic_Array[0][j] = Dynamic_Array[0][j - 1];
            continue;
        }
        Dynamic_Array[0][j] = cij;

    }

    for(int i = 1; i < m1; i++){//C(i,j)
        for(int j = 1; j < m2; j++){////MAX BETWEEN MIN{C(i-1,j),C(i-1,j-1),C(i,j-1)} AND C(i, j)
            double c1 = Dynamic_Array[i - 1][j];
            double c2 = Dynamic_Array[i - 1][j - 1];
            double c3 = Dynamic_Array[i][j - 1];
            double cij = abs(x[i] - y[j]);
            double c;
            if(c1 < c2){
                if(c1 < c3)
                    c = c1;
                else
                    c = c3;
            }else{ 
                if(c2 < c3)
                    c = c2;
                else
                    c = c3;
            }

            if(c > cij)
                Dynamic_Array[i][j] = c;
            else    
                Dynamic_Array[i][j] = cij;

            if(i == j)
                distance += Dynamic_Array[i][j];

        }
    }

    vector<double> mean_curve;
    int Pi = m1 - 1;
    int Qi = m2 - 1;
    
    mean_curve.insert(mean_curve.begin(), (x[Pi] + y[Qi])/2);
    
    while(Pi > 0 && Qi > 0){ //BACKTRUCKING IN DYNAMIC TABLE AND CREATE THE MEAN CURVE 
        
        double a = Dynamic_Array[Pi - 1][Qi];
        double b = Dynamic_Array[Pi][Qi - 1];
        double c = Dynamic_Array[Pi - 1][Qi -1];
        
        if(a < b){//FIND THE MINIMUM BETWEEN C(i-1,j), C(i,j-1), C(i-1,j-1)
            if(a < c){
                Pi = Pi - 1;
                mean_curve.insert(mean_curve.begin(), (x[Pi] + y[Qi])/2);
            }else{
                Pi = Pi - 1;
                Qi = Qi - 1;
                mean_curve.insert(mean_curve.begin(), (x[Pi] + y[Qi])/2);
            }
        }else{
            if(b < c){
                Qi = Qi - 1;
                mean_curve.insert(mean_curve.begin(), (x[Pi] + y[Qi])/2);
            }else{
                Pi = Pi - 1;
                Qi = Qi - 1;
                mean_curve.insert(mean_curve.begin(), (x[Pi] + y[Qi])/2);
            }
        }
        
       
    }

    mean_curve.insert(mean_curve.begin(), (x[0] + y[0])/2);
    cluster_filtering(mean_curve, d);//FILTERING

    return mean_curve;//RETURN THE dSIZE MEAN CURVE 
}

data* Clustering::Frechet_new_centroid(int i, int d){//CALCULATE NEW CENTROID
   
    vector<vector<double>> current_array;
    vector<vector<double>> new_array;
 
    int cluster_size = cluster_Table[i].size();    
    
    if(cluster_size == 1)//IF THE CLUSTER HAS ONLY THE CENTROID
        return cluster_Table[i].begin()->p;
  
    for(auto v: cluster_Table[i])
        current_array.push_back(v.p->p_data);
        
   
    while(current_array.size() > 1){//LOOP UNTIL WE FOUND THE ROOT OF THE TREE 
        
        int j = 0;
        int odd_flag = 0;
        int current_size = current_array.size();

        if(current_size % 2 != 0){//IF THE NUMBER OF THE NODES IN THE CURRENT LEVEL ITS ODD
            current_size--;
            odd_flag = 1;
        }
        
        while(j < current_size){//CREATE THE NEXT LEVEL WITH THE MEAN CURVES BETWEEN THE CURVES IN THE CURRENT LEVEL
            new_array.push_back(Mean_Curve(current_array[j], current_array[j + 1], d));
            j += 2;
        }

        if(odd_flag == 1)//IF THE ODD FLAG IS 1 ADD THE EXTRA CURVE IN THE NEXT LEVEL
            new_array.push_back(current_array[current_size]);
        
        current_array = new_array;//CURRENT LEVEL = NEXT LEVEL
        
        new_array.clear();
    }
    
    data* new_centroid_data = new data();
    new_centroid_data->ID = "centroid";
    new_centroid_data->d = current_array[0].size();
    new_centroid_data->p_data = current_array[0];
    
    int new_size = new_centroid_data->p_data.size();
    
    return new_centroid_data;//RETURN THE DATA OF THE NEW CENTROID

}

void Clustering::Frechet_Di_first_centroid(int d){//INITIALIZE THE Di(DISTANCE BETWEEN CURVE AND CENTROID) FOR EACH CURVE

    clusterList centroid = cluster_Table[0].front();    

    double max = Discrete_Frechet_Distance(cluster_Table->begin()->p->p_data, centroid.p->p_data); //MAX DISTANCE BETWEEN CENTROID AND CURVE

    for(auto v = cluster_Table->begin(); v != cluster_Table->end(); v++){      
        v->Di = Discrete_Frechet_Distance(v->p->p_data, centroid.p->p_data); 

        if( v->Di > max)
            max = v->Di;
    }

    current_max_dist = max;
    
}

void Clustering::Frechet_k_means(int d, int n, vector<data*> dataList){//K_MEANS++

    Frechet_Di_first_centroid(d);
    
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
        
        for(int i = 0; i <= t; i++){//ERASE ALL THE CURVES
            if(cluster_Table[i].size() < 2)
                continue;
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());
        }
        
        Frechet_table_Push(dataList, t + 2);//ASSIGN ALL THE CURVES TO THE NEAREST CENTROID

        P_array.clear(); 

    }

}

void Clustering::Frechet_Lloyds(int d, vector<data*> dataList){//LLOYDS
    
    int iteration_flag = 0;
    data* new_centro = new data();
    while(iteration_flag < 7){//STOP AFTER 7 ITERATIONS
        
        for(int i = 0; i < K; i++){//FOR EACH CLUSTER
            new_centro = Frechet_new_centroid(i, d);//CALCULATE THE NEW CENTROID                
            
            cluster_Table[i].begin()->p = new_centro;
            
            cluster_Table[i].begin()->Di = -1;
            
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());
        
        }
        iteration_flag++;
        
        Frechet_table_Push(dataList, K);  //ASSING ALL THE CURVES TO THE NEAREST CENTROID    
    }

}

void Clustering::Frechet_table_Push(vector<data*> dataList, int t){//ASSIGN ALL THE CURVES FROM DATALIST TO THE NEAREST CENTROID
    
    for(auto point: dataList){//FOR EACH CURVE
        double min = Discrete_Frechet_Distance(point->p_data, cluster_Table[0].begin()->p->p_data);
        int centroid = 0;
        
        double max = Discrete_Frechet_Distance(point->p_data, cluster_Table[0].begin()->p->p_data);
        for(int i = 1; i < t; i++){//CHECK THE DISTANCE BETWEEN THE CURRENT CURVE AND THE CENTROIDS
            double current = Discrete_Frechet_Distance(point->p_data, cluster_Table[i].begin()->p->p_data);
            if(current < min){//FIND THE NEAREST CENTROID (WITH MIN DISTANCE)
                min = current;
                centroid = i;
            }

            if(current > max)
                max = current;
        }

        if(min != 0){//IF THE CURVE ITS NOT A CENTROID ASSIGN IT TO THE NEAREST CENTROID
            clusterList node;
            node.Di = min;
            node.p = point;
            cluster_Table[centroid].push_back(node);
        }

        current_max_dist = max;
    }

}


void Clustering::Frechet_table_Push_reverse(vector<data*> dataList, int t){//ASSIGN THE CURVES TO CENTROIDS FOR THE REVERSE LSH/HYPERCUBE

    for(auto point: dataList){//FOR EACH CURVE
        if(point->cluster_counter == 1){//IF THE CURVE LIES IN ONE CLUSTER
            point->cluster_counter = 0;
            continue;
        }

        double min = Discrete_Frechet_Distance(point->p_data, cluster_Table[0].begin()->p->p_data);
        int centroid = 0;

        point->cluster_counter = 0;
        for(int i = 1; i < t; i++){//CHECK THE DISTANCE BETWEEN THE CURRENT CURVE AND THE CENTROIDS

            double current = Discrete_Frechet_Distance(point->p_data, cluster_Table[i].begin()->p->p_data);
            if(current < min){//FIND THE NEAREST CENTROID (WITH MIN DISTANCE)
                min = current;
                centroid = i;
            }
        }
        
        if(min != 0){//IF THE CURVE ITS NOT A CENTROID ASSIGN IT TO THE NEAREST CENTROID
            clusterList node;
            node.Di = min;
            node.p = point;
            
            cluster_Table[centroid].push_back(node);
        }

    }

}

void Clustering::reverse_range_search_LSH_Frechet(LSH* objLSH, int L, int R, vector<data*> dataList){//REVERS RANGE SEARCH LSH_FRECHET

    int minR = 0;
    int maxR = R;

    int g_array[K][L];

    int iteration_flag = 0;
    data* new_centro = new data();

    while(iteration_flag < 7){//STOP THE CENTROIDS AFTER 7 ITERATIONS
        
        for(int i = 0; i < K; i++)
            cluster_Table[i].erase(cluster_Table[i].begin().operator++(), cluster_Table[i].end());


        for(int i = 0; i < K; i++)//FIND THE G FOR EACH HASH TABLES
            for(int j = 0; j < L; j++)
                g_array[i][j] = objLSH[j].DiscreteFrechet_gp_cluster(cluster_Table[i].begin()->p);
            
        int radii_flag = 0;

        while(radii_flag < 1){//STOP INCREASING RANGE IF THE ALL CLUSTERS GET NO NEW CURVE
            radii_flag = 0;

            for(int j = 0; j < K; j++){//FOR EACH CENTROID FIND THE CURVES BETWEEN MIN AND MAX RANGE

                vector<vector<data*>> R_arr;
                for(int i = 0; i < L; i++){
                    objLSH[i].current_q_bucket = g_array[j][i];
                    R_arr.push_back(objLSH[i].DiscreteFrechet_range_search(cluster_Table[j].front().p->p_data, cluster_Table[j].front().p->d, minR, maxR));
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

                for(int i = 1; i < Rarr.size(); i++){// ASSAGN THE CURVES BETWEEN MIN AND MAX RANGE IN CURRENT CLUSTER ONLY ONCE
                
                    if(Rarr[i]->ID != prev->ID){
                        clusterList node;
                        node.Di = 0.0;
                        node.p = prev; 
                        node.p->cluster_counter++;//INCREASE THE COUNTER EACH TIME THE CURVE ASSIGN TO A CLUSTER

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

        for(int i = 0; i < K; i++){//REMOVE THE CURVES THAT LIES MORE THAN ONE CLUSTER
            cluster_Table[i].begin()->p->cluster_counter = 0;
            cluster_Table[i].remove_if(cond);
        }
       
        Frechet_table_Push_reverse(dataList, K);//ASSIGN THE REMAIN CURVES TO THE NEAREST CLUSTER
       
        for(int i = 0; i < K; i++){//FOR EACH CLUSTER FIND NEW  CENTROID
           
            new_centro = Frechet_new_centroid(i, dataList[0]->d);
       
            cluster_Table[i].begin()->p = new_centro;
            
            cluster_Table[i].begin()->Di = -1.0;        
        }  
        iteration_flag++;   
    }
    
}

void Clustering::cluster_filtering(vector<double>& mean_curve, int d){//FILTERING FOR CLUSTER

    double e = 1.0;//START WITH e=1.0;
    while(mean_curve.size() > d){//UNTIL THE SIZE OF THE CURVE IS d
        int current_element = 1;
        
        while(current_element < mean_curve.size()){
            double a = mean_curve[current_element - 1];
            double b = mean_curve[current_element];
            double c = mean_curve[current_element + 1];

            if(abs(a-b) <= e && abs(b-c) <= e){//FILTERING CONDITION
                mean_curve.erase(mean_curve.begin() + current_element); 
                if(mean_curve.size() <= d)//IF THE SIZE IS OF THE CURVE IS d
                    break;
        
            }else{
                current_element++;
            
            }
        }
        e = e + 0.1;//INCREASE BY 0.1
    }

    if(mean_curve.size() < d)//PADDING
        for(int j = mean_curve.size(); j < d; j++)
            mean_curve.push_back(33554432.0);
}

void Clustering::curve_silhouttes(ofstream& outfile){//SILHOUTTES

    outfile << "Silhouette: [";
    
    double s[K] = {0};//S(I) ARRAY
    double stotal = 0;
    
    for(int i = 0; i < K; i++){//FOR EACH CLUSTER
        for(auto x: cluster_Table[i]){//FOR EACH CURVE IN CURRENT CLUSTER
            double a = 0, b = 0;
            
            for(auto y: cluster_Table[i])//CALCULATE AVERAGE DISTANCE BETWEEN ALL CURVES IN SAME CLUSTER
                a = a + Discrete_Frechet_Distance(x.p->p_data, y.p->p_data);
            a = a / cluster_Table[i].size();

            double min;
            int centroid;
            if( i != K - 1){//FIND THE NEXT NEAREST CENTROID
                min = Discrete_Frechet_Distance(x.p->p_data, cluster_Table[i + 1].begin()->p->p_data);
                centroid = i + 1;
            }
            else{
                min = Discrete_Frechet_Distance(x.p->p_data, cluster_Table[i - 1].begin()->p->p_data);
                centroid = i - 1;
            }

            for(int j = 0; j < K; j++){
                if(i == j)
                    continue;

                double current = Discrete_Frechet_Distance(x.p->p_data, cluster_Table[j].begin()->p->p_data);
                if(current < min){
                    min = current;
                    centroid = j;
                }
            }

            for(auto y: cluster_Table[centroid])//CALCULATE AVERAGE DISTANCE BETWEEN ALL CURVES IN SECOND BEST CLUSTER
                b = b + Discrete_Frechet_Distance(x.p->p_data, y.p->p_data);
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