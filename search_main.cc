#include "LSH.hpp"
#include "hypercube.hpp"
#include <fstream>
#include <sstream>
#include <time.h>

int main(int argc, char *argv[]){
    
    ifstream infile(argv[2]);
    ofstream outfile;
    ifstream infileq(argv[4]);
    
    int k;
    int L;
    int N;
    int M;
    int probes;
    string algorithm;
    string metric;
    double delta;

    if(argc == 21){
        outfile.open(argv[14]);        //if taken in right position with other given argvs
        k = atoi(argv[6]);
        L = atoi(argv[8]);
        M = atoi(argv[10]);
        probes = atoi(argv[12]);
        algorithm = argv[16];
        metric = argv[18];
        delta = atoi(argv[20]);
    }
    else if(argc == 19){
        outfile.open(argv[14]);        //if taken in right position with other given argvs
        k = atoi(argv[6]);
        L = atoi(argv[8]);
        M = atoi(argv[10]);
        probes = atoi(argv[12]);
        algorithm = argv[16];
        delta = atoi(argv[18]); 
    }
    else{
        outfile.open(argv[6]);      //if taken without other argvs is placed in 6th position
        algorithm = argv[8];
        
        if(algorithm == "Frechet"){//ARGS FOR FRECHET
            metric = argv[10];
            delta = atoi(argv[12]);

            k = 4;
            L = 1;
        }
        else
            delta = atoi(argv[10]);

        if(algorithm == "LSH"){  //ARGS FOR LSH 
            k = 4;
            L = 5;
        }
        else if(algorithm == "Hypercube"){//ARGS FOR HYPERCUBE
            k = 4;
            M = 10;
            probes = 2;
        } 
    }

    N = 1;
    int w = 100;
    
    ifstream countfile(argv[2]);
    int file_lines = count(istreambuf_iterator<char>(countfile), istreambuf_iterator<char>(), '\n');
    int tableSize = file_lines/16;
    
    struct timeval start, end;

    string line;
    string id;
    double element;
    vector<data*> dataList;//LIST WITH ALL THE CURVES DATA

    LSH* obj_LSH;
    HyperCube obj_Hyper;

    if(algorithm == "LSH" || algorithm == "Frechet"){
        obj_LSH = new LSH[L];

        for(int i = 0; i < L; i++)//LSH CONSTRUCTORS
            obj_LSH[i] = LSH(k, w, tableSize, delta); 
    }
    else{
        obj_Hyper = HyperCube(k, w);   
        srand(time(NULL));  
    }

    ////////////////////// filling the hashtables ///////////////////////
    int flag = 0;
    while(getline(infile, line)){//FOR EACH INPUT CURVE
        
        vector<double> current_data;
        data* p = new data();
        stringstream ss(line);

        ss >> id;
        p->ID = id;//CURVE ID

        if(algorithm == "LSH" || algorithm == "Hypercube" || metric == "Discrete"){
            while( ss >> element )//SPLIT THE LINE AND GET THE ELEMENTS OF CURVE
                p->p_data.push_back(element);
            
            p->d = p->p_data.size();
        }
        else if(metric == "Continuous"){
            while( ss >> element )//SPLIT THE LINE AND GET THE ELEMENTS OF CURVE
                current_data.push_back(element);
            
            p->d = current_data.size();
            
            p->true_curve = current_data;
            filtering(current_data, p);//FILTERING FOR CONTINUOUS
        }            

        dataList.push_back(p);//PUSH CURVE IN DATALIST

        if(algorithm == "LSH"){
            if(flag == 0){//CALCULATE ALL THE V FOR EACH H FUNCTION ONLY ONCE
                for(int i = 0; i < L; i++)
                    obj_LSH[i].v_calculator(p->d);
                flag = 1;
            }
        }
        else if(algorithm == "Hypercube"){

            if(flag == 0){//CALCULATE ALL THE V FOR EACH H FUNCTION ONLY ONCE
                obj_Hyper.v_calculator(p->d);
                flag = 1;
            }
        }
        else{
            if(flag == 0){//CALCULATE ALL THE V FOR EACH H FUNCTION ONLY ONCE FOR
                for(int i = 0; i < L; i++)
                    obj_LSH[i].v_calculator(2 * p->d);
                flag = 1;
            }
        }
       
        if(algorithm == "LSH")
            for(int i = 0; i < L; i++)//FOR EACH HASH TABLE PUSH THE CURVE
                obj_LSH[i].hashPush(dataList[dataList.size() - 1]);
        if(metric == "Discrete")
            for(int i = 0; i < L; i++)//FOR EACH HASH TABLE PUSH THE P
                obj_LSH[i].DiscreteFrechet_hashPush(dataList[dataList.size() - 1]);
        else if(metric == "Continuous")
            for(int i = 0; i < L; i++)//FOR EACH HASH TABLE PUSH THE CURVE
                obj_LSH[i].ContinuousFrechet_hashPush(dataList[dataList.size() - 1]);
        if(algorithm == "Hypercube")
                obj_Hyper.hashPush(dataList[dataList.size() - 1]);//PUSH CURVE IN HASH TABLE
    
    }

    //////////////////////////////////////////////////////////////////

    string lineq;
    string idq;
    double elementq;
    
    double tApproximateAverage = 0.0;
    double tTrueAverage = 0.0;

    double MAF = 0.0;

    ifstream countfile_q(argv[4]);
    int file_lines_q = count(istreambuf_iterator<char>(countfile_q), istreambuf_iterator<char>(), '\n');
        int c = 0;
    while(getline(infileq, lineq)){//FOR EACH QUERY CURVE
        vector<double> q_data;
        data* q = new data();
        stringstream ss(lineq);
        
        ss >> idq;
        string q_id = idq;//QUERY ID
        while( ss >> elementq )//SPLIT THE LINE AND GET ALL THE ELEMENTS
            q_data.push_back(elementq);

        int d = q_data.size();

        if(metric == "Continuous")
            filtering(q_data, q);
  
        //////////////////////////// finding the kNN /////////////////////
        
        vector<vector<data>> knn_arr;
        data kNNarr[N*L];
        
        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);

        if(algorithm == "LSH")
            for(int i = 0; i < L; i++)// FROM EACH HASH TABLE TAKE N NEAREST NEIGHBORS
                knn_arr.push_back(obj_LSH[i].kNNsearch(q_data, d, N));
        else if(metric == "Discrete")
            for(int i = 0; i < L; i++)// FROM EACH HASH TABLE TAKE N NEAREST NEIGHBORS
                knn_arr.push_back(obj_LSH[i].DiscreteFrechet_kNNsearch(q_data, d, N));
        else if(metric == "Continuous")
            for(int i = 0; i < L; i++)// FROM EACH HASH TABLE TAKE N NEAREST NEIGHBORS
                knn_arr.push_back(obj_LSH[i].ContinuousFrechet_kNNsearch(q->p_data, d, N));
        else if(algorithm == "Hypercube"){
            knn_arr.push_back(obj_Hyper.kNNsearch(q_data, d, N, M, probes));//GET N NEAREST NEIGHBORS
        }
        
        gettimeofday(&end, NULL);

        double exec_time;
  
        exec_time = (end.tv_sec - start.tv_sec) * 1e6;
        exec_time = (exec_time + (end.tv_usec - start.tv_usec)) * 1e-6;
       
        int knn_size = 0;
        for(auto k: knn_arr){
            for(auto l: k){
                kNNarr[knn_size] = l;
                knn_size++;
            }
        }       

        if(knn_size == 0)
            continue;
        
        sort(&kNNarr[0], &kNNarr[knn_size - 1], compare_distance);

        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);
        
        if(algorithm == "LSH" || algorithm == "Hypercube")
            for(int i = 0; i < dataList.size(); i++)//FOR EACH CURVE GET TRUE DISTANCE
                dataList[i]->true_distance = distance(dataList[i]->p_data, q_data, dataList[i]->d);
        else if(metric == "Discrete")
            for(int i = 0; i < dataList.size(); i++)//FOR EACH CURVE GET TRUE DISTANCE
                dataList[i]->true_distance = Discrete_Frechet_Distance(dataList[i]->p_data, q_data);
        else if(metric == "Continuous")
            for(int i = 0; i < dataList.size(); i++)//FOR EACH CURVE GET TRUE DISTANCE
                dataList[i]->true_distance = ContinuousFrechet_distance(dataList[i]->p_data, q->p_data);
       
        gettimeofday(&end, NULL);

        double exec_time_2;
  
        exec_time_2 = (end.tv_sec - start.tv_sec) * 1e6;
        exec_time_2 = (exec_time_2 + (end.tv_usec - start.tv_usec)) * 1e-6;
        
        sort(&dataList[0], &dataList[dataList.size() - 1], compare_dist);
       
        string prev = kNNarr[0].ID;
       
        outfile << "Query: " << q_id << endl;
        if(metric == "Continuous")
            outfile << "Algorithm: " << algorithm << " with Metric: " << metric << endl;
        else if(metric == "Discrete")
            outfile << "Algorithm: " << algorithm << " with Metric: " << metric << endl;
        else
            outfile << "Algorithm: " << algorithm << endl;
        outfile << "Approximate Nearest neighbour: " << kNNarr[0].ID << endl;
        outfile << "True Nearest neighbour: " << dataList[0]->ID << endl;
        outfile << "distanceApproximate: " << kNNarr[0].current_distance << endl;
        outfile << "distanceTrue: " << dataList[0]->true_distance << endl << endl;
        c++;
        
        double currentMAF = kNNarr[0].current_distance/dataList[0]->true_distance;
        if(currentMAF > MAF)
            MAF = currentMAF;
        
        tApproximateAverage += exec_time;
        tTrueAverage += exec_time_2;
   
    }

    tApproximateAverage /= file_lines_q;
    tTrueAverage /= file_lines_q;

    outfile << "tApproximateAverage: " << tApproximateAverage << endl;
    outfile << "tTrueAverage: " << tTrueAverage << endl;
    outfile << "MAF: " << MAF << endl;

    outfile.close(); 
    
    return 0;
}