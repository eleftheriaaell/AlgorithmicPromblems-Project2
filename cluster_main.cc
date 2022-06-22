#include "clustering.hpp"
#include <fstream>
#include <sstream>

int main(int argc, char *argv[]){
    
    ifstream infile(argv[2]);
    ifstream configfile(argv[4]);
    ofstream outfile(argv[6]);

    int K;
    int k_LSH;
    int k_HYPER;
    int L;
    int M;
    int probes; 
    double delta = 2.0;
   
    string update = argv[8];
    string method = argv[10];
    
    bool option_complete;
    bool option_silh;
    if(argc == 15){             //with optional included in command
        option_complete = 1;
        option_silh = 1;
    }
    else if(argc == 13){
        string option = argv[11];
        if(option == "-complete"){
            option_complete = 1;
            option_silh = 0;
        }
        else if(option == "-silhouette"){
            option_silh = 1;
            option_complete = 0;
        }
    }
    else{
        option_complete = 0;
        option_silh = 0;
    }
 
    int w = 150;
    int R = 10000;
    
    ifstream countfile(argv[2]);
    int file_lines = count(istreambuf_iterator<char>(countfile), istreambuf_iterator<char>(), '\n');
    int tableSize = file_lines/16;
 
    struct timeval start, end; 

    /////////////////////////// LSH and HYPER //////////////////////////  
    
    string line2;
    string val;
    int value;

    int i = 0;
    while(getline(configfile, line2)){//GET VALUES FROM CONFIGURATION FILE 
        stringstream ss(line2);
        
        ss >> val;
        ss >> value;

        if(i == 0)
            K = value;
        else if(i == 1)
            L = value;
        else if(i == 2)
            k_LSH = value;
        else if(i == 3)
            M = value;
        else if(i == 4)
            k_HYPER = value;
        else
            probes = value;

        i++;
        
    }
   
    LSH* objLSH = new LSH[L];//L OBJECTS FOR LSH

    for(int i = 0; i < L; i++)
        objLSH[i] = LSH(k_LSH, w, tableSize, delta);

    HyperCube objHYPER = HyperCube(k_HYPER, w);//HYPERCUBE OBJECT

    srand(time(NULL));

    ////////////////////// CLUSTERING ///////////////////////////////

    string line;
    string id;
    double element;
    vector<data*> dataList;

    Clustering obj = Clustering(K);  //CLUSTERING OBJECT

    ////////////////////// filling the tables ///////////////////////
    int flag = 0;
    while(getline(infile, line)){//FOR EACH CURVE
        data* p = new data();
        stringstream ss(line);
         
        ss >> id;//CURVE ID
        p->ID = id;
        while( ss >> element )//SPLIT THE LINE AND GET THE ELEMENTS OF VECTORS
            p->p_data.push_back(element);
           
        p->d = p->p_data.size();
            
        dataList.push_back(p);//PUSH THE CURVE DATA IN DATALIST

        obj.init_Push(p);//INITIALIZE CLUSTER

        if(flag == 0){
            if(method == "LSH_Frechet")//V VECTOR FOR LSH FRECHET (2*d)
                for(int i = 0; i < L; i++)
                    objLSH[i].v_calculator(2 * p->d);
            else
                for(int i = 0; i < L; i++)//V VECTOR FOR LSH
                    objLSH[i].v_calculator(p->d);
            
            objHYPER.v_calculator(p->d);//V VECTOR FOR HYPER
            flag = 1;
        }

        //PUSH THE CURVE IN HASTABLE FOR LSH
        if(method == "LSH"){//FOR LSH
            for(int i = 0; i < L; i++)
                objLSH[i].hashPush(dataList[dataList.size() - 1]);
        }

        if(method == "LSH_Frechet")//FOR FRECHET
            for(int i = 0; i < L; i++)
                objLSH[i].DiscreteFrechet_hashPush(dataList[dataList.size() - 1]);

        if(method == "Hypercube")//FOR HYPERCUBE
            objHYPER.hashPush(dataList[dataList.size() - 1]);

    }
    
    if(method == "LSH_Frechet" || update == "Mean_Frechet") //KMEANS FOR LSH FRECHET 
        obj.Frechet_k_means(dataList.front()->d, dataList.size(), dataList);
    else//CLASSIC KMEANS
        obj.k_means(dataList.front()->d, dataList.size(), dataList);//KMEANS
    
    if (method == "Classic"){//CLASSIC MATHOD(LLOYDS)
        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);

        if(update == "Mean_Frechet")//WITH FRECHET DISTANCE
            obj.Frechet_Lloyds(dataList.front()->d, dataList);
        else//CLASSIC
            obj.Lloyds(dataList.front()->d, dataList);

        gettimeofday(&end, NULL);

        outfile << "Algorithm: Lloyds" << endl << endl;
    }
    else if(method == "LSH"){//LSH METHOD
        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);
        
        obj.reverse_range_search_LSH(objLSH, L, R, dataList);
        outfile << "Algorithm: Range Search LSH" << endl << endl;

        gettimeofday(&end, NULL);
    }
    else if(method == "Hypercube"){//HYPERCUBE METHOD
        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);

        obj.reverse_range_search_HYPER(objHYPER, R, dataList, M, probes);
        outfile << "Algorithm: Range Search Hypercube" << endl << endl;

        gettimeofday(&end, NULL);
    }
    else if(method == "LSH_Frechet"){//LSH FRECHET
        gettimeofday(&start, NULL);
        ios_base::sync_with_stdio(false);

        obj.reverse_range_search_LSH_Frechet(objLSH, L, R, dataList);
        outfile << "Algorithm: Range Search Discrete Frechet" << endl << endl;

        gettimeofday(&end, NULL);
    }

    obj.output(outfile, dataList.front()->d);

    double exec_time;
  
    exec_time = (end.tv_sec - start.tv_sec) * 1e6;
    exec_time = (exec_time + (end.tv_usec - start.tv_usec)) * 1e-6;

    outfile << "clustering_time: " << exec_time << endl << endl;
    
    if(option_silh == 1)//OPTIONAL OUTPUT 
        if(update == "Mean_Frechet")//FRECHET SILHOUTTES
            obj.curve_silhouttes(outfile);  
        else//CLASSIC
            obj.silhouttes(outfile);  

    if(option_complete == 1)//OPTIONAL OUTPUT 
        obj.optional(outfile, dataList.front()->d);
        
}