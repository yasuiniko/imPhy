#include <iostream>
#include <fstream>
#include "gurobi_c++.h"
#include <vector>

using namespace std;

// COMMENTED TO SIGNAL LATEST VERSION -- 07/12/2016 10:46pm CST, 07/13/2016 12:46pm (?) JST

int translatePosToSpecies(int nSp, int *nInds, int i){
    int *nStart=new int[nSp];
    int *nEnd=new int[nSp];
    nStart[0]=0;
    nEnd[0]=nInds[0]-1;
    for(int j=1;j<nSp;j++){
        nStart[j]=nEnd[j-1]+1;
        nEnd[j]=nStart[j]+nInds[j]-1;
    }
    
    for(int j=0;j<nSp;j++){
        if (nStart[j]<=i && i<=nEnd[j]){
            return j;
        }
    }
    return -1;
}


int main(int argc, const char * argv[]) {
    char inFile[50];
    sprintf(inFile, "./data/%s.txt", argv[1]);
    ifstream inp(inFile); // open the argument file
    
    int nSpecies, nLeaves, nGenes=0;
    
    inp >> nSpecies >> nLeaves; // read number of species and number of leafs (global)
    
    int *nIndividuals=new int[nSpecies]; // create array to hold the number of individuals per species
    
    for(int i=0;i<nSpecies;i++){
        inp >> nIndividuals[i]; // read number of individuals per species
    }
    
    
    vector<float**>TreeDistances;
    
    float m=-1;
    float *eps=new float[nSpecies];
    for(int s=0;s<nSpecies;s++){
        eps[s]=0;
    }
    
    bool done=false;
    do{
        float temp;
        float **TreeD=new float*[nLeaves];
        for(int i=0;i<nLeaves;i++){
            TreeD[i]=new float[nLeaves];
        }
        for(int i=0;i<nLeaves;i++){
            for(int j=i+1;j<nLeaves;j++){
                inp >> temp;
                if(temp>m){
                    m=temp; // Just to find maximum distance in all trees
                }
                cout << i << " " << j << " " << translatePosToSpecies(nSpecies,nIndividuals,i) << " " << translatePosToSpecies(nSpecies,nIndividuals,j) << " " << temp << endl;
                if(translatePosToSpecies(nSpecies, nIndividuals, i)==translatePosToSpecies(nSpecies, nIndividuals, j)){
                    if(temp>eps[translatePosToSpecies(nSpecies, nIndividuals, i)]){
                        cout << i << " " << j << " " << temp << endl;
                        eps[translatePosToSpecies(nSpecies, nIndividuals, i)]=temp;
                    }
                }
                TreeD[i][j]=temp;
            }
        }
        if(inp.eof()){
            done=true;
        }
        else{
            nGenes++;
            TreeDistances.push_back(TreeD);
        }
    }while(!done);
    
    
    
    
    try{
        GRBEnv env;
        GRBModel model(env);
        
        model.getEnv().set("BarHomogeneous", "1");
        
        GRBVar ***x=new GRBVar**[nGenes];
        for(int t=0;t<nGenes;t++){
            x[t]=new GRBVar*[nLeaves];
            for(int i=0;i<nLeaves;i++){
                x[t][i]=new GRBVar[nLeaves];
                for(int j=i+1;j<nLeaves;j++){
                    if (TreeDistances.at(t)[i][j]==-1){
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(0,1,0,'C',name);
                    }
                    else{
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(TreeDistances.at(t)[i][j]/m,TreeDistances.at(t)[i][j]/m,0,'C',name);
                    }
                }
            }
        }
        
        model.update();
        
        GRBQuadExpr obj, const1;
        GRBLinExpr const2;
        
        // ALL SPECIES OBJECTIVE
        
        /*for(int t1=0;t1<nGenes;t1++){
            for(int t2=t1+1;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                        if(translatePosToSpecies(nSpecies, nIndividuals, i)!=translatePosToSpecies(nSpecies, nIndividuals, j))
                            obj+=(x[t1][i][j]-x[t2][i][j])*(x[t1][i][j]-x[t2][i][j]);
                    }
                }
            }
        }*/
        
        // ONLY DIFFERENT SPECIES OBJECTIVE
        
        for(int t1=0;t1<nGenes;t1++){
            for(int t2=t1+1;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                        if(translatePosToSpecies(nSpecies, nIndividuals, i)!=translatePosToSpecies(nSpecies, nIndividuals, j))
                            obj+=(x[t1][i][j]-x[t2][i][j])*(x[t1][i][j]-x[t2][i][j]);
                    }
                }
            }
        }
        
        model.setObjective(obj);
        
        for(int t1=0;t1<nGenes;t1++){
            for(int t2=0;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                        if(translatePosToSpecies(nSpecies, nIndividuals, i)==translatePosToSpecies(nSpecies, nIndividuals, j)){
                            float b=eps[translatePosToSpecies(nSpecies, nIndividuals, i)]*eps[translatePosToSpecies(nSpecies, nIndividuals, i)];
                            model.addQConstr((x[t1][i][j]-x[t2][i][j])*(x[t1][i][j]-x[t2][i][j])<=b/(m*m));
                        }
                    }
                }
            }
        }
        
        
        model.update();
        model.write("./model/test.lp");
        model.optimize();
        char fname[20];
        sprintf(fname, "./sol/%s.sol", argv[1]);
        ofstream f(fname);
        for(int t=0;t<nGenes;t++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    f << x[t][i][j].get(GRB_DoubleAttr_X)*m << "\t";
                }
            }
            f << endl;
        }
        
        sprintf(fname, "./data/%s_true.txt", argv[1]);
        
        ifstream in2(fname);
        float trueDistance;
        for(int t=0;t<nGenes;t++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    in2 >> trueDistance;
                    if(TreeDistances.at(t)[i][j]==-1){
                        cout << t << " " << i << " " << j << "\t: " << x[t][i][j].get(GRB_DoubleAttr_X)*m << " vs. " << trueDistance << endl;
                    }
                }
            }
        }
        
        
        
    }catch(GRBException e){ // catching exception of solver and printing error code and error message
        cout << e.getErrorCode() << " " << e.getMessage() << endl;
    }catch(...){
        cout << "Error during optimization.";
    }
    
    
    
    return 0;
}
