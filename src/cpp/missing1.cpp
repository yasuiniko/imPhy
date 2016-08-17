#include <iostream>
#include <fstream>
#include "gurobi_c++.h"
#include <vector>

using namespace std;

// UPDATED at 10.09pm (Italy Time) 08/10/2016

class TreeDistance{
private:
    vector<float> distances;
    int nLeaves;
    int translateLeavesToIndex(int l1, int l2){
        int firstIndex=0;
        for(int i=0;i<l1;i++){
            firstIndex+=nLeaves-1-i;
        }
        return firstIndex+(l2-l1-1);
    }
public:
    TreeDistance(){
        nLeaves=-1;
    }
    TreeDistance(int n){
        nLeaves=n;
    }
    TreeDistance(int n, vector<float> v){
        nLeaves=n; distances=v;
    }
    vector<float> getDistances(){
        return distances;
    }
    float getDistanceAt(int index){
        if(index<nLeaves*(nLeaves-1)/2){
            return distances.at(index);
        }
        else{
            return -1;
        }
    }
    float getDistanceAt(int leaf1, int leaf2){
        if(leaf1==leaf2){
            return -1;
        }
        else if(leaf1>leaf2){
            int temp=leaf1;
            leaf1=leaf2;
            leaf2=temp;
        }
        return getDistanceAt(translateLeavesToIndex(leaf1, leaf2));
    }
    
};

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
    
    
    vector<TreeDistance>TreeDistances;
    
    float m=-1;
    float *eps=new float[nSpecies];
    for(int s=0;s<nSpecies;s++){
        eps[s]=0;
    }
    
    bool done=false;
    do{
        float temp;
        vector<float> TreeD;
        for(int i=0;i<nLeaves;i++){
            for(int j=i+1;j<nLeaves;j++){
                inp >> temp;
                if(temp>m){
                    m=temp; // Just to find maximum distance in all trees
                }
                if(translatePosToSpecies(nSpecies, nIndividuals, i)==translatePosToSpecies(nSpecies, nIndividuals, j)){
                    if(temp>eps[translatePosToSpecies(nSpecies, nIndividuals, i)]){
                        eps[translatePosToSpecies(nSpecies, nIndividuals, i)]=temp;
                    }
                }
                TreeD.push_back(temp);
            }
        }
        if(inp.eof()){
            done=true;
        }
        else{
            nGenes++;
            TreeDistance x(nLeaves,TreeD);
            TreeDistances.push_back(x);
        }
    }while(!done);
    
    for(int i=0;i<nSpecies;i++){
        cout << "Species " << i << ": " << eps[i] << endl;
    }
    
    cout << endl << m << endl;
    
    
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
                    if (TreeDistances.at(t).getDistanceAt(i,j)==-1){
                        if(translatePosToSpecies(nSpecies, nIndividuals, i)!=translatePosToSpecies(nSpecies, nIndividuals, j)){
                            char name[10];
                            sprintf(name, "x(%d,%d,%d)", t, i, j);
                            x[t][i][j]=model.addVar(0,m,0,'C',name);
                        }
                        else{
                            char name[10];
                            sprintf(name, "x(%d,%d,%d)", t, i, j);
                            x[t][i][j]=model.addVar(0,eps[translatePosToSpecies(nSpecies, nIndividuals, i)],0,'C',name);
                        }
                    }
                    else{
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(TreeDistances.at(t).getDistanceAt(i,j),TreeDistances.at(t).getDistanceAt(i,j),0,'C',name);
                    }
                }
            }
        }
        
        model.update();
        
        GRBQuadExpr obj, const1;
        GRBLinExpr const2;
        
        // ALL SPECIES OBJECTIVE
        
        for(int t1=0;t1<nGenes;t1++){
            for(int t2=t1+1;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                            obj+=(x[t1][i][j]-x[t2][i][j])*(x[t1][i][j]-x[t2][i][j]);
                    }
                }
            }
        }
        // ONLY DIFFERENT SPECIES OBJECTIVE
        /*
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
         */
        
        model.setObjective(obj);
        
        /*for(int t1=0;t1<nGenes;t1++){
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
         }*/
        
        
        model.update();
        
        // Uncomment below to write the model file.
        // model.write("./model/test.lp");
        
        model.optimize();
        char fname[20];
        sprintf(fname, "./sol/%s.sol", argv[1]);
        ofstream f(fname);
        for(int t=0;t<nGenes;t++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    f << x[t][i][j].get(GRB_DoubleAttr_X) << "\t";
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
                    if(TreeDistances.at(t).getDistanceAt(i,j)==-1){
                        cout << t << " " << i << " " << j << "\t: " << x[t][i][j].get(GRB_DoubleAttr_X) << " vs. " << trueDistance << endl;
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