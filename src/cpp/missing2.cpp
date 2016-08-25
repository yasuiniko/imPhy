#include <iostream>
#include <fstream>
#include "gurobi_c++.h"
#include <vector>

using namespace std;

/*
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright (C) 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
// UPDATED at 2.25pm (CST) 07/19/2016

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
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(0,1,0,'C',name);
                    }
                    else{
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(TreeDistances.at(t).getDistanceAt(i,j)/m,TreeDistances.at(t).getDistanceAt(i,j)/m,0,'C',name);
                    }
                }
            }
        }
        
        model.update();
        
        GRBQuadExpr obj;
        
        for(int t1=0;t1<nGenes;t1++){
            for(int t2=t1+1;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                        for(int k=j+1;k<nLeaves;k++){
                            obj+=((x[t1][i][j]-x[t2][i][j])-(x[t1][i][k]-x[t2][i][k]))*((x[t1][i][j]-x[t2][i][j])-(x[t1][i][k]-x[t2][i][k]));
                        }
                    }
                }
            }
        }
        
        model.setObjective(obj);
        
        for(int t1=0;t1<nGenes;t1++){
            for(int t2=t1+1;t2<nGenes;t2++){
                for(int i=0;i<nLeaves;i++){
                    for(int j=i+1;j<nLeaves;j++){
                        if(translatePosToSpecies(nSpecies, nIndividuals, i)==translatePosToSpecies(nSpecies, nIndividuals, j)){
                            model.addConstr(x[t1][i][j]-x[t2][i][j]<=eps[translatePosToSpecies(nSpecies, nIndividuals, i)]/m);
                            model.addConstr(x[t1][i][j]-x[t2][i][j]>=-eps[translatePosToSpecies(nSpecies, nIndividuals, i)]/m);
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
                    f << m*x[t][i][j].get(GRB_DoubleAttr_X) << "\t";
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
                        cout << t << " " << i << " " << j << "\t: " << m*x[t][i][j].get(GRB_DoubleAttr_X) << " vs. " << trueDistance << endl;
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
