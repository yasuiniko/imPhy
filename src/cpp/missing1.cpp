#include <iostream>
#include <fstream>
#include "gurobi_c++.h"
#include <vector>

using namespace std;

// UPDATED at 5.20pm (Central Time) 04/07/2017
// Adding two stage optimization part

/*
 This file is part of imPhy, a pipeline for evaluating the quality of
 phylogenetic imputation software.
 Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis
 
 imPhy uses GTP, which is Copyright (C) 2008, 2009  Megan Owen, Scott Provan
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
    int nmax=-1;
    
    inp >> nSpecies >> nLeaves; // read number of species and number of leafs (global)
    
    int *nIndividuals=new int[nSpecies]; // create array to hold the number of individuals per species
    
    for(int i=0;i<nSpecies;i++){
        inp >> nIndividuals[i]; // read number of individuals per species
        if(nIndividuals[i]>nmax){
            nmax=nIndividuals[i];
        }
    }
    
    cout << nmax << endl;
    
    vector<TreeDistance>TreeDistances;
    
    float **avg=new float*[nLeaves];
    int **num=new int*[nLeaves];
    
    for(int i=0;i<nLeaves;i++){
        avg[i]=new float[nLeaves];
        num[i]=new int[nLeaves];
        for(int j=i+1;j<nLeaves;j++){
            avg[i][j]=0;
            num[i][j]=0;
        }
    }
    
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
                if(temp>=0){
                    avg[i][j]+=temp;
                    num[i][j]++;
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
    
    float overallAvg=0;
    int overallCount=0;
    
    for(int i=0;i<nLeaves;i++){
        for(int j=i+1;j<nLeaves;j++){
            avg[i][j]/=num[i][j];
            //cout << i << "\t" << j << ": " << avg[i][j] << "\t" << num[i][j] << endl;
            overallAvg+=num[i][j]*avg[i][j];
            overallCount+=num[i][j];
        }
    }
    
    overallAvg/=overallCount;
    
    //cout << endl << endl << overallAvg/overallCount << endl;
    
    //for(int i=0;i<nSpecies;i++){
    //    cout << "Species " << i << ": " << eps[i] << endl;
    //}
    
    cout << endl << m << endl;
    
    int **yVal=new int*[nLeaves];
    
    for(int i=0;i<nLeaves;i++){
        yVal[i]=new int[nSpecies];
        for(int k=0;k<nSpecies;k++){
            yVal[i][k]=0;
        }
    }
    
    int ***wVal=new int**[nLeaves];
    
    for(int i=0;i<nLeaves;i++){
        wVal[i]=new int*[nLeaves];
        for(int j=0;j<nLeaves;j++){
            wVal[i][j]=new int[nSpecies];
            for(int k=0;k<nSpecies;k++){
                wVal[i][j][k]=0;
            }
        }
    }
    
    try{
        GRBEnv env;
        GRBModel model(env);
        
        GRBVar **y=new GRBVar*[nLeaves];
        for(int i=0;i<nLeaves;i++){
            y[i]=new GRBVar[nSpecies];
            for(int k=0;k<nSpecies;k++){
                char name[10];
                sprintf(name, "y(%d,%d)", i, k);
                y[i][k]=model.addVar(0,1,0,'B',name);
            }
        }
        cout << "y variables DONE\n";
        
        GRBVar ***w=new GRBVar**[nSpecies];
        for(int k=0;k<nSpecies;k++){
            w[k]=new GRBVar*[nLeaves];
            for(int i=0;i<nLeaves;i++){
                w[k][i]=new GRBVar[nLeaves];
                for(int j=i+1;j<nLeaves;j++){
                        char name[20];
                        sprintf(name, "w(%d,%d,%d)", k, i, j);
                        w[k][i][j]=model.addVar(0,1,avg[i][j],'B',name);
                }
            }
        }
        
        cout << "w variables DONE\n";
        
        
        
        model.update();
        
        for(int k=0;k<nSpecies;k++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    model.addConstr(w[k][i][j]<=y[i][k]);
                    model.addConstr(w[k][i][j]<=y[j][k]);
                    model.addConstr(w[k][i][j]>=y[i][k]+y[j][k]-1);
                }
            }
        }
        
        cout << "Linearization DONE\n";
        
        GRBLinExpr assign;
        
        for(int i=0;i<nLeaves;i++){
            assign.clear();
            for(int k=0;k<nSpecies;k++){
                assign+=y[i][k];
            }
            model.addConstr(assign==1);
        }
        
        
        for(int k=0;k<nSpecies;k++){
            assign.clear();
            for(int i=0;i<nLeaves;i++){
                assign+=y[i][k];
            }
            model.addConstr(assign>=1);
            model.addConstr(assign<=nmax);
        }
        
        cout << "Assignment DONE\n";
        
        model.optimize();
        
        for(int k=0;k<nSpecies;k++){
            for(int i=0;i<nLeaves;i++){
                yVal[i][k]=y[i][k].get(GRB_DoubleAttr_X);
            }
        }
        
        for(int k=0;k<nSpecies;k++){
            cout << "Species " << k << endl;
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    wVal[i][j][k]=w[k][i][j].get(GRB_DoubleAttr_X);
                    if(wVal[i][j][k]>0){
                        cout << i << ", " << j << endl;
                    }
                }
            }
        }
        
        
    }catch(GRBException e){ // catching exception of solver and printing error code and error message
        cout << e.getErrorCode() << " " << e.getMessage() << endl;
        }catch(...){
            cout << "Error during optimization.";
        }
    
    try{
        GRBEnv env;
        GRBModel model(env);
        
        model.getEnv().set("BarHomogeneous", "1");
        model.getEnv().set("NumericFocus", "3");
        
        
        GRBVar **xbar=new GRBVar*[nLeaves];
        for(int i=0;i<nLeaves;i++){
            xbar[i]=new GRBVar[nLeaves];
            for(int j=i+1;j<nLeaves;j++){
                char name[20];
                sprintf(name, "xbar(%d,%d)", i, j);
                xbar[i][j]=model.addVar(0,GRB_INFINITY,0,'C',name);
            }
        }
        
        
        GRBVar ***x=new GRBVar**[nGenes];
        for(int t=0;t<nGenes;t++){
            x[t]=new GRBVar*[nLeaves];
            for(int i=0;i<nLeaves;i++){
                x[t][i]=new GRBVar[nLeaves];
                for(int j=0;j<nLeaves;j++){
                    if(i!=j){
                    if (TreeDistances.at(t).getDistanceAt(i,j)==-1){
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(0,m,0,'C',name);
                    }
                    else{
                        char name[10];
                        sprintf(name, "x(%d,%d,%d)", t, i, j);
                        x[t][i][j]=model.addVar(TreeDistances.at(t).getDistanceAt(i,j),TreeDistances.at(t).getDistanceAt(i,j),0,'C',name);
                    }
                    }
                }
            }
        }
        
        cout << "x variables DONE\n";
        
        /*GRBVar **y=new GRBVar*[nLeaves];
        for(int i=0;i<nLeaves;i++){
            y[i]=new GRBVar[nSpecies];
            for(int k=0;k<nSpecies;k++){
                char name[10];
                sprintf(name, "y(%d,%d)", i, k);
                y[i][k]=model.addVar(0,yVal[i][k],0,'B',name);
                /*if(i==0 && k==3){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==1 && k==0){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==2 && k==1){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==3 && k==2){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==4 && k==0){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==5 && k==2){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==6 && k==4){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else if(i==7 && k==4){
                    y[i][k]=model.addVar(1,1,0,'B',name);
                }
                else{
                    y[i][k]=model.addVar(0,1,0,'B',name);
                }*//*
            }
        }
        cout << "y variables DONE\n";
        
        GRBVar ***w=new GRBVar**[nSpecies];
        for(int k=0;k<nSpecies;k++){
            w[k]=new GRBVar*[nLeaves];
            for(int i=0;i<nLeaves;i++){
                w[k][i]=new GRBVar[nLeaves];
                for(int j=i+1;j<nLeaves;j++){
                    if(avg[i][j]<overallAvg){
                        char name[20];
                        sprintf(name, "w(%d,%d,%d)", k, i, j);
                        w[k][i][j]=model.addVar(0,1,0,'B',name);
                    }
                    else{
                        char name[20];
                        sprintf(name, "w(%d,%d,%d)", k, i, j);
                        w[k][i][j]=model.addVar(0,0,0,'B',name);
                    }
                }
            }
        }
        cout << "w variables DONE\n";
        
        */
        
        model.update();
        
        GRBQuadExpr def1;
        
        
        for(int i=0;i<nLeaves;i++){
            for(int j=i+1;j<nLeaves;j++){
                def1.clear();
                for(int l=0;l<nLeaves;l++){
                    if(l!=i && l!=j){
                        for(int t=0;t<nGenes;t++){
                            def1+=((x[t][i][l]-x[t][j][l])*(x[t][i][l]-x[t][j][l]));
                        }
                    }
                }
                char name[20];
                sprintf(name, "QC[%d,%d]", i, j);
                model.addQConstr(xbar[i][j]>=def1,name);
            }
        }
        
        cout << "Quadratic constraint DONE\n";
        
        /*for(int k=0;k<nSpecies;k++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    model.addConstr(w[k][i][j]<=y[i][k]);
                    model.addConstr(w[k][i][j]<=y[j][k]);
                    model.addConstr(w[k][i][j]>=y[i][k]+y[j][k]-1);
                }
            }
        }
        
        cout << "Linearization DONE\n";*/
        
        /*GRBLinExpr assign;
        
        for(int i=0;i<nLeaves;i++){
            assign.clear();
            for(int k=0;k<nSpecies;k++){
                assign+=y[i][k];
            }
            model.addConstr(assign==1);
        }
        
        
        for(int k=0;k<nSpecies;k++){
            assign.clear();
            for(int i=0;i<nLeaves;i++){
                assign+=y[i][k];
            }
            model.addConstr(assign>=1);
            model.addConstr(assign<=nmax);
        }
        
        cout << "Assignment DONE\n";
        */
        
        for(int t=0;t<nGenes;t++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    model.addConstr(x[t][i][j]==x[t][j][i]);
                }
            }
        }
        
        
        
        GRBQuadExpr obj, const1;
        GRBLinExpr const2;
    
        for(int k=0;k<nSpecies;k++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    //cout << wVal[i][j][k] << endl;
                    obj+=wVal[i][j][k]*xbar[i][j];
                }
            }
        }
        
        cout << "Objective DONE\n";
        
        model.setObjective(obj);
        
        
        
        model.update();
        
        // Uncomment below to write the model file.
        // model.write("./model/test.lp");
        
        // Uncomment below to add TimeLimit (in seconds)
        // model.getEnv().set(GRB_DoubleParam_TimeLimit, 200);
        
        model.optimize();
        char fname[20];
        sprintf(fname, "./sol/%s.sol", argv[1]);
        ofstream f(fname);
        for(int t=0;t<nGenes;t++){
            for(int i=0;i<nLeaves;i++){
                for(int j=i+1;j<nLeaves;j++){
                    if(TreeDistances.at(t).getDistanceAt(i,j)==-1)
                        f << x[t][i][j].get(GRB_DoubleAttr_X) << "\t";
                    else
                        f << TreeDistances.at(t).getDistanceAt(i,j) << "\t";
                }
            }
            f << endl;
        }
        
        
        
        
    }catch(GRBException e){ // catching exception of solver and printing error code and error message
        cout << e.getErrorCode() << " " << e.getMessage() << endl;
    }catch(...){
        cout << "Error during optimization.";
    }
    
    
    
    return 0;
}
