
#include <iostream>
#include <fstream>
#include "gurobi_c++.h"

using namespace std;

/********* Need to provide a file as an argument to the executable in the format:
 //******** nSpecies, nLeafs, nIndividuals[0], nIndividuals[1], ..., nIndividuals[nSpecies-1]
 //******** d_12  d_13  d_14 ... d_{nLeafs-1, nLeafs}
 //******** d_12  d_13  d_14 ... d_{nLeafs-1, nLeafs}
 //******** .......
 //******** .......
 //******** NOTE: ASSIGNING A NEGATIVE NUMBER (e.g., -1) MEANS DISTANCE UNAVAILABLE -- MISSING ONE OR BOTH LEAFS
 **************/

int main(int argc, const char * argv[]) {
    ifstream inp(argv[1]); // open the argument file
    
    int nSpecies, nLeafs;
    
    inp >> nSpecies >> nLeafs; // read number of species and number of leafs (global)
    
    int *nIndividuals=new int[nSpecies]; // create array to hold the number of individuals per species
    
    for(int i=0;i<nSpecies;i++){
        inp >> nIndividuals[i]; // read number of individuals per species
    }
    
    float maxDist, epsilon=1000000; // disregard epsilon; maxDist is there to have an upper bound on the maximum possible distance for unknown leafs
    
    float *max=new float[nSpecies]; // the max array (per species) holds the biggest possible "disagreement" between distances within the species
    
    maxDist=-1;
    
    float ****d=new float***[nSpecies];
    for(int i=0;i<nSpecies;i++){
        max[i]=-1;
        d[i]=new float**[nIndividuals[i]];
        for(int j=0;j<nIndividuals[i];j++){
            d[i][j]=new float*[nLeafs];
            for(int k=0;k<nLeafs;k++){
                d[i][j][k]=new float[nLeafs];
                for(int l=k+1;l<nLeafs;l++){
                    inp >> d[i][j][k][l]; // read distance per species per individual per leafs
                    if(maxDist<d[i][j][k][l]){
                        maxDist=d[i][j][k][l]; // find maximum distance currently available (for upper bounding)
                    }
                }
            }
        }
    }
    
    for(int i=0;i<nSpecies;i++){
        for(int j1=0;j1<nIndividuals[i];j1++){
            for(int j2=j1+1;j2<nIndividuals[i];j2++){
                for(int k=0;k<nLeafs;k++){
                    for(int l=k+1;l<nLeafs;l++){
                        if(d[i][j1][k][l]>=0 && d[i][j2][k][l]>=0){
                            if((d[i][j1][k][l]-d[i][j2][k][l])*(d[i][j1][k][l]-d[i][j2][k][l])>=max[i]){
                                max[i]=(d[i][j1][k][l]-d[i][j2][k][l])*(d[i][j1][k][l]-d[i][j2][k][l]); // find maximum within each species (in lieu of epsilon)
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Full model (no decomposition)
    
    try{
        GRBEnv env; // define environment
        GRBModel model(env); // define model within environment
        
        
        // Definitions of quadratic objective and quadratic constraints
        GRBQuadExpr Obj;
        GRBQuadExpr QConstr;
        
        
        // Dynamically define variable $x$
        GRBVar  ****x=new GRBVar***[nSpecies];
        for(int i=0;i<nSpecies;i++){
            x[i]=new GRBVar**[nIndividuals[i]];
            for(int j=0;j<nIndividuals[i];j++){
                x[i][j]=new GRBVar*[nLeafs];
                for(int k=0;k<nLeafs;k++){
                    x[i][j][k]=new GRBVar[nLeafs];
                    for(int l=k+1;l<nLeafs;l++){
                        char name[20];
                        sprintf(name, "x(%d, %d, %d, %d)", i, j, k, l);
                        x[i][j][k][l]=model.addVar(0, maxDist, 0, 'C', name); // add the new variable to the model and assign name in the format of x(i,j,k,l)
                    }
                }
            }
        }
        
        model.update(); // updating model to contain the variables
        
        for(int i=0;i<nSpecies;i++){
            for(int j=0;j<nIndividuals[i];j++){
                for(int k=0;k<nLeafs;k++){
                    for(int l=k+1;l<nLeafs;l++){
                        if(d[i][j][k][l]>=0){
                            model.addConstr(x[i][j][k][l]==d[i][j][k][l]); // if distance is available, x[i][j][k][l] is SET
                        }
                    }
                }
            }
        }
        
        for(int i=0;i<nSpecies;i++){
            for(int j1=0;j1<nIndividuals[i];j1++){
                for(int j2=j1+1;j2<nIndividuals[i];j2++){
                    for(int k=0;k<nLeafs;k++){
                        for(int l=k+1;l<nLeafs;l++){
                            model.addQConstr((x[i][j1][k][l]-x[i][j2][k][l])*(x[i][j1][k][l]-x[i][j2][k][l])<=max[i]/*epsilon*/); // within the species the max 2-norm is constrained -- no constraint outside species
                        }
                    }
                }
            }
        }
        
        for(int i1=0;i1<nSpecies;i1++){
            for(int i2=i1+1;i2<nSpecies;i2++){
                for(int j1=0;j1<nIndividuals[i1];j1++){
                    for(int j2=0;j2<nIndividuals[i2];j2++){
                        for(int k=0;k<nLeafs;k++){
                            for(int l=k+1;l<nLeafs;l++){
                                Obj+=(x[i1][j1][k][l]-x[i2][j2][k][l])*(x[i1][j1][k][l]-x[i2][j2][k][l]); // setting the objective to be the summation of 2-norms for two different species (the distance within the species is constrained from before)
                            }
                        }
                    }
                }
            }
        }
        
        model.setObjective(Obj); // setting the objective function
        model.update(); // updating the model
        model.optimize(); // optimizing the model
        
        char fname[20];
        sprintf(fname, "FINAL%s", argv[1]);
        ofstream out(fname);
        out.precision(3); // set output precision to 2 fractional digits (for simplicity)
        
        for(int i=0;i<nSpecies;i++){
            for(int j=0;j<nIndividuals[i];j++){
                for(int k=0;k<nLeafs;k++){
                    for(int l=k+1;l<nLeafs;l++){
                        out << x[i][j][k][l].get(GRB_DoubleAttr_X) << "\t"; // output the distances in file
                    }
                }
                out << endl;
            }
            out << endl << endl;
        }
        
        
    }catch(GRBException e){ // catching exception of solver and printing error code and error message
        cout << e.getErrorCode() << " " << e.getMessage() << endl;
    }catch(...){
        cout << "Error during optimization.";
    }
    
    
    return 0;
}
