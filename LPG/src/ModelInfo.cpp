#include <float.h>

#include "ModelInfo.hpp"

using namespace std;
using namespace LPG;

ModelInfo::ModelInfo(char *stub){
	
    abound=10000;                           //artificial bound 10^4
    Infinity=DBL_MAX;                       //infinity
    negInfinity=-Infinity;                  //negative infinity
    name=stub;                              //set the name
    asl = ASL_alloc(ASL_read_pfgh);         //allocate asl structure
    nl = jac0dim(stub,(fint)strlen(stub));  //read in model file
    pfgh_read(nl,0);                        //get all the data
};


void ModelInfo::print(){

    cout << "\t - Model Details - " << endl;
    cout << "\t -----------------------" << endl;
    cout << "\t model: " << name << endl;
    cout << "\t n_var: " << n_var << endl;
    cout << "\t n_con: " << n_con << endl;
    cout << "\t   nlc: " << nlc << endl;
};

double ModelInfo::getBound(){
    return abound;
};

double ModelInfo::getInfty(){
    return Infinity;
};

double ModelInfo::getNegInfty(){
    return negInfinity;
};

