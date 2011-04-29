#include <cmath> 	//pow
#include "FeasibilityVector.hpp"

using namespace std;
using namespace LPG;

FeasibilityVector::FeasibilityVector(int s){

    //cout << "Fv constructor" << endl;
    id = -1;                            //default to -1
    size = s;
    //cout << "size: " << size << endl;
    v = new double[size];               //create the actual vector
    clear();                            //set all elements to zero
    //cout << "Fv constructor end" << endl;
};

FeasibilityVector::~FeasibilityVector(){

    //cout << "Fv deconstructor" << endl;
    delete [] v;
    v=NULL;
    //cout << "Fv deconstructor end" << endl;
};

double FeasibilityVector::getLengthSqr(){

    return lengthSqr;
};

void FeasibilityVector::clear(){
    
    lengthSqr=0.0;              //set the length to zero
    for(int i=0; i<size; i++)
            v[i]=0;             //set all elements to zero
    //for(int i=0; i<size; i++)
    //	cout << Fv[i] << endl;
};

double FeasibilityVector::calc(Point *P, int i){

    //standard linear feasibility vector calculation
    //fv = constraint violation * gradient of constraint / ||gradient of constraint||^2

    id=i;               //set the id number
    cgrad *cg;          //iterator
    clear();            //start clean
    double tempLength;	//sum of (jac values squared)

    tempLength=0.0;     //initialize to zero

    for(cg=P->m_pModel->asl->i.Cgrad_[i]; cg; cg = cg->next)
    {
            v[cg->varno]=P->getConVio(i)*P->getJac(cg->goff);
            tempLength+=pow(P->getJac(cg->goff),2);
    }
    for (cg=P->m_pModel->asl->i.Cgrad_[i]; cg; cg = cg->next)
    {
            v[cg->varno]=v[cg->varno]/tempLength;
            lengthSqr+=pow(v[cg->varno],2);
    }
    return lengthSqr;
};

int FeasibilityVector::getId(){

    if(id>-1)
	return id;
    else{
	cout << "Error: FeasibilityVector::getId(): Fv not calculated yet" << endl;
	return -1;
    }
};
