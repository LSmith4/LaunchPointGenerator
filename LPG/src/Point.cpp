#include <cmath> //abs, pow
#include <sstream>
#include "Point.hpp"

using namespace std;
using namespace LPG;


Point::Point(ModelInfo *pModel){
	
    //cout << "point constructor " << endl;
    m_pModel = pModel;
    x = new double[m_pModel->n_var];    //an element for each variable
    c = new double[m_pModel->n_con];    //an element for each constraint
    v = new double[m_pModel->n_con];    //an element for each constraint
    j = new double[m_pModel->nzc];      //an element for each non zero constraint
    randLocation();                     //set to random location in search space
};

Point::~Point(){

    //cout << "point destructor "<< endl;
    delete [] x;
    x=NULL;
    delete [] c;
    c=NULL;
    delete [] v;
    v=NULL;
    delete [] j;
    j=NULL;
    //cout << "point destructor end "<< endl;
};

Point::Point(const Point &P){
	
    //cout << "copy constructor" << endl;
	
    m_pModel = P.m_pModel;
    x = new double[m_pModel->n_var];
    for(int i=0; i< m_pModel->n_var; i++)
    	x[i]=P.x[i];
    j = new double[m_pModel->nzc];
    for(int i=0; i< m_pModel->nzc; i++)
    	j[i]=P.j[i];
    feasible = P.feasible;
    c = new double[m_pModel->n_con];
    v = new double[m_pModel->n_con];
    for(int i=0; i< m_pModel->n_con; i++){
    	c[i]=P.c[i];
    	v[i]=P.v[i];
    }
    MaxVio = P.MaxVio;		
    SumVio = P.SumVio;
	
    //cout << "copy constructor end" << endl;
};

Point& Point::operator=(const Point &P){
    //cout << "assignment operator" << endl;
    if(this != &P)
    {
        delete [] x;
        delete [] c;
        delete [] v;
        delete [] j;
        m_pModel = P.m_pModel;
        x = new double[m_pModel->n_var];
        for(int i=0; i< m_pModel->n_var; i++)
            x[i]=P.x[i];
        j = new double[m_pModel->nzc];
        for(int i=0; i< m_pModel->nzc; i++)
            j[i]=P.j[i];
        feasible = P.feasible;
        c = new double[m_pModel->n_con];
        v = new double[m_pModel->n_con];
        for(int i=0; i< m_pModel->n_con; i++){
            c[i]=P.c[i];
            v[i]=P.v[i];
        }
        MaxVio = P.MaxVio;
        SumVio = P.SumVio;
    }
    return *this;
}; 

bool Point::operator<(const Point& P) const{

    return (MaxVio < P.MaxVio);
}; 

void Point::randLocation(){
    
    double tempL;                               //temporary lower bonud
    double tempU;                               //temporary upper bound
    for(int i=0; i<m_pModel->n_var; i++){	//for each variable

        if(m_pModel->LUv[2*i] > m_pModel->getNegInfty() && m_pModel->LUv[2*i+1] < m_pModel->getInfty()){
            //both known
            tempU = m_pModel->LUv[2*i+1];       //set upper
            tempL = m_pModel->LUv[2*i];         //set lower
        }
        else if(m_pModel->LUv[2*i] > m_pModel->getNegInfty()){
            //only upper undetermined
            tempL = m_pModel->LUv[2*i];         //set lower
            tempU = tempL+m_pModel->getBound(); //artificial upper
        }
        else if(m_pModel->LUv[2*i+1] < m_pModel->getInfty()){
            //only lower undetermined
            tempU = m_pModel->LUv[2*i+1];       //set upper
            tempL = tempU-m_pModel->getBound(); //artificial lower
        }
        else{
            //both bounds undetermined
            tempL = -m_pModel->getBound();      //artificial lower
            tempU = m_pModel->getBound();       //artificial upper
        }
        
        //calculate random number within bounds
        x[i]=(tempU-tempL)*((double)rand()/(double)RAND_MAX)+tempL;		
        //cout << m_pModel->LUv[2*i] << " " << x[i] << " " << m_pModel->LUv[2*i+1] << endl;
    }
    setToBounds();	//make sure point is in search space
    refresh();		//refresh the violation data
};

void Point::refresh(){
	
    //cout << "refresh" << endl;

    //initialize violations
    MaxVio=0.0;			
    SumVio=0.0;
    feasible=false;		//default to infeasible

    fint nerror;		//used for returning error from asl calculations
    nerror=0;

    //get constraint function values at current point
    conval(x, c, &nerror);
    if(nerror){
        cout << "Error in Point::refresh()::Conval" << endl;
        cout << "Problem is ignored, TODO: something intelligent instead!" << endl;
    }

    //for(int i=0; i<m_pModel->n_con; i++){
    //	cout << m_pModel->LUrhs[2*i] << " " << c[i] << " " << m_pModel->LUrhs[2*i+1] << endl;
    //}

    //calculate constraint violations
    //constraints of form: LUrhs[2*i] < c(i) < LUrhs[2*i+1]
    for(int i=0; i<m_pModel->n_con; i++){

        if(c[i] < m_pModel->LUrhs[2*i])
            v[i]=(double)m_pModel->LUrhs[2*i]-c[i];
        else if(c[i] > m_pModel->LUrhs[2*i+1])
            v[i]=(double)m_pModel->LUrhs[2*i+1]-c[i];
        else
            v[i]=0.0;

        if(MaxVio<abs(v[i]))
            MaxVio=abs(v[i]);           //if exceeded, set L0 violation
        SumVio+=abs(v[i]);		//sum up L1 violation

        //cout << v[i] <<  endl;
    }
    //cout << "Max: " << MaxVio << " Taxi: " << TaxiVio << endl;

    if(MaxVio<FeasThres)
        feasible=true;

    //get jacobian values at current point
    jacval(x,j,&nerror);
    if(nerror!=0){
        cout << "Error in Point::refresh()::Jacval" << endl;
        cout << "Problem is ignored, TODO: something intelligent instead!" << endl;
    }
    //for(int i=0; i<m_pModel->n_var; i++){		//for each variable
    //	cout << j[i] <<  endl;
    //}
    //cout << "refresh end" << endl;
};

void Point::setToBounds(){
	
    //cout << "setToBounds" << endl;
    double badElement;
    for(int i=0; i<m_pModel->n_var; i++){		//for each variable
        badElement=x[i];
        if(x[i] < m_pModel->LUv[2*i]) {
            x[i] = m_pModel->LUv[2*i];			//reset to lower bound if violated
            //cout << "setToBounds(): " << i << "th variable to " << x[i] << " from " << badElement << endl;
        }
        else if(x[i] > m_pModel->LUv[2*i+1]){
            x[i] = m_pModel->LUv[2*i+1];		//reset to upper bound if violated
            //cout << "setToBounds(): " << i << "th variable to " << x[i] << " from " << badElement << endl;
        }
        //cout << m_pModel->LUv[2*i] << " " << x[i] << " " << m_pModel->LUv[2*i+1] << endl;
    }
    //cout << "setToBounds end" << endl;
};

void Point::setLocation(double *new_x){

    for(int i=0; i<m_pModel->n_var; i++){		//for each variable
        x[i]=new_x[i];
        //cout << m_pModel->LUv[2*i] << " " << x[i] << " " << m_pModel->LUv[2*i+1] << endl;
    }
    setToBounds();	//make sure point is in search space
    refresh();		//refresh the violation data
};

double Point::getMaxVio(){
    return MaxVio;
};

double Point::getSumVio(){
    return SumVio;
};

double Point::getConVio(int i){
    if(i<m_pModel->n_con && i > -1)
	return v[i];
    else{
	cout << "Out of range: Point::getConVio()" << endl;
	return 1;
    }
};


double Point::getJac(int i){
    if(i<m_pModel->nzc && i > -1)
	return j[i];
    else{
	cout << "Out of range: Point::getJac()" << endl;
	return 1;
    }
};

int Point::addVec(double *nx, int n){
	
    //cout << "addVec" << endl;
    if(n!=m_pModel->n_var){		//if new vector isn't right size
	cout << "Error: wrong number of elements in nx" << endl;
	return 1;			//return unsuccessfully
    }
    else{
        for(int i=0; i<n; i++)          //add each element
            x[i]+=nx[i];

        setToBounds();			//make sure point is in search space
        refresh();			//refresh the violation data
    }
	
    //cout << "addVec end" << endl;
    return 0;
};

void Point::getLocation(double* loc){

    for(int i=0; i<m_pModel->n_var; i++)
        loc[i]=x[i];
};

double Point::getDist(Point* P){

    double dist = 0;
    double *loc;
    loc = new double[m_pModel->n_var];
    P->getLocation(loc);
    for(int i=0; i<m_pModel->n_var; i++){
        dist+=pow(x[i]-loc[i], 2);
    }
    dist = sqrt(dist);

    delete [] loc;
    return dist;
};

bool Point::morePromising(Point P){

    //returns true if P is more promising than this
    bool mp=false;

    if(this->feasible && P.feasible){

        //both feasible, until obj implemented, current is more promising
    }
    else if(P.feasible)
        mp=true;
    else if(this->getMaxVio()>P.getMaxVio())
        mp=true;

    return mp;
};

string Point::getLocation(){
  stringstream ss;
  ss << "[";
  for(int i=0; i<m_pModel->n_var-1; i++){
     ss << x[i] << ", ";
  }
  ss << x[m_pModel->n_var-1];
  ss << "];" << endl;
  return ss.str();
};