#include <cmath> 	//pow
#include <float.h>      //DBL_MAX

#include "ConstraintConsensus.hpp"


using namespace std;
using namespace LPG;

ConstraintConsensus::ConstraintConsensus(double a, double b, int c, int g, int t, ModelInfo *pModel, double mt){
	
    //set the basic CC variables
    alpha=a;
    beta=b;
    mu=c;
    gamma=g;
    maxTime=mt;
    cType=t;

    int nvar = pModel->n_var;
    n=new int[nvar];            //initialize counters
    s=new double[nvar];
    cv=new double[nvar];        //and consensus vector
    for(int i=0; i<nvar; i++){
        n[i]=0;
        s[i]=0.0;
        cv[i]=0.0;
    }

    lPoints.clear();            //clear list, size now zero
    //lPoints.push_back(P);     //add initial point as first element
    lTimes.clear();             //clear list, first entry added after first point calculated
};

ConstraintConsensus::~ConstraintConsensus(){
	
    delete [] n;
    n=NULL;
    delete [] s;
    s=NULL;
    delete [] cv;
    cv=NULL;
};


int ConstraintConsensus::Run(Point &S){

    //cout << "ConstraintConsensus: Run" << endl;

    lPoints.clear();
    lPoints.push_back(S);
    lTimes.clear();

    TimeKeeper Tk;                  //start measuring elapsed time
    lTimes.push_back(0.0);          //give inital point time zero
    
    Point *P = &lPoints.back();     //get the initial point
    int nvar=P->m_pModel->n_var;    //store number of variables locally
    int ncon=P->m_pModel->n_con;    //store number of constraints locally
    int n_inf;                      //number of infeasible constraints
    FeasibilityVector Fv(nvar);     //create feasibiltiy vector instance
    FeasibilityVector *pFv;         //ptr to feasibility vector
    pFv = &Fv;                      //assign pointer
	

    if(P->getMaxVio()<FeasThres){   //check if initial point is feasible
        //cout << "initial point is feasible" << endl;
        return 0;                   //return successfully without doing anything more
    }
    else{
        //cout << "MaxVio of initial point is: " << P->getMaxVio() << endl;
    }
		
    for(int u=0; u<mu; u++){        //MAIN LOOP: repeat mu times

        //cout << "u: " << u << endl;

        n_inf=0;                    //initialize loop
        if(u>0)                     //if not initial loop
            P=&lPoints.back();      //get the most recent point

        for(int i=0; i<P->m_pModel->nlc; i++){      //for every nonlinear constraint
            if(abs(P->getConVio(i))>0.0){           //if constraint i violated
                if(Fv.calc(P,i) > pow(alpha,2)){    //if feasibility distance > than feasibility tolerance
                    n_inf++;                        //constraint is infeasible at current point
                    updateCounters(pFv,P);
                }
            }
        }
        if(n_inf==0){                   //no violated constraints?
            //cout << "No violated nonlinear constraints!" << endl;
            return 0;                   //then exit successfully
        }
        if(calcCv(nvar)<pow(beta,2)){   //||cv||^2 < beta^2, is the concensus vector less than the movement tolerance
            return 1;			//return unsuccessfully due to lack of movement
        }
        else{
            //cout << "make copy of point" << endl;
            Point nP(*P);		//make copy of current point
            if(nP.addVec(cv, nvar)>0){	//add the consensus vector to the new point
                //cout << "Error: Basic::Run(): calculating new iterate" << endl;
                return 1;
            }
            //cout << "maxVio: " << nP.getMaxVio() << " time: " << Tk.getElapsedTimeSec() << endl;
            lPoints.push_back(nP);                      //add new point to list
            lTimes.push_back(Tk.getElapsedTimeSec());
            if(Tk.getElapsedTimeSec()>maxTime){
                //cout << "CC out of time!" << endl;
                return 3;
            }
        }
    }
	
    //cout << "Run() end" << endl;
    //cout << "CC out of iterations!" << endl;
    return 2;   //did not exit within main loop
};

void ConstraintConsensus::updateCounters(FeasibilityVector *Fv, Point *P){
	
    //cout << "updateCounters" << endl;
    cgrad *cg;			//iterator

    //update this for each type of constraint consensus!
    switch(cType){
        case 1:
            //cout << "Consensus: BASIC" << endl;
            //for every variable in constraint i
            for (cg=cgStart(P,Fv); cg; cg = cg->next){
		n[cg->varno]++;
		s[cg->varno]=s[cg->varno]+Fv->v[cg->varno];
            }
            break;
        case 2:
            //cout << "Consensus: SUM" << endl;
            //for every variable in constraint i
            for (cg=cgStart(P,Fv); cg; cg = cg->next){
		n[cg->varno]++;
		s[cg->varno]=s[cg->varno]+Fv->v[cg->varno];
            }
            break;
        default:
            cout << "Consensus out of range!" << endl;
            cout << "This should never happen! If it does there is a serious problem." << endl;
            exit(0);
            break;
    }
    //cout << "updateCounters end" << endl;
};

double ConstraintConsensus::calcCv(int nvar){

    //cout << "calcCv" << endl;
    lengthSqr=0;                            //reset to zero
    //update this for each type of constraint consensus!
    switch(cType){
        case 1:
            //cout << "Consensus: Basic" << endl;
            for(int j=0; j<nvar; j++){      //for all variables
                if(n[j]!=0){                //if variable part of sub-set of violated constraints
                    cv[j]=s[j]/n[j];        //calculate consensus vector component
                    lengthSqr+=pow(cv[j],2);//update length of cv
                }
                else                        //variable not part of sub-set of violated constraints
                    cv[j]=0;                //variable not part of consensus
            }
            break;
        case 2:
            //cout << "Consensus: Basic" << endl;
            for(int j=0; j<nvar; j++){      //for all variables
                if(n[j]!=0){                //if variable part of sub-set of violated constraints
                    cv[j]=s[j];        //calculate consensus vector component
                    lengthSqr+=pow(cv[j],2);//update length of cv
                }
                else                        //variable not part of sub-set of violated constraints
                    cv[j]=0;                //variable not part of consensus
            }
            break;
        default:
            cout << "Consensus out of range!" << endl;
            cout << "This should never happen! If it does there is a serious problem." << endl;
            break;
    }
    //cout << "calcCv end" << endl;
    return lengthSqr;                       //return the squared length
};

void ConstraintConsensus::printPoints(){

    cout << endl;
    cout << " - print points - " << endl;
    cout << endl;
    int i=0;
    vector<double>::iterator lj = lTimes.begin();
    cout << "point \t MaxVio \t SumVio \t time(s)" << endl;
    for (vector<Point>::iterator li = lPoints.begin(); li != lPoints.end();++li){
        cout << i << " \t " << li->getMaxVio() << " \t " << li->getSumVio() << " \t " << *lj << endl;
        i++;
        ++lj;
    }
};

void ConstraintConsensus::getBestPoint(Point *P){

    //Point B(lPoints.begin()->m_pModel);
    double MinVio = DBL_MAX;
    for (vector<Point>::iterator li = lPoints.begin(); li != lPoints.end();++li){
        if(li->getMaxVio() < MinVio){
            MinVio=li->getMaxVio();
            //cout << MinVio << endl;
            *P=*li;
            //cout << P->getMaxVio() << endl;
        }
    }
    
};

double ConstraintConsensus::getTotalTime(){

  return lTimes[lTimes.size()-1];
};