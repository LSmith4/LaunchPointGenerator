#include <cmath> //abs, pow
#include "Cluster.hpp"

using namespace std;
using namespace LPG;


Cluster::Cluster(Point P){
//cout << "cluster constructor" << endl;
    lPoints.push_back(P);
    lpIndex=0;
};

Cluster::~Cluster(){
//cout << "cluster destructor" << endl;

};

Cluster::Cluster(const Cluster &C){

    //cout << "cluster copy constructor" << endl;

    lPoints = C.lPoints;
    lpIndex = C.lpIndex;

    //cout << "cluster copy constructor end" << endl;
};

Cluster& Cluster::operator=(const Cluster &C){
    //cout << "cluster assignment operator" << endl;
    if(this != &C)
    {
        lPoints = C.lPoints;
        lpIndex = C.lpIndex;
    }
    //cout << "cluster assignment operator end" << endl;
    return *this;
};

void Cluster::addPoint(Point P){
    //cout << "add point " << lpIndex << " " << lPoints.size() << endl;
    lPoints.push_back(P);
    if(lPoints[lpIndex].morePromising(P))
        lpIndex=lPoints.size()-1;
};

Point Cluster::getLaunchPoint(){

    if(lpIndex>lPoints.size() || lpIndex<0){
        cout << "lpIndex out of Range!" << endl;
        exit(0);
    }
  return lPoints[lpIndex];
};