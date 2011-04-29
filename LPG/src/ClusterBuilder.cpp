#include <vector>
#include <iomanip>
#include <float.h> //DBL_MAX
#include <math.h> //floor

#include "ClusterBuilder.hpp"
#include "Point.hpp"
#include "Cluster.hpp"

using namespace std;
using namespace LPG;

ClusterBuilder::ClusterBuilder(vector<Point> *lPoints, int windowSize, int maxLaunches){

    cout << "\n**************************************************************************\n";
    cout << "\tCreate the inter-point distance distribution.\n";
    
    Tk.reset();                     //set the timer
    time = 0.0;                     //initialize time
    w = windowSize;
    tau = maxLaunches;
    allPoints.clear();
    for (vector<Point>::iterator li = lPoints->begin(); li != lPoints->end();++li){
        allPoints.push_back(*li);
    }
    promPeaks.clear();
    calcDistances();                //calcualte the distances between each point, dmin, dmax, and dwidth
    calcFrequencies();              //calculate the interpoint frequency distribution

    cout << endl;
    cout << "\n**************************************************************************\n";
    cout << "\tCalculate prominent peaks and find clusters.\n";

    Tk.reset();                     //reset timer

    //algorithm 1 in thesis
    bool finished = false;
    int numClusters = 0;
    for(w; w>0; w--){               //try different window sizes
        extractPP();                //using w, extract the prominent peaks from the interpoint frequency distribution

        
        bool write = true;
        
        for(int i=0; i<promPeaks.size()-1; i++){
            if(i==0)
                cDist=(promPeaks[0].d+dmin)/2;              //initial minimum
            else
                cDist=(promPeaks[i-1].d+promPeaks[i].d)/2;  //minimum between peaks
            numClusters = cluster();                        //calcualte the clusters

            if(numClusters>0 && write){
                cout    << "\tIdentify Clusters:\n";
                cout    << "\t" << setw(14)   << "Window Size"
                        << setw(22)   << "Critical Distance"
                        << setw(18)  << "# of Clusters"
                        << endl;
                cout    << "\t-------------------------------------------------------\n";
                write=false;
            }
  
            cout    << "\t" << setw(14) << w
                            << setw(22) << cDist
                            << setw(18) << numClusters
                            << endl;

            if(numClusters<=tau && numClusters>0){          //number of clusters meets criteria?
                finished=true;                              //yes: exit loop
                break;
            }                                               //no: try next critical distance
        }
        if(finished)                                        //was w successful?
            break;                                          //yes: exit loop
    }                                                       //no: decrement w, try again
    if(!finished){
        cout << "ClusterBuilder Failed! w==" << w << endl;
        exit(0);
    }
    cout    << "\t-------------------------------------------------------\n";
    cout << "\t"    << "\t Total time to find clusters: " << Tk.getElapsedTimeSec() << endl;
    cout << "\t\t " << numClusters << " clusters found with critical distance: " << cDist << endl;
    cout << endl;
    time+=Tk.getElapsedTimeSec();


};

ClusterBuilder::~ClusterBuilder(){

};

double ClusterBuilder::getTime(){
    return time;
};

void ClusterBuilder::calcDistances(){

    Dist dtemp;
    dmin = DBL_MAX;
    dmax = 0.0;
    for(int i=0; i<allPoints.size()-1; i++){
        for(int j=i+1; j<allPoints.size(); j++){

            dtemp.dist=allPoints[i].getDist(&allPoints[j]);
            dtemp.p1=i;
            dtemp.p2=j;

            d.push_back(dtemp);             //add distance struct to vector

            if(dtemp.dist<dmin)
                dmin=dtemp.dist;            //the minimum distance in the vector
            if(dtemp.dist>dmax)
                dmax=dtemp.dist;            //the maximum distance in the vector
        }
    }
    dwidth=(dmax-dmin)/allPoints.size();    //the width of each bin

    cout << endl;
    cout << "\t- Distribution Details - " << endl;
    cout << "\t---------------------- " << endl;
    cout << "\t      dmin: " << dmin << endl;
    cout << "\t      dmax: " << dmax << endl;
    cout << "\t    dwidth: " << dwidth << endl;
};

void ClusterBuilder::calcFrequencies(){

    Freq temp;
    temp.F=0;
    temp.d=0.0;
    for(int j=0; j<allPoints.size(); j++){
        F.push_back(temp);
    }

    double lb = dmin;
    double ub = lb+dwidth;
    for(int i=0; i<F.size(); i++){
        lb=dmin+i*dwidth;
        ub=lb+dwidth;
        F[i].d=dmin+(i+0.5)*dwidth;
        //cout << "lb: " << lb << endl;
        //cout << "ub: " << ub << endl;
        for(int j=0; j<d.size(); j++){
            if(d[j].dist<ub && d[j].dist>lb)
                F[i].F++;
        }
    }

    cout    << endl;
    cout    << "\t" << setw(25)  << "Bin Distance"
                    << setw(3)   << ""
                    << setw(14)  << "Frequency"
                    << endl;
    cout    << "\t" << setw(14)  << "Lower"
                    << setw(14)  << "Upper"
                    << setw(14)  << ""
                    << endl;
    cout    << "\t" << "------------------------------------------------------\n";
    for (int j=0; j<F.size(); j++){
        cout    << "\t" << setw(14)  << setiosflags(ios::fixed) <<setprecision(3) << F[j].d-dwidth/2
                        << setw(14)  << setiosflags(ios::fixed) <<setprecision(3) << F[j].d+dwidth/2
                        << setw(14) << setiosflags(ios::fixed) <<setprecision(3) << F[j].F
                        << endl;
    }

    cout    << "\t" << "------------------------------------------------------\n";
    time = Tk.getElapsedTimeSec();
    cout    << "\t" << "\t Total time creating distribution: " << time << endl;

    //for(int j=0; j<F.size(); j++){
    //    cout << "F[" << j << "]=" << F[j] << endl;
    //}
    
};

vector<Freq> ClusterBuilder::getFreq(){

    if(F.size()<2){
        cout << "ClusterBuilder: F isn't defined yet!" << endl;
        exit(0);
    }
    return F;
};

vector<Cluster> ClusterBuilder::getClusters(){

    if(lClusters.size()<1){
        cout << "ClusterBuilder: no clusters available yet!" << endl;
        exit(0);
    }
    return lClusters;
};

vector<Point> ClusterBuilder::getLaunchPoints(){
  
    vector<Point> lp;
    lp.clear();
    for(int i=0; i<lClusters.size(); i++){
        lp.push_back(lClusters[i].getLaunchPoint());
    }
    if(lp.size()==0)
        cout << "There are no launch points!" << endl;
    sort(lp.begin(), lp.end());
    return lp;
};

void ClusterBuilder::extractPP(){

    promPeaks.clear();
    PromPeak tempPP;
    //cout << "\nDetermine the set of prominent peaks\n";
    //cout << "w: " << w << endl;
    //cout << "number of PP: " << promPeaks.size() << endl;
    
    while(promPeaks.size()==0 && w>0){
        if(2*w>=F.size()){
            w=floor(F.size()/2);
            //cout << "w reduced to: " << w << endl;
        }
        if(w<1){
            cout << "Critical distance extraction fail! w=" << w << endl;
            exit(0);
        }
        for(int i=w; i< F.size()-w; i++){
            bool pp=true;
            for(int j=-w; j<w; j++){
                if(i!=i+j && F[i+j].F>=F[i].F){
                    pp=false;
                    break;
                }
            }
            if(pp==true){
                tempPP.F=F[i].F;
                tempPP.bin=i;
                tempPP.d=F[i].d;
                promPeaks.push_back(tempPP);
            }
        }
        if(promPeaks.size()==0 && w>1){
            w--;
            //cout << "w reduced to: " << w << endl;
        }
        else if(promPeaks.size()==0 && w==1){
            cout << "Failed to extract critical distance!" << endl;
            exit(0);
        }
    }
    //cout << "number of PP: " << promPeaks.size() << endl;

    cout << endl;
    cout << "\tProminent Peaks: w=" << w << endl;
    cout    << "\t" << setw(8)  << "Peak"
                    << setw(14) << "Bin Center"
                    << setw(14) << "Frequency"
                    << endl;
    cout << "\t-------------------------------------" << endl;
    for(int i=0; i<promPeaks.size(); i++){

        cout    << "\t" << setw(8)  << i+1
                        << setw(14) << promPeaks[i].d
                        << setw(14) << promPeaks[i].F
                        << endl;
    }
    cout << endl;


};

int ClusterBuilder::cluster(){

    lClusters.clear();                      //empty list of clusters
    Cluster C(allPoints[0]);                //start a cluster with the first point
    lClusters.push_back(C);                 //add it to the list
    //cout << "Add Cluster to list" << endl;
    for(int k=1; k<allPoints.size(); k++){   //for each point

        //cout << "k: " << k << endl;
        bool assigned = false;
        while(assigned==false){
            //cout << "Enter while loop" << endl;
            for(int i=0; i<lClusters.size(); i++){  //for each cluster
                //cout << "i: " << i << " size: " << lClusters.size() << endl;
                for(int j=0; j<lClusters[i].lPoints.size(); j++){ //for each point in cluster
                    //cout << "j: " << j << " size: " << lClusters[i].lPoints.size() << endl;
                    if(allPoints[k].getDist(&lClusters[i].lPoints[j])<cDist){
                        lClusters[i].addPoint(allPoints[k]);
                        //cout << "Add Point " << k << " to lClusters " << i << endl;
                        assigned=true;
                    }
                    if(assigned)
                        break;
                }
                if(assigned)
                    break;
            }
            if(assigned==false){
                Cluster nC(allPoints[k]);
                lClusters.push_back(nC);
                //cout << "Add Cluster to list" << endl;
                assigned==true;
            }
        }
        //cout << lClusters.size() << endl;
    }
    //double *loc;
    //loc = new double[2];
    //for(int i=0; i<lClusters.size(); i++){
    //    cout << "most promising: " << lClusters[i].lpIndex << endl;
     //   for(int j=0; j<lClusters[i].lPoints.size(); j++){
     //       lClusters[i].lPoints[j].getLocation(loc);
     //       cout << "cluster: " << i << " point: " << loc[0] << ", " << loc[1] << " " << lClusters[i].lPoints[j].getMaxVio() << endl;
    //    }
    //}
    //delete [] loc;

    return lClusters.size();
};