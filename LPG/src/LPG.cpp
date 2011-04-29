#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

#include "LPG.hpp"
#include "ModelInfo.hpp"
#include "Point.hpp"
#include "ConstraintConsensus.hpp"
#include "MatlabScriptWriter.hpp"
#include "ClusterBuilder.hpp"
#include "TimeKeeper.hpp"

using namespace std;
using namespace LPG;

int main(int argc, char* argv[])
{
    TimeKeeper Tk_all;                  //keep track of total time
    TimeKeeper Tk_util;                 //used periodically
    double CCTime;                      //constraint consensus time
    double CBTime;                      //cluster builder time
    double LPTime;                      //calculating launch point time

    cout << endl;
    cout << "\n**************************************************************************\n";
    cout << "CCopt - Generates a list of high-quality launch points for a local solver.\n";
    cout << "\n";

    srand(time(NULL));              //initialize random seed

    if(argv[1]==NULL){
        cout << "\nNo model supplied!" << endl;
        cout << "Exiting CCopt\n" << endl;
        return 1;
    }
    stub=argv[1];                   //get model name from command line
    ModelInfo model(stub);          //create model
    ModelInfo *pModel = &model;     //pointer
    pModel->print();                //print the basic model info


    //default parameter values
    double alpha = 0.00001;         //CC: feasibility tolerance
    double beta = 0.001;            //CC: movement tolerance
    int mu = 100;                   //CC: maximum iterations
    int CCtype = 2;                 //CC: concensus type 1=basic, 2=sum
    double maxTime = 0.01;          //CC: max run time
    int gamma =0;                   //CC: augmentaion period !!!not implemented yet
    int w = 3;                      //ClusterBuilder: window size
    int tau = 25;                   //ClusterBuilder: max number of clusters
    int p = 150;                     //number of initial sample points

    cout << endl;
    cout << "\t - CC Details - " << endl;
    cout << "\t ----------------------" << endl;
    cout << "\t     alpha: " << alpha << endl;
    cout << "\t      beta: " << beta << endl;
    cout << "\t        mu: " << mu << endl;
    cout << "\t   maxTime: " << maxTime << endl;
    switch(CCtype){
        case 1:
            cout << "\t Consensus: BASIC" << endl;
            break;
        case 2:
            cout << "\t Consensus: SUM" << endl;
            break;
        default:
            cout << "Error: Consensus out of range!\n";
            exit(0);
            break;
    }
    cout << endl;
    cout << "\t - Clustering Details - " << endl;
    cout << "\t -----------------------" << endl;
    cout << "\t   w: " << w << endl;
    cout << "\t tau: " << tau << endl;
    cout << "\t   p: " << p << endl;
    cout << endl;

    cout << "\n**************************************************************************\n";
    cout << "\tRun Constraint Consensus from each start point\n";
    cout << endl;

    //create a list of start points
    vector<Point> sPoints;
    sPoints.clear();
    Point X(pModel);
    for(int i=0; i<p; i++){
        X.randLocation();
        sPoints.push_back(X);
    }
    
    //run constraint consensus from start points
    vector<Point> bPoints;
    bPoints.clear();
    vector<double> bTimes;
    bTimes.clear();
    ConstraintConsensus CC(alpha, beta, mu, gamma, CCtype, pModel, maxTime);
    Point B(pModel);
    for (vector<Point>::iterator li = sPoints.begin(); li != sPoints.end();++li){
        CC.Run(*li);
        CC.getBestPoint(&B);
        bPoints.push_back(B);
        bTimes.push_back(CC.getTotalTime());
    }

    cout    << "\t" << setw(8)  << "CC Run"
                    << setw(22)  << "Violation"
                    << setw(20)  << "Time(s)"
                    << endl;
    cout    << "\t" << setw(8)  << ""
                    << setw(14)  << "Start"
                    << setw(14)  << "End"
                    << setw(14)  << ""
                    << endl;
    cout    << "\t" << "------------------------------------------------------\n";


    CCTime=0.0;
    for (int i=0; i<bPoints.size(); i++){
        cout    << "\t" << setw(8)  << setiosflags(ios::fixed) <<setprecision(3) << i
                        << setw(14) << setiosflags(ios::fixed) <<setprecision(3) << sPoints[i].getMaxVio()
                        << setw(14) << setiosflags(ios::fixed) <<setprecision(3) << bPoints[i].getMaxVio()
                        << setw(14) << setiosflags(ios::fixed) <<setprecision(3) << bTimes[i] << endl;
        CCTime+=bTimes[i];
    }
    cout    << "\t" << "------------------------------------------------------\n";
    cout    << "\t" << "\t Total time running CC: " << CCTime << endl;


    //write m files for plotting the start and end points
    MatlabScriptWriter MLS(pModel);
    MLS.writePlotPoints(&sPoints, "startPoints.m");
    MLS.writePlotPoints(&bPoints, "bestPoints.m");

    CBTime=0.0;
    Tk_util.reset();
    //create clusters
    ClusterBuilder CB(&bPoints, w, tau);
    CBTime = Tk_util.getElapsedTimeSec();
        
    //write m files for the inter-points frequency distribution and launch points
    MLS.writePlotFreqDist(CB.getFreq(), "freqDist.m");
    MLS.writePlotPoints(&CB.getLaunchPoints(), "launchPoints.m");

    LPTime=0.0;
    Tk_util.reset();
    cout << "\n**************************************************************************\n";
    cout << "\tCalculate and sort the launch points.\n";
    cout << endl;
    cout << "\t" << setw(14) << "Point"
                 << setw(14) << "Violation"
                 << endl;
    cout << "\t---------------------------------\n";
    for(int i=0; i<CB.getLaunchPoints().size(); i++){
        cout << "\t"    << setw(14) << i+1
                        << setw(14) << CB.getLaunchPoints()[i].getMaxVio()
                        << endl;

    }
    cout << endl;
    cout << "\t\t There are " << CB.getLaunchPoints().size() << " suggested launch points" << endl;
    LPTime=Tk_util.getElapsedTimeSec();

    cout << "\n**************************************************************************\n";
    cout << endl;
    cout << "\t Time Summary" << endl;
    cout << "\t---------------------------------------------" << endl;
    cout << "\t\t      Constraint Consensus: " << CCTime << endl;
    cout << "\t\t           Cluster Builder: " << CBTime << endl;
    cout << "\t\t Calculating Launch Points: " << LPTime << endl;
    cout << "\t---------------------------------------------" << endl;
    cout << "\t\t                Total Time: " << Tk_all.getElapsedTimeSec() << endl;
    
    //write launch points out to file
    string outFile = "../tmp/LaunchPoints.data";
    ofstream mfile(outFile.c_str());                //open stream
    if (mfile.is_open())
    {
        for(int i=0; i<CB.getLaunchPoints().size(); i++){
            //to console
            //cout << "Point" << i << " = " << CB.getLaunchPoints()[i].getLocation() << endl;
            //to file
            mfile << "Point" << i << " = " << CB.getLaunchPoints()[i].getLocation() << endl;
        }
        mfile.close();
    }
    else cout << "Unable to open file: \n" << outFile;
    cout << endl;
    cout << "\t Wrote launch points to file: " << outFile << endl;
    cout << "\n**************************************************************************\n";
    

    return 0;
};

