#include <iostream> 
#include <vector>

#include "ModelInfo.hpp"
#include "Point.hpp"
#include "FeasibilityVector.hpp"
#include "TimeKeeper.hpp"

using namespace std;

namespace LPG{

#ifndef CONSTRAINTCONSENSUS_HPP
#define CONSTRAINTCONSENSUS_HPP
	
    class ConstraintConsensus{		//abstract class
        private:
            vector<Point> lPoints;	//list of points, first is provided as aurgument
            vector<double> lTimes;	//list of times in secs, elapsed time to each point
            double alpha;		//feasibility tolerance
            double beta;		//movement tolerance
            int gamma;			//augmentation period
            int mu;			//maximum iterations
            int cType;                  //type of consensus
            double *s;			//counter
            int *n;			//counter
            double *cv;			//consensus vector
            double lengthSqr;		//length^2 of cv
            double maxTime;             //maximum run time for CC
        public:
            ConstraintConsensus(double a,   //feasibility tolerance
                    double b,               //movement tolerance
                    int c,                  //maximum iterations
                    int g,                  //augmentation period
                    int t,                  //concensus type
                    ModelInfo *pModel,      //model info
                    double mt);             //max run time
            ~ConstraintConsensus();                                         //virtual so it gets executed
            //ConstraintConsensus(const ConstraintConsensus &CC);           //copy constructor
            //ConstraintConsensus& operator=(const ConstraintConsensus &CC);//assignment operator

            void printPoints();                         //output list of points
            int Run(Point &P);                          //executes main loop from P
            void updateCounters(FeasibilityVector *Fv,  //feasibilty vector
                Point *P);                              //current point
            double calcCv(int n);                       //calculate the consensus vector using counters
                                                        //number of variables, return lengthSqr
            void getBestPoint(Point *P);                //return the point in the list in terms of MaxVio
            double getTotalTime();                      //return the total time taken by CC run
    };


	
	
#endif
};
