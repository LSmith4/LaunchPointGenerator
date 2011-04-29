#include <iostream> 
#include <list>

#include "ModelInfo.hpp"
#include "Point.hpp"

using namespace std;

namespace LPG{

#ifndef FEASIBILITYVECTOR_HPP
#define FEASIBILITYVECTOR_HPP
	
    class FeasibilityVector{
        private:
            double *v;                          //the actual vector
            int size;                           //number of elements
            double lengthSqr;			//the norm2 length squared aka ||feasibility distance||^2
            int id;                             //the constraint number
        public:
            FeasibilityVector(int s);           //create Fv of size s
            ~FeasibilityVector();
            //FeasibilityVector(const FeasibilityVector &Fv);			//copy constructor
            //FeasibilityVector& operator=(const FeasibilityVector &Fv);	//assignment operator

            friend class ConstraintConsensus;   // constraint consensus is a friend class to access v in ConstraintConsensus::updateCounters()

            double getLengthSqr();              //return the norm2 length squared aka ||feasibility distance||^2
            void clear();                       //set all values to 0
            double calc(Point *P, int i);	//calculate Fv to constraint i from Point P
            int getId();                        //return the constraint number Fv was calculated for, -1 if not set
    };


	
#endif
};
