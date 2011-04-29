#include <iostream> 

#include "ModelInfo.hpp"

using namespace std;

namespace LPG{

#ifndef POINT_HPP
#define POINT_HPP

    const double FeasThres = 1e-6;              //point is feasible if max violation is less than FeasThres
    typedef long fint;
	
    class Point{
        private:
            double *x;                          //location
            bool feasible;                      //in feasible region?
            double *c;                          //conval return values
            double *v;                          //violation of each constraint
            double MaxVio;                      //max violation
            double SumVio;                      //sum of absolute violations
            double *j;                          //jacval return values
            
            void refresh();                     //refresh violation data
            void setToBounds();			//reset out of bound variables

        public:
            Point(ModelInfo *pModel);           //makes a point, initializes
            ~Point();                           //destructor
            Point(const Point &P);              //copy constructor
            Point& operator=(const Point &P);	//assignment operator
            bool operator<(const Point &P) const;   //for sorting points

            ModelInfo *m_pModel;		//copy of model pointer

            void setLocation(double *new_x);	//set point to specific location
            double getMaxVio();                 //return L0 constraint violation
            double getSumVio();                 //return L1 constraint violation
            double getConVio(int i);            //get violation of constraint i
            double getJac(int i);               //get jacobian value at offset i
            int addVec(double *nx, int n);      //add nx to x component wise, n is number of elements of nx
            void randLocation();		//set to random location in search space
            void getLocation(double *loc);      //copy location to loc
            double getDist(Point *P);           //return the distance to other point
            bool morePromising(Point P);        //is P more promising?
            string getLocation();               //returns location in string format for output to file
    };

#endif
	
};
