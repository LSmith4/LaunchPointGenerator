#include <iostream>
#include <vector>

#include "Point.hpp"

using namespace std;

namespace LPG{

#ifndef CLUSTER_HPP
#define CLUSTER_HPP

    class Cluster{
        private:
            vector<Point> lPoints;                  //list of points in cluster
            int lpIndex;                            //index of launch point

        public:
            Cluster(Point P);                       //makes a point, initializes
            ~Cluster();                             //destructor
            Cluster(const Cluster &C);              //copy constructor
            Cluster& operator=(const Cluster &C);   //assignment operator

            friend class ClusterBuilder;

            void addPoint(Point P);                  //add a point to the cluster
            Point getLaunchPoint();                 //get the launch point 

    };

#endif

};
