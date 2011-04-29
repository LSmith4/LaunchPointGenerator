#include <iostream>
#include <vector>
#include <fstream>

#include "ModelInfo.hpp"
#include "Point.hpp"
#include "ClusterBuilder.hpp"
#include "Cluster.hpp"

using namespace std;

namespace LPG{

#ifndef MATLABSCRIPTWRITER_HPP
#define MATLABSCRIPTWRITER_HPP


    class MatlabScriptWriter{
        private:
            bool write;             //to write or not
            string MatlabScriptDir; //the directory to put files
        public:
            MatlabScriptWriter(ModelInfo *pModel);
            ~MatlabScriptWriter();

            ModelInfo *m_pModel;    //copy of model pointer

            void writePlotPoints(vector<Point> *lPoints, string fileTitle); //write PlotPoints file
            void writePlotFreqDist(vector<Freq> F, string fileTitle);       //write FreqDist file
    };

#endif

};
