#include <vector>

#include "MatlabScriptWriter.hpp"

using namespace std;
using namespace LPG;

MatlabScriptWriter::MatlabScriptWriter(ModelInfo *pModel){

    m_pModel = pModel;
    write=false;
    //cout << "n_var: " << m_pModel->n_var << endl;
    if(m_pModel->n_var == 2){
        write=true;
        //cout << "two variables: make matlab plots!" << endl;
    }
    MatlabScriptDir="../tmp/";
};

MatlabScriptWriter::~MatlabScriptWriter(){

};

void MatlabScriptWriter::writePlotPoints(vector<Point> *lPoints, string fileTitle){

    string dir = MatlabScriptDir;                 //directory to place matlab scripts
    dir.append(fileTitle);                          //append the file title
    //cout << "Title: " << dir << endl;
    if(write==true){                                //only true for 2d models

        ofstream mfile(dir.c_str());                //open stream
        if (mfile.is_open())
        {
            mfile << "%M file created by: MatlabScriptWriter.cpp\n";
            mfile << "%This script plots a list of points\n";
            mfile << "close all;\n";
            mfile << "clc;\n";
            mfile << "clear all;\n";

            mfile << "x=[";                         //define 2d matrix

            double *loc;
            loc = new double[m_pModel->n_var];
            for (vector<Point>::iterator li = lPoints->begin(); li != lPoints->end();++li){
                li->getLocation(loc);
                //cout << loc[0] << " , " << loc[1] << endl;
                mfile << loc[0] << " " << loc[1] << ";  ";
            }

            mfile << "];\n";
            mfile << "scatter(x(:,1), x(:,2));\n";  //scatter plot
            mfile << "axis([" <<m_pModel->LUv[0]<<" "<<m_pModel->LUv[1]<<" "<<m_pModel->LUv[2]<<" "<<m_pModel->LUv[3]<<"]);\n";
            mfile.close();
        }
        else cout << "Unable to open file: \n" << dir;


    }
};

void MatlabScriptWriter::writePlotFreqDist(vector<Freq> F, string fileTitle){

    string dir = MatlabScriptDir;                   //directory to place matlab scripts
    dir.append(fileTitle);                          //append the file title
    //cout << "Title: " << dir << endl;
    //if(write==true){                                //only true for 2d models

        ofstream mfile(dir.c_str());                //open stream
        if (mfile.is_open())
        {
            mfile << "%M file created by: MatlabScriptWriter.cpp\n";
            mfile << "%This script plots the inter-point frequency distribution\n";
            mfile << "close all;\n";
            mfile << "clc;\n";
            mfile << "clear all;\n";

            mfile << "x=[";                         //define 2d matrix

            for(int i=0; i<F.size(); i++){
                mfile << F[i].d << " " << F[i].F << ";  ";
            }

            mfile << "];\n";
            mfile << "bar(x(:,1), x(:,2),'g');\n";  //histogram
            mfile << "axis([min(x(:,1)) max(x(:,1)) min(x(:,2)) max(x(:,2))]);\n";
            mfile.close();
        }
        else cout << "Unable to open file: \n" << dir;


    //}
};

