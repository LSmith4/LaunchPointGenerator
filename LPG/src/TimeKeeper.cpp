#include "TimeKeeper.hpp"

using namespace std;
using namespace LPG;

TimeKeeper::TimeKeeper(){

    reset();
	//getrusage(RUSAGE_SELF, &ru);
	//startTime=ru.ru_utime;		//set start time
	
};

double TimeKeeper::getElapsedTimeSec(){

    getrusage(RUSAGE_SELF, &ru);
    tempTime=ru.ru_utime;		//get time now

    double tS = startTime.tv_sec + (startTime.tv_usec)*1e-6;
    double tE = tempTime.tv_sec  + (tempTime.tv_usec)*1e-6;

    return tE-tS;

};

void TimeKeeper::reset(){

    getrusage(RUSAGE_SELF, &ru);
    startTime=ru.ru_utime;		//set start time
};