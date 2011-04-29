#include <iostream> 
#include <sys/resource.h>

using namespace std;

namespace LPG{

#ifndef TIMEKEEPER_HPP
#define TIMEKEEPER_HPP
	
	class TimeKeeper{		
	  private:
		struct timeval startTime;	//time clock starts
		struct timeval tempTime;	//general time variable
		struct rusage ru;		//to query system
	  public:
		TimeKeeper();
		
		double getElapsedTimeSec(); 	//get the time elapsed from startTime til now
                void reset();                   //set the timer
	};


	
#endif
};
