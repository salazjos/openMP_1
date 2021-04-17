/* 
 * Joseph Salazar
 * salazjos@oregonstate.edu
 * openMPproject1.cpp
 * Base line instructor provided Monte Carlo simulation of laser beam reflection. Circle expands randomly to
 * different sizes and at different locations. Use OpenMP with 4 threads to compute the probabilty
 * of the laser being reflected off the circle and back on to the x plane given at 10 tries with
 * 1,000,000 trials each. 
 * 
 * */
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <stdio.h>
// setting the number of threads:
#ifndef NUMT
#define NUMT        4
#endif

// setting the number of trials in the monte carlo simulation:
#ifndef NUMTRIALS
#define NUMTRIALS    1000000
#endif

// how many tries to discover the maximum performance:
#ifndef NUMTRIES
#define NUMTRIES    10
#endif

// ranges for the random numbers:
const float XCMIN =    -1.0;
const float XCMAX =     1.0;
const float YCMIN =     0.0;
const float YCMAX =     2.0;
const float RMIN  =     0.5;
const float RMAX  =     2.0;

// function prototypes:
float        Ranf( float, float );
int          Ranf( int, int );
void         TimeOfDaySeed( );
float        calculateIntersection(const float &, const float &, const float &);
float        calculateInterSectionTime(const float &, const float &, const float &);
double       findShortestTime(double *);

// main program:
int main(int argc, const char * argv[]) {
    
    #ifndef _OPENMP
        fprintf( stderr, "No OpenMP support!\n" );
        return 1;
    #endif
    
    TimeOfDaySeed( );        // seed the random number generator
    
    omp_set_num_threads( NUMT );    // set the number of threads to use in the for-loop:`
	
	int num = omp_get_num_procs();
	printf("Cores = %d\n", num); //print the number of cores the running processor has.
    
    //define outside and fill outside of thread timing
    float *xcs = new float[NUMTRIALS];
    float *ycs = new float[NUMTRIALS];
    float *rs  = new float[NUMTRIALS];
    
    // fill the random-value arrays:
    for( int n = 0; n < NUMTRIALS; n++ )
    {
        xcs[n] = Ranf(XCMIN, XCMAX);
        ycs[n] = Ranf(YCMIN, YCMAX);
        rs[n]  = Ranf(RMIN, RMAX);
    }
    
    // get ready to record the maximum performance and the probability:
    float maxPerformance = 0.0;      // must be declared outside the NUMTRIES loop
    float currentProb;              // must be declared outside the NUMTRIES loop
    
    //array to hold the elapsed time
    double elapsedTimeArray[NUMTRIES];
    
    printf("Threads   Trials   Probability    MegaTrials/sec\n");  
    // looking for the maximum performance:
    for(int t = 0; t < NUMTRIES; t++)
    {
        double time0 = omp_get_wtime();
        int numHits = 0;
        
        #pragma omp parallel for default(none) shared(xcs,ycs,rs) reduction(+:numHits)
        for(int n = 0; n < NUMTRIALS; n++)
        {
            // randomize the location and radius of the circle:
            float xc = xcs[n];
            float yc = ycs[n];
            float r  =  rs[n];
            
            // solve for the intersection using the quadratic formula
            float tempA = 2.0;
            float tempB = -2.0 * ( xc + yc );
            float tempC = xc * xc + yc * yc - r * r;
            float tempD = tempB * tempB - 4.0 * tempA * tempC;
            
            if(tempD < 0){
                continue;
            }
            //else d > 0. hits the circle, can take square root
            // hits the circle:
            // get the first intersection:
            float t_min = calculateInterSectionTime(tempD, tempA, tempB);
            
            if(t_min < 0){ //circle engulfs laser pointer
                continue;
            }
            //else tmin > 0
            float tt = calculateIntersection(t_min, xc, yc);
            
            if(tt < 0){ //beam went up instead of down
                continue;
            }
            //else t > 0, beam hits the infinite plate
            numHits++;
        } //end inner for-loop openmp
        
        double time1 = omp_get_wtime();
        double megaTrialsPerSecond = (double)NUMTRIALS / ( time1 - time0 ) / 1000000.;
        if( megaTrialsPerSecond > maxPerformance )
            maxPerformance = megaTrialsPerSecond;
        currentProb = (float)numHits/(float)NUMTRIALS;
		
		// Print number of threads, number of trials, probability, megatrialsperseconds.
		printf("%d      %d     %2.4f      %8.2lf", NUMT, NUMTRIALS, currentProb, megaTrialsPerSecond);
		printf("\n");
   		 
		double elapsedTime = 1000000.*(time1-time0);
		elapsedTimeArray[t] = elapsedTime;
    
    }//end outer for-loop 
    
    //determine the shortest run time
    double shortestTime = findShortestTime(elapsedTimeArray);
	//print the shortest run time
    printf("Shortest run time = %10.2lf microseconds.\n", shortestTime);
    //print the maxPerformance value
    printf("Max performance = %8.2lf\n", maxPerformance);
    
	delete []xcs;
	delete []ycs;
	delete []rs;
	
    return 0;
}

double findShortestTime(double *elapsedTimeArr){
    double shortestTime = elapsedTimeArr[0];
    for(int i = 1; i < NUMTRIES; i++){
        if(elapsedTimeArr[i] < shortestTime)
            shortestTime = elapsedTimeArr[i];
    }
    return shortestTime;
}


//Instructor provided
float calculateInterSectionTime(const float &d, const float &a, const float &b){
    float temp_d = sqrt( d );
    float t1 = (-b + temp_d ) / ( 2.0*a );    // time to intersect the circle
    float t2 = (-b - temp_d ) / ( 2.0*a );    // time to intersect the circle
    float tmin = t1 < t2 ? t1 : t2;     // only care about the first intersection
    return tmin;
}

//Instructor provided
float Ranf( float low, float high ){
    float r = (float) rand();               // 0 - RAND_MAX
    float t = r  /  (float) RAND_MAX;       // 0. - 1.
    return   low  +  t * ( high - low );
}

//Instructor provided
int Ranf( int ilow, int ihigh ){
    float low = (float)ilow;
    float high = ceil( (float)ihigh );
    
    return (int) Ranf(low,high);
}

//Instructor provided
void TimeOfDaySeed( ){
    struct tm y2k = { 0 };
    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
    
    time_t  timer;
    time( &timer );
    double seconds = difftime( timer, mktime(&y2k) );
    unsigned int seed = (unsigned int)( 1000.*seconds );    // milliseconds
    srand( seed );
}

//Instructor provided
float calculateIntersection(const float &tmin, const float &xc, const float &yc){
    // where does it intersect the circle?
    float xcir = tmin;
    float ycir = tmin;
    
    // get the unitized normal vector at the point of intersection:
    float nx = xcir - xc;
    float ny = ycir - yc;
    float _n = sqrt( nx*nx + ny*ny );
    nx /= _n;    // unit vector
    ny /= _n;    // unit vector
    
    // get the unitized incoming vector:
    float inx = xcir - 0.;
    float iny = ycir - 0.;
    float in = sqrt( inx*inx + iny*iny );
    inx /= in;    
    iny /= in;    
    
    // get the outgoing (bounced) vector:
    float dot = inx*nx + iny*ny;
    float outx = inx - 2.*nx*dot;    
    float outy = iny - 2.*ny*dot;    
    
    // find out if it hits the infinite plate:
    float t = ( 0. - ycir ) / outy;
    return t;
}

