#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <time.h>
#include <algorithm>
#include <thread>
#include <numeric>
#include <chrono>
#include <iomanip>

using namespace std;

#include "parSolver.hpp"
#include "structUsed.hpp"

//========================================
// Flow field parameters
//========================================
double A = 1.0;                 // double gyre flow parameters
double mu = 1.0;
double eps = 0.0;
double s = 1.0;                 // scale for the workspace
double w = 2*M_PI*2;

//================================
// simulation parameters
//================================

// number of paths to simulate
int nPaths = 1000;

// integration step
double dt = 0.0001;

// noise intensity
double D = 1/30.0;

// gyre centre
double xC = -0.443780476044467;
double yC = -0.545824753793719;

// desired escape time
double Tesc_desired = 13.0; // the uncontrolled escape time for D=1/30 is about 5.7s

//================================
// phase space boundaries
//================================

//for general DG flow
double xmax = 0;
double xmin = -s;
double ymax = 0;
double ymin = -s;
double tmin = 0;
double tmax;

//===========================
// parameters for histogram
//===========================

// number of bins along each axis
int nXBins = 200;
int nYBins = 200;

// size of the complete histogram
double xSpBin = (xmax-xmin)/nXBins;
double ySpBin = (ymax-ymin)/nYBins;

// number of total histogram bins
int nHashBins = nXBins*nYBins;

//===========================
// hashbin function
//===========================
int getHashBinNumber(double x, double y){
    int xBin = (x-xmin)/xSpBin;
    int yBin = (y-ymin)/ySpBin;

    return (xBin + yBin*nXBins);
}

//=============================================
// flow field and diffusion functions
//=============================================
// no control
vector<double> DG_flow( double t, vector<double> x )
{
  vector<double> x_dot(2);

  x_dot[0] = (-M_PI*A*sin(M_PI*x[0]/s)*cos(M_PI*(x[1]/s)) - mu*x[0]);
  x_dot[1] = ( M_PI*A*cos(M_PI*x[0]/s)*sin(M_PI*(x[1]/s)) - mu*x[1]);

  return x_dot;
    
}

vector<double> g ( double t, vector<double> u )

// This should return a matrix, in general, but for additive noise,
// we can get by with return unit vector. 

{
  vector<double> uprime(2);

  uprime[0] = 1;
  uprime[1] = 1;

  return uprime; 
}

// with control
vector<double> DG_flow_control(double t, vector<double> x, double c)
{
      vector<double> x_dot = DG_flow(t,x);
      vector<double> F(2);
      F[0] = x_dot[0] + mu*x[0];
      F[1] = x_dot[1] + mu*x[1];

      double sqrtF1F2 = sqrt( F[0]*F[0] + F[1]*F[1] );

      x_dot[0] = x_dot[0] - c*F[1]/sqrtF1F2;
      x_dot[1] = x_dot[1] + c*F[0]/sqrtF1F2;
      
      return x_dot;

}

// iomanip function
string retTruncDbl(double val, int nDec){
    stringstream ss;
    ss << fixed << setprecision(nDec) << val;
    return ss.str();
}

//========================================
// function to read data files as vectors
//========================================

vector<double> readDataToVecs(string inFname){
    vector<double> dataMat;
    dataMat.clear();

    ifstream inF(inFname);
    cout << "Reading File " << inFname << endl;
    int q = 0;
    double val = 0;
    while( inF >> val ){
        q++;
        dataMat.push_back(val);
    }
    cout << "Done reading file " << inFname << endl;
    return dataMat;
}

//==========================================
// Main Code
//==========================================
int main(int argc, char *argv[]){

    // select integrator type (only EU_MAR with control)
    sdeIntegratorType sdeIntegrator(EU_MAR);
    
    stringstream tempStream;
    if (argc > 1){
        tempStream << argv[1] << " "; 
        tempStream >> nPaths;
    }

    time_t stTime, lvlStTime;
    stTime = time(NULL);

    // setup gsl random number generator
    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    // init point of simulation
    double x0 = xC;
    double y0 = yC;
    double t0 = 0;
    
    // number of threads supported
    unsigned int nThreads = thread::hardware_concurrency(); 
    int pathsPerSolver = nPaths/nThreads + 1;
    int totalPaths = pathsPerSolver*nThreads;

    // histogram data
    vector<int> histData(nHashBins,0);

    // avgEscapeTime
    double avgCtrlEffort;
    double avgEscapeTime;
   
    // progress counter variables
    int pIterPaths, cmpltPaths;
    double avgCmptTime;
    
    // output file for average escape times
    ofstream avgInfoOut;
    avgInfoOut.open("avgStatistics.txt", ofstream::out | ofstream::trunc);

    // compute paths for each noise level
    lvlStTime = time(NULL);
     
    // clear variables
    fill(histData.begin(), histData.end(), 0);
    avgEscapeTime = 0.0;
    avgCtrlEffort = 0.0;

    // set sigma value for this run
    double sigma = sqrt(2*D);
    
    cout << "computing paths for noise level D = 1/"<<( round(1/D) ) << endl; 
    cout << "number of threads used " << nThreads << endl;
    
    vector<parallelPathSolver> pathSolverVect;
    
    // setup solvers
    for( int i=0; i<nThreads; i++){
        pathSolverVect.push_back( parallelPathSolver(i,nHashBins,xmin,xmax,ymin,ymax,dt) );
    }
    
    // start solver threads on solvers
    vector<int> nPathCmplt(nThreads,0);
    for( int i=0; i<nThreads; i++){
        pathSolverVect[i].startThread(pathsPerSolver, x0, y0, t0, sdeIntegrator, sigma, xC, yC, r, DG_flow_control, g, getHashBinNumber, &nPathCmplt[i], Tesc_desired, s);
    }
    
    this_thread::sleep_for(chrono::milliseconds(100));
    pIterPaths = totalPaths;
    bool printStatus = false;
    while (true){
        cmpltPaths = 0;
        printStatus = false;
        for (int i=0; i<nPathCmplt.size(); i++)
          cmpltPaths += nPathCmplt[i];
          
        if(totalPaths>100 && totalPaths<20000){
            if (cmpltPaths%100==0)
                printStatus = true;
        }
        if(totalPaths>=20000){
            if (cmpltPaths%(totalPaths/100)==0)
                printStatus = true;
        }

        if ( printStatus ){
            if (pIterPaths!=cmpltPaths){  
                avgCmptTime = (double)(time(NULL)-lvlStTime)/cmpltPaths;
                cout << "\t" << retTruncDbl(cmpltPaths/(double)totalPaths*100,2) << "\% complete. ";
                cout << "Time remaining: " << round( (totalPaths-cmpltPaths)*avgCmptTime ) << "s."<< endl;
            }
            pIterPaths = cmpltPaths;
        }
        if (cmpltPaths==totalPaths)
            break;
    } 
    // join solvers
    for( int i=0; i<nThreads; i++){
        pathSolverVect[i].joinThread();
    }
    
    cout << "\tsynchronizing all threads...\n";
    
    // combine data
    for( int i=0; i<nThreads; i++){
        transform(histData.begin(), histData.end(), pathSolverVect[i].histData.begin(), histData.begin(), plus<int>() );    // add histogram data
        avgEscapeTime = accumulate(pathSolverVect[i].escapeTime.begin(),pathSolverVect[i].escapeTime.end(),avgEscapeTime);  // add upescape times
        avgCtrlEffort = accumulate(pathSolverVect[i].ctrlEffort.begin(),pathSolverVect[i].ctrlEffort.end(),avgCtrlEffort);  // add up control efforts
    }
    avgEscapeTime = avgEscapeTime/(totalPaths);
    avgCtrlEffort = avgCtrlEffort/(totalPaths);

    // save data for noise Level D[nLvlCnt]
    cout << "\tsaving data... " << endl;
    
    // save histogram data
    ofstream histDataOut;
    histDataOut.open("histData_nLevel_1_"+to_string( round(1/D) )+".txt", ofstream::out | ofstream::trunc);
    for (int i=0; i<histData.size(); i++)
        histDataOut << histData[i] << endl;

    histDataOut.close();
    
    //// save escape time data
    //ofstream escapeTimeOut;
    //escapeTimeOut.open("escapeTimes_nLevel_1_"+to_string( round(1/D) )+".txt", ofstream::out|ofstream::trunc);
    //for(int i=0; i<nThreads; i++)
    //    for(int j=0; j<pathSolverVect[i].escapeTime.size(); j++)
    //        escapeTimeOut << pathSolverVect[i].escapeTime[j] << endl;
    //
    //escapeTimeOut.close();
    
    //// save control effort data
    //ofstream ctrlEffOut;
    //ctrlEffOut.open("ctrlEffort_nLevel_1_"+to_string( round(1/D) )+".txt", ofstream::out|ofstream::trunc);
    //for(int i=0; i<nThreads; i++)
    //    for(int j=0; j<pathSolverVect[i].ctrlEffort.size(); j++)
    //        ctrlEffOut << pathSolverVect[i].ctrlEffort[j] << endl;

    //ctrlEffOut.close();

    //// save noise effort data
    //ofstream noiseEffOut;
    //noiseEffOut.open("noiseEffort_nLevel_1_"+to_string( round(1/D) )+".txt", ofstream::out|ofstream::trunc);
    //for(int i=0; i<nThreads; i++)
    //    for(int j=0; j<pathSolverVect[i].noiseEffort.size(); j++)
    //        noiseEffOut << pathSolverVect[i].noiseEffort[j] << endl;

    //noiseEffOut.close();

    // print average escape time for this noise level
    avgInfoOut << D << " " << avgEscapeTime << " " << avgCtrlEffort << endl;

    cout << "\tdata saving complete" << endl;
    cout << "\tStats : " << endl;
    cout << "\t\tNumber of paths computed        : " << totalPaths << endl;

    cout << "Total running time : " << ( time(NULL)-stTime ) << "s" << endl;
    avgInfoOut.close();
     
}


