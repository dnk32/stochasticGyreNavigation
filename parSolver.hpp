//parSolver.hpp
#ifndef PAR_SOLVER_HPP
#define PAR_SOLVER_HPP

#include <vector>
#include <thread>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "structUsed.hpp"
using namespace std;
typedef vector<double> dVec;


//=========================================
// Paraller Path Solver Class 
//=========================================

class parallelPathSolver{
    private:
        int solverID;
    public:
    // x and y trajectories
    dVec X, Y, T;

    // histogram data
    vector<int> histData;

    dVec escapeTime;
    dVec ctrlEffort;
    dVec noiseEffort;

    // limits of workspace
    double xmin, xmax, ymin, ymax;

    // time step
    double dt;

    // placeholder variables
    double t, xn, yn, escapeT, vX, vY;
    int pInd;
    int hashBin, hashBin_3D;
    bool escapeRight;
    double sigma;
    double cntrlCost, noiseCost;

    // Counters for number of paths and iterations
    int pathCnt, iterCnt;

    // thread for solving
    thread th;

    parallelPathSolver(int nSolver, int nBins, double xL, double xH, double yL, double yH, double timeStep) : 
        solverID(nSolver), 
        histData(nBins,0), 
        xmin(xL), xmax(xH), ymin(yL), ymax(yH),
        dt(timeStep)
    { }

    void runPaths (int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double, dVec, double), dVec g(double, dVec), int hBinFnc(double, double), int *nPathsCmplt, double T_des, double gyre_width);

    void startThread(int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double, dVec, double), dVec g(double, dVec), int hBinFnc(double, double), int *nPathsCmplt, double T_des, double gyre_width);
    
    void joinThread();
    double computeDist2Boundary(double, double);
    double computeControlParam(double, double, double, double);
};

#endif
