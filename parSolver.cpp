// parSolver.cpp
#include <iostream>
#include <vector>
#include <thread>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <fstream>

#include "parSolver.hpp"

using namespace std;

//=========================================
// Paraller Path Solver Clasis functions 
//=========================================

void parallelPathSolver::runPaths (int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double , dVec, double), dVec g(double, dVec), int getHashBinNumber(double, double), int *nPathsCmplt, double T_des, double gyre_width){
    pathCnt = 0; iterCnt = 0;
    escapeTime.clear();
    ctrlEffort.clear();
    noiseEffort.clear();

    while(true){
        
        // initialize path at start coordinate and time
        X.clear();
        Y.clear();
        T.clear();
        X.push_back(x0);
        Y.push_back(y0);
        T.push_back(t0);
        
        t = t0;
       
        // initialize path coordinanate dounter and escape indicator
        pInd = 0;
        escapeRight = false;
        cntrlCost = 0.0;
        noiseCost = 0.0;

        // Spectral energy paramater 
        double eta_x, eta_y;
    
        // placeholder vectors for flow velocities and current coordinates
        vector<double> vFlow;
        vector<double> xt(2);

        double controlParam = 0.0;
        double dist2Boundary = 0.0;

        while (true){
            // update position
            // current 
            xt[0] = X[pInd];
            xt[1] = Y[pInd];
            vFlow = f(t, xt, controlParam);
            
            // noise components
            eta_x = gsl_ran_gaussian(r,1.0);
            eta_y = gsl_ran_gaussian(r,1.0);
            
            xn = X[pInd] + vFlow[0]*dt + sigma*sqrt(dt)*eta_x;           
            yn = Y[pInd] + vFlow[1]*dt + sigma*sqrt(dt)*eta_y;
            
            // check if path escapes the gyre
            if (xn>xmax || xn < xmin || yn > ymax || yn < ymin){
                // check if the escape happens from the right boundary
                //if (xn>xmax)
                    escapeRight = true;
    
                break;
            }
                
            // update time
            t += dt;
            cntrlCost += 0.0; 
            noiseCost += 0.5*( eta_x*eta_x + eta_y*eta_y )*dt;
            X.push_back(xn);
            Y.push_back(yn);
            T.push_back(t);
            pInd++;

            dist2Boundary = computeDist2Boundary(xn,yn);
            controlParam = computeControlParam(dist2Boundary, gyre_width, T_des, t);
        }
        // update hitogram if escape happens from right boundary
        if (escapeRight){
            //ofstream pathSave;
            //pathSave.open("Solver"+to_string(solverID)+"_Path"+to_string(pathCnt)+".txt", ofstream::out | ofstream::trunc );
            //for (int k=pInd; k>=0; k--){
            for (int k=pInd; k>=pInd*0.9; k--){
                hashBin = getHashBinNumber(X[k], Y[k]);
                histData[hashBin] = histData[hashBin] + 1;
                //pathSave << X[k] << " " << Y[k] << endl;
            }
            // log time of escape
            escapeTime.push_back(t);
            ctrlEffort.push_back(cntrlCost);
            noiseEffort.push_back(noiseCost);
            pathCnt++;
            *nPathsCmplt = pathCnt;
            //pathSave.close();
        }
        iterCnt++;
        if(pathCnt==nPaths)
            break;
    }// end of running nPaths for noise Level D[nLvlCnt]

}// end function runpaths

void parallelPathSolver::startThread(int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double, dVec, double ), dVec g(double, dVec), int hBinFnc(double, double), int *nPathsCmplt, double T_des, double gyre_width){
    cout << "\tSolver started on thread " << solverID << endl;
    th = thread(&parallelPathSolver::runPaths, this, nPaths, x0, y0, t0, sdeIntegrator, sigma, xC,  yC, r, f, g, hBinFnc, nPathsCmplt, T_des, gyre_width);
    //th.join();
}

void parallelPathSolver::joinThread(){
    th.join();
    cout << "\tSolver completed on thread " << solverID << endl;
}

double parallelPathSolver::computeDist2Boundary(double x, double y){
    double minDist = fabs(x-xmin);
    double temp = fabs(x-xmax);
    minDist = ( minDist > temp )? temp : minDist;
    temp = fabs(y-ymin);
    minDist = ( minDist > temp )? temp : minDist;
    temp = fabs(y-ymax);
    minDist = ( minDist > temp )? temp : minDist;
    return minDist;
}
            
double parallelPathSolver::computeControlParam(double dist2Boundary, double gyre_width, double t_des, double t_current){
    double t_left = t_des - t_current;   //what if time_left is negative?
    double contPar = 0.0;
    if( dist2Boundary > ( gyre_width*3/4.0 ) ){ // if the particle is close to the center
        if( t_left < ( dist2Boundary*2/gyre_width*t_des ) ){
           contPar = 0.1;
        }
    }
    else if( dist2Boundary < ( gyre_width/4.0 ) ){
        if( t_left > ( dist2Boundary*2/gyre_width*t_des ) ){
            contPar = -0.5;
        }
        else if( t_left < ( dist2Boundary*2/gyre_width*t_des )*0.1 ){
            contPar = 0.1;
        }
    }
    return contPar;
}
