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

void parallelPathSolver::runPaths (int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double , dVec, double), dVec g(double, dVec), int getHashBinNumber(double, double), int *nPathsCmplt, double T_des, double gyre_width, double dR_A, double dR_B){
    pathCnt = 0; iterCnt = 0;
    escapeTime.clear();
    ctrlEffort.clear();
    noiseEffort.clear();

    ofstream pathSave; bool savedOnePath = false;
    ofstream tCurrOut, dist2BndryOut, contParamOut;
    if (nPaths ==1 ){
        tCurrOut.open("tCurr"+to_string(solverID)+".txt", ofstream::out|ofstream::trunc);
        dist2BndryOut.open("dist2Bnry"+to_string(solverID)+".txt", ofstream::out|ofstream::trunc);
        contParamOut.open("contParam"+to_string(solverID)+".txt", ofstream::out|ofstream::trunc);
        pathSave.open("Path"+to_string(solverID)+".txt", ofstream::out | ofstream::trunc );
    }
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

        // tracking if threshold has been crossed
        bool threshCrossed = false; 
        double thresholdTime = 0.0;
        double lambda_thresh = 0.6;
        vector<double> contPar_all;
        contPar_all.clear();
        contPar_all.push_back(0.0);
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
            X.push_back(xn);
            Y.push_back(yn);
            T.push_back(t);
            pInd++;

            dist2Boundary = computeDist2Boundary(xn,yn);
            if (!threshCrossed && dist2Boundary < gyre_width/2*lambda_thresh){
                threshCrossed = true;
                thresholdTime = t;
            }
            if ( dist2Boundary > gyre_width/2*lambda_thresh)
                threshCrossed = false;
            controlParam = computeControlParam(dist2Boundary, gyre_width, T_des, t, dR_A, dR_B,thresholdTime, lambda_thresh);
            contPar_all.push_back(controlParam);
            
            cntrlCost += 0.5*controlParam*controlParam*dt; 
            noiseCost += 0.5*( eta_x*eta_x + eta_y*eta_y )*dt;
            
            if( nPaths == 1){
                tCurrOut << t << endl;
                dist2BndryOut << dist2Boundary << endl;
                contParamOut << controlParam << endl;
                pathSave << xn << " " << yn << " " << t << endl;
            }
        }
        // update hitogram if escape happens from right boundary
        if (escapeRight){
            if (!savedOnePath && t < 10.0){
                dist2BndryOut.open("dist2Bnry"+to_string(solverID)+".txt", ofstream::out|ofstream::trunc);
                contParamOut.open("contParam"+to_string(solverID)+".txt", ofstream::out|ofstream::trunc);
                pathSave.open("Path"+to_string(solverID)+".txt", ofstream::out | ofstream::trunc );
            }
            for (int k=pInd; k>=0; k--){
            //for (int k=pInd; k>=pInd*0.9; k--){
                hashBin = getHashBinNumber(X[k], Y[k]);
                histData[hashBin] = histData[hashBin] + 1;
                if (!savedOnePath && t < 10.0){
                    pathSave << X[k] << " " << Y[k] << " " << T[k] << endl;
                    dist2Boundary = computeDist2Boundary(X[k],Y[k]);
                    contParamOut << contPar_all[k] << endl;
                    dist2BndryOut << dist2Boundary << endl;
                }
            }
            // log time of escape
            escapeTime.push_back(t);
            ctrlEffort.push_back(cntrlCost);
            noiseEffort.push_back(noiseCost);
            pathCnt++;
            *nPathsCmplt = pathCnt;
            if (!savedOnePath && t < 10.0){
                pathSave.close();
                dist2BndryOut.close();
                contParamOut.close();
                savedOnePath = true;
            }
        }
        iterCnt++;
        if(pathCnt==nPaths)
            break;
    }// end of running nPaths for noise Level D[nLvlCnt]

    tCurrOut.close();
    dist2BndryOut.close();
    contParamOut.close();
}// end function runpaths

void parallelPathSolver::startThread(int nPaths, double x0, double y0, double t0, sdeIntegratorType sdeIntegrator, double sigma, double xC, double yC, gsl_rng *r, dVec f(double, dVec, double ), dVec g(double, dVec), int hBinFnc(double, double), int *nPathsCmplt, double T_des, double gyre_width, double dR_A, double dR_B){
    cout << "\tSolver started on thread " << solverID << endl;
    th = thread(&parallelPathSolver::runPaths, this, nPaths, x0, y0, t0, sdeIntegrator, sigma, xC,  yC, r, f, g, hBinFnc, nPathsCmplt, T_des, gyre_width, dR_A, dR_B);
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
            
double parallelPathSolver::computeControlParam(double dist2Boundary, double gyre_width, double t_des, double t_current, double dR_A, double dR_B, double t_thresh, double lambda_thresh){

    double s = gyre_width/2;
    //double lambda_s = 0.2;
    //double lambda_t = 2;  // (1-lambda) is the fraction of time the final transit takes to go from  lambda_s*s to boundary
    double dR_ratio = 400;
    if ( dist2Boundary > lambda_s*s && t_current < (1-lambda_t)*t_des )
        return 0.0;
    if ( dist2Boundary < lambda_s*s && t_current > (1-lambda_t)*t_des )
        return c_max;//0.0;
    if ( dist2Boundary >= lambda_s*s && t_current >= (1-lambda_t)*t_des ){
        //if( t_des < 30)   
            return c_max;
        //else
        //    return 0.0;
    }
    else{ // dist < g_w/8 && t_cur < lambda*t_des
        double predicted_t_4_esc = t_current/( 1 - (lambda_t)*dist2Boundary/(lambda_s*s) );
        //double predicted_t_4_esc = ( lambda_thresh*t_current - lambda_s*t_thresh)/(lambda_thresh - lambda_s);
        if (predicted_t_4_esc > t_des)
            return 0.0;
        double dR = D*log( t_des/predicted_t_4_esc )*dR_ratio;
        double contPar =  ( dR_B - sqrt( dR_B*dR_B + 4*dR_A*dR ) )/(2*dR_A);
        if (isnan(contPar) )
            double dummyParam = 1;
        return (contPar<c_min)?c_min:contPar;
        //return contPar;
    }
}

//==================
//Old code
//==================
//double parallelPathSolver::computeControlParam(double dist2Boundary, double gyre_width, double t_des, double t_current, double dR_A, double dR_B){
//    double t_left = t_des - t_current;   //what if time_left is negative?
//    double predicted_t_4_esc = dist2Boundary*2/gyre_width*5.7;
//    //double dR = (t_left>0)?10*D*log( t_left/predicted_t_4_esc ):c_max; // if t_left is bigger dR >0;
//    double dR = ((t_left>0)?1:-1)*10*D*log( fabs(t_left)/predicted_t_4_esc ); // if t_left is bigger dR >0;
//    double contPar = 0.0;
//    double temp = ( (dR_B*dR_B + 4*dR_A*dR) < 0 )?( dR_B/(2*dR_A) ):(dR_B - sqrt( dR_B*dR_B + 4*dR_A*dR ) )/(2*dR_A);
//    
//    if( dist2Boundary > ( gyre_width*3/4.0 ) ){ // if the particle is close to the center
//        if( t_left > predicted_t_4_esc ){
//            contPar = temp;
//        }
//        else if( t_left < predicted_t_4_esc ){
//            contPar = temp; 
//        }
//    }
//    else if( dist2Boundary < ( gyre_width*3/4.0 ) ){
//        if( t_left > predicted_t_4_esc ){
//            contPar = temp;
//        }
//        else if( t_left < predicted_t_4_esc*0.1 ){
//            contPar = temp;
//        }
//    }
//    //if (t_left < 0)
//    //    contPar = c_max;
//    return contPar;
//}
