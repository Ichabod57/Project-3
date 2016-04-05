//#pragma once

//#ifndef SOLARSYSTEM_H
//#define SOLARSYSTEM_H

#include <iostream>
#include <cmath>
#include <fstream>
#include "Planet.h"



//class SolarSystem
////{
//public:

double Seperation(double, double);
double Radius(double , double , double );
double FunctionX(int, int,  double *, double *, double *, double *);
double FunctionY(int, int, double *, double *, double *, double *);
double Kinetic_Energy(int, double *,double *, double *, double *);
double Potential_Energy(int, double *,double *,double *, double *);
double AngulerMomentum(int, double *, double *, double *, double *, double *, double *, double *);
void Verlet(double, int, int, double *, double *, double *, double *, double *, double *, double *);
void RungeKutta4(double, int, int, double *, double *, double *, double *, double *, double *, double *);
double Precession(int, double *, double *, double*, double *, double *, double *, double *);

//#endif // SOLARSYSTEM_H