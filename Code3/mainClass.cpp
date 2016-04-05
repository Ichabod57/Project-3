

#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>


using namespace std;

#include "Planet.h"
#include "SolarSystem.h"
#include "lib.h"
//#include "Radius.h"
//#include "Function.h"
//#include "Solvers.h"
void main() {

	
//initial conditions
double step;
cout <<"How big is the time step? If you don't tell me then I cannot do anything ";
cin >> step;
cin.get();

int iterations;
cout << "How many iterations ";
cin >> iterations;
cin.get();

int PlanetNum=10;
Planet Plan[10];
//SolarSystem thing;
double *Mass, *X, *Y, *Z, *VX, *VY, *VZ;


X = new double[PlanetNum];
Y = new double[PlanetNum];
Z = new double[PlanetNum];
VX = new double[PlanetNum];
VY = new double[PlanetNum];
VZ = new double[PlanetNum];
Mass= new double[PlanetNum];
double pi=3.14159265358979323;

/*
///////////// Whole Solar System//////////////
Planet Sun(1,0.00425,0,0,0,-0.002546,0);
//Planet Sun(1,0,0,0,0,0,0);
Planet Mercury(0.00000012,0.39,0,0,0,9.96,0);
Planet Venus(0.0000024,0.72,0,0,0,7.36,0);
Planet Earth(0.0000015,1,0,0,0,2*3.1415926535,0);
Planet Mars(0.00000033,1.52,0,0,0,5.06,0);
Planet Jupiter(0.00095,5.2,0,0,0,(2*pi/sqrt(5.2028)),0);
Planet Saturn(0.000275,9.54,0,0,0,2.04,0);
Planet Uranus(0.000044,19.19,0,0,0,1.43,0);
Planet Neptune(0.000051,30.06,0,0,0,1.14,0);
Planet Pluto(0.0000000056,39.53,0,0,0,0.99,0);
Plan[0]=Sun;
Plan[1]=Mercury;
Plan[2]=Venus;
Plan[3]=Earth;
Plan[4]=Mars;
Plan[5]=Jupiter;
Plan[6]=Saturn;
Plan[7]=Uranus;
Plan[8]=Neptune;
Plan[9]=Pluto;
*/
//////////Earth and Sun/////////////////

Planet Sun(1,0,0,0,0,0,0);
Planet Mercury(0.00000012,0.3075,0,0,0,12.44,0);
///Planet Jupiter(0.00095,5.2028,0,0,0,(2*pi/sqrt(5.2028)),0);
Plan[0]=Sun;
Plan[1]=Mercury;
//Plan[2]=Jupiter;




for (int j=0; j<PlanetNum; j++){
	Mass[j]=Plan[j].mas;
	X[j]=Plan[j].position[0];
	Y[j]=Plan[j].position[1];
	Z[j]=Plan[j].position[2];
	VX[j]=Plan[j].velocity[0];
	VY[j]=Plan[j].velocity[1];
	VZ[j]=Plan[j].velocity[2];}
	
//cout << Position[0][1] <<" "<<Position[1][1]<<" "<<Velocity[0][1]<<" "<<Velocity[1][1]<<" "<<Mass[1] <<endl;

Verlet(step,iterations,2,X,Y,Z,VX,VY,VZ,Mass);
//RungeKutta4(step,iterations,2,X,Y,Z,VX,VY,VZ,Mass);

cin.get();

}