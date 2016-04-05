#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstdlib>
using namespace std;

#include "time.h"
#include "lib.h"
#include "Radius.h"
#include "Function.h"
#include "Solvers.h"
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

double xin, yin, xvin, yvin;

cout << "initial x";
cin >> xin;
cin.get();

cout << "initial y";
cin >> yin;
cin.get();

cout << "initial x velocity ";
cin >> xvin;
cin.get();

cout << "initial y velocity ";
cin >> yvin;
cin.get();

//Euler(step,iterations,xin,yin,xvin,yvin);
Verlet(step,iterations,xin,yin,xvin,yvin);
//RunK4(step,iterations,xin,yin,xvin,yvin);
cin.get();

}