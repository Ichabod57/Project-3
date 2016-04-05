#include "Planet.h"
#include "SolarSystem.h"

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <iomanip>

#include <cmath>

using namespace std;
double Seperation(double x1, double x2){

return x1-x2;
}




double Radius(double X, double Y, double Z){
	return sqrt(X*X + Y*Y + Z*Z);}




double FunctionX(int PlanetNum, int planet, double *X, double *Y, double *Z, double *Mass){
	double pi=3.14159265358979323;
	double a=0;
   
	for (int i=0; i<PlanetNum; i++){
	
	if (i != planet){
	a += -(4*pi*pi*(Seperation(X[planet],X[i]))/pow(Radius(Seperation(X[planet],X[i]), Seperation(Y[planet],Y[i]), 0),3))*(Mass[i]/Mass[0])*(1+(0.000000459*0.000000459)/(63198*63198*Radius(X[1],Y[1],0)*Radius(X[1],Y[1],0)));
	//cout <<" "<<X[1]<<" "<<X[0]<<endl;
	}
	//cout << a << endl;
	}
return a;

}

double FunctionY(int PlanetNum, int planet, double *X, double *Y, double *Z, double *Mass){
	double pi=3.14159265358979323;
	double a=0;
    
	for (int i=0; i<PlanetNum; i++){
	if (i != planet){
	a += -(4*pi*pi*(Seperation(Y[planet],Y[i]))/pow(Radius(Seperation(X[planet],X[i]), Seperation(Y[planet],Y[i]), 0),3))*(Mass[i]/Mass[0])*(1+(0.000000459*0.000000459)/(63198*63198*Radius(X[1],Y[1],0)*Radius(X[1],Y[1],0)));
	//cout <<Y[1]<<" "<<Y[0]<<endl;
	}
	}
return a;

}

double Kinetic_Energy(int n, double *VX, double *VY, double *VZ, double *Mass){
	double a;
	a=0;
    for (int i=0; i<n; i++){
		a += Mass[i]*0.5*(pow(VX[i],2) + pow(VY[i],2) + pow(VZ[i],2));

	}
return a;
}

double Potential_Energy(int n, double *X, double *Y, double *Z, double *Mass){
	double a;
	double Pi=3.14159265358979323;
	a=0;

	for (int i=0; i<n; i++){
		for (int j=i+1; j<n-1; j++){
		a += -((4*Pi*Pi)/Radius(Seperation(X[i],X[j]),Seperation(Y[i],Y[j]),Seperation(Z[i],Z[j])))*(Mass[j]*Mass[i]);
		}
	}
	return a;

}


double Angular_Momentum(int PlanetNum, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ, double *Mass){
double momentum=0;

for (int i=0; i<PlanetNum; i++){
	momentum = momentum+Mass[i]*(X[i]*VY[i] - Y[i]*VX[i]);
}

return momentum;
}






void Verlet(double step, int numstep, int PlanetNum, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ, double *Mass){
int v=0;
int vv=0;
double time2=0;
double *X2, *Y2, *Z2, *VX2, *VY2, *VZ2;
double a, b, c, d;
double inEnergy, finalEnergy, time, inMomentum, finalMomentum;
double year=0;
X2 = new double[PlanetNum];
Y2 = new double[PlanetNum];
Z2 = new double[PlanetNum];
VX2 = new double[PlanetNum];
VY2 = new double[PlanetNum];
VZ2 = new double[PlanetNum];

inEnergy = Kinetic_Energy(PlanetNum,VX,VY,VZ,Mass)+Potential_Energy(PlanetNum,X,Y,Z,Mass);
inMomentum = Angular_Momentum(PlanetNum,X,Y,Z,VX,VY,VZ,Mass);
//ofstream VerletTest;
ofstream VerletData;

//VerletTest.open("VT_Mercury.txt");
VerletData.open("VD_Mercury.txt");


VerletData<<0<<" ";
	for (int j=0; j<PlanetNum; j++){
		VerletData<<X[j]<<" "<<Y[j]<<" ";
	}
	VerletData<<endl;
	


for (int i=0; i < numstep; i++){
	
	//verlet
for (int j=0; j<PlanetNum; j++){
 b= FunctionX(PlanetNum,j,X,Y,Z,Mass);
  d=FunctionY(PlanetNum,j,X,Y,Z,Mass); 
 X2[j] = X[j] + step*VX[j]+0.5*step*step*b;
 Y2[j] = Y[j] + step*VY[j]+0.5*step*step*d;
}

//set position of the sun to 0
X2[0]=0;
Y2[0]=0;

for (int k=0; k<PlanetNum; k++){
 b= FunctionX(PlanetNum,k,X,Y,Z,Mass);
  d=FunctionY(PlanetNum,k,X,Y,Z,Mass);
 a=FunctionX(PlanetNum,k,X2,Y2,Z2,Mass);
 c=FunctionY(PlanetNum,k,X2,Y2,Z2,Mass);
 VX2[k] = VX[k] + 0.5*step*(a + b);
 VY2[k] = VY[k] + 0.5*step*(c + d);
}
//set velocity of the sun to 0
VX2[0]=0;
VY2[0]=0;

for (int L=0; L<PlanetNum; L++){
 X[L]=X2[L];
 Y[L]=Y2[L];
 VX[L]=VX2[L];
 VY[L]=VY2[L];
}
//cout <<Radius(X[1],Y[1],0)<<endl;
finalEnergy=Kinetic_Energy(PlanetNum,VX,VY,VZ,Mass)+Potential_Energy(PlanetNum,X,Y,Z,Mass);
finalMomentum=Angular_Momentum(PlanetNum,X,Y,Z,VX,VY,VZ,Mass);
//output the change inenergy and the positions of planets

v++;
time=step*(i+1);
time2=step + time2;
//cout << time2 << endl;
////////puts energy difference, time and earth radius into a text file/////////////////
while (time2>0.2407){
cout << time2 << endl; 
cout << Y[1]/X[1] << endl;
time2=0;
}



if (v==1/step){
//if (v==20){
	//VerletTest<<time <<" "<<fabs(inEnergy-finalEnergy)<<" "<<fabs(inMomentum-finalMomentum)<<" "<<Radius(X[1],Y[1],0)<<endl;
	
	
	VerletData<<time<<" ";
	for (int j=0; j<PlanetNum; j++){
		VerletData<<X[j]<<" "<<Y[j]<<" ";
	}
	VerletData<<endl;
	

	//v=0;
}
	
}

}


void RungeKutta4(double step, int numstep, int PlanetNum, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ, double *Mass){

	double *kx1, *ky1, *kx2, *ky2, *kx3, *ky3, *kx4, *ky4,*kvx1, *kvy1, *kvx2, *kvy2, *kvx3, *kvy3, *kvx4, *kvy4;
	int v=0;
	double inEnergy, finalEnergy, time, inMomentum, finalMomentum;
kx1=new double[PlanetNum]; ky1=new double[PlanetNum]; kx2=new double[PlanetNum];
ky2=new double[PlanetNum]; kx3=new double[PlanetNum]; ky3=new double[PlanetNum];
kx4=new double[PlanetNum]; ky4=new double[PlanetNum]; kvx1=new double[PlanetNum];
kvy1=new double[PlanetNum]; kvx2=new double[PlanetNum]; kvy2=new double[PlanetNum];
kvx3=new double[PlanetNum]; kvy3=new double[PlanetNum]; kvx4=new double[PlanetNum];
kvy4=new double[PlanetNum];

//////////////////////////////////calculate initial energy
inEnergy = Kinetic_Energy(PlanetNum,VX,VY,VZ,Mass)+Potential_Energy(PlanetNum,X,Y,Z,Mass);
inMomentum = Angular_Momentum(PlanetNum, X,Y,Z,VX,VY,VZ,Mass);
ofstream KuttaTest;
//ofstream KuttaData;
KuttaTest.open("HD_Planet.txt");
//KuttaData.open("RKD_AP_S_500_0.001.txt");

/*
KuttaData<<0<<" ";
	for (int j=0; j<PlanetNum; j++){
		KuttaData<<0<<" "<<X[j]<<" "<<Y[j]<<" ";
	}
	KuttaData<<endl;
	*/
////////////////////////Run RungeKutta//////////////////
for (int i=0; i<numstep; i++){
	//find K1's
for (int j=0; j<PlanetNum; j++){
	kvx1[j] = FunctionX(PlanetNum, j, X, Y, Z, Mass);
	kvy1[j] = FunctionY(PlanetNum, j, X, Y, Z, Mass);
	kx1[j] = VX[j];
	ky1[j] = VY[j];
}

//find K2's
for (int j=0; j<PlanetNum; j++){
	X[j]=X[j]+step*0.5*kx1[j];
	Y[j]=Y[j]+step*0.5*ky1[j];
}

for (int j=0; j<PlanetNum; j++){
	kvx2[j] = FunctionX(PlanetNum, j, X, Y, Z, Mass);
	kvy2[j] = FunctionY(PlanetNum, j, X, Y, Z, Mass);
	kx2[j] = VX[j]+step*0.5*kvx1[j];
	ky2[j] = VY[j] + step*0.5*kvy1[j];
}

//find K3's
for (int j=0; j<PlanetNum; j++){
	X[j]=X[j]-step*0.5*kx1[j];
	Y[j]=Y[j]-step*0.5*ky1[j];
	X[j]=X[j]+step*0.5*kx2[j];
	Y[j]=Y[j]+step*0.5*ky2[j];
}

for (int j=0; j<PlanetNum; j++){
	kvx3[j] = FunctionX(PlanetNum, j, X, Y, Z, Mass);
	kvy3[j] = FunctionY(PlanetNum, j, X, Y, Z, Mass);
	kx3[j] = VX[j] + step*0.5*kvx2[j];
	ky3[j] = VY[j] + step*0.5*kvy2[j];
}

//find k4
for (int j=0; j<PlanetNum; j++){
    X[j]=X[j]-step*0.5*kx2[j];
	Y[j]=Y[j]-step*0.5*ky2[j];
	X[j]=X[j]+step*kx3[j];
	Y[j]=Y[j]+step*ky3[j];
}

for (int j=0; j<PlanetNum; j++){
	kvx4[j] = FunctionX(PlanetNum, j, X, Y, Z, Mass);
	kvy4[j] = FunctionY(PlanetNum, j, X, Y, Z, Mass);
	kx4[j] = VX[j] + step*kvx3[j];
	ky4[j] = VY[j] + step*kvy3[j];
}

//find the new values of x, y, vx and vy
for (int j=0; j<PlanetNum; j++){
	X[j]=X[j]-step*kx3[j];
	Y[j]=Y[j]-step*ky3[j];
}

for (int j=0; j<PlanetNum; j++){
	X[j] = X[j] + (step/6)*(kx1[j] + 2*kx2[j] + 2*kx3[j] + kx4[j]);
	Y[j] = Y[j] + (step/6)*(ky1[j] + 2*ky2[j] + 2*ky3[j] + ky4[j]);
    VX[j] = VX[j] + (step/6)*(kvx1[j] + 2*kvx2[j] + 2*kvx3[j] + kvx4[j]);
    VY[j] = VY[j] + (step/6)*(kvy1[j] + 2*kvy2[j] + 2*kvy3[j] + kvy4[j]);
}

//set the suns velocity and coordinates to zero
X[0]=0; Y[0]=0; VX[0]=0; VY[0]=0;

//output the energy
finalEnergy=Kinetic_Energy(PlanetNum,VX,VY,VZ,Mass)+Potential_Energy(PlanetNum,X,Y,Z,Mass);
finalMomentum = Angular_Momentum(PlanetNum,X,Y,Z,VX,VY,VZ,Mass);


v++;
time=step*(i+1);
//////////puts energy difference, and earth radius in a text file/////////////
if (v==1/step){

//if (v==20){
	KuttaTest<<time<<" "<<fabs(inEnergy-finalEnergy)<<" "<<fabs(inMomentum-finalMomentum)<<" "<<Radius(X[1],Y[1],0)<<endl;
	/*
	//cout<<Radius(X[1],Y[1],Z[1])<<endl;
	KuttaData<<time<<" ";
	for (int j=0; j<PlanetNum; j++){
		KuttaData<<X[j]<<" "<<Y[j]<<" ";
	}
	KuttaData<<endl;
	*/
v=0;}


//cout << Radius(X[1],Y[1],0)<<endl;
}

}

double Precession(int PlanetNum, double *X, double *Y, double *Z, double *VX, double *VY, double *VZ, double *Mass){
	double a;
	
	a = 1 + (0.000000459*0.000000459)/(Radius(X[1],Y[1],0)*Radius(X[1],Y[1],0)*63198*63198);
	//cout << a << endl;
	//cout <<Radius(X[1],Y[1],0)<<endl;
	return a;
}