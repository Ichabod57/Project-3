//#pragma once
#ifndef PLANET_H
#define PLANET_H
class Planet
{
public:
	
double position[3];
double velocity[3];
double mas;
Planet(double mass, double x, double y, double z, double vx, double vy, double vz);
Planet();




};




#endif