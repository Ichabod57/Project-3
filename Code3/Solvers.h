void Euler(double step, int numstep, double inx, double iny, double invx, double invy){

	double x_earth, vx_earth, y_earth, vy_earth, x2_earth, y2_earth, vx2_earth, vy2_earth, radius;
	double pi =3.14159265358979323;
	x_earth = inx;
	y_earth = iny;
	vx_earth = invx;
	vy_earth = invy;
	
	
	
	for (int i=0; i<numstep; i++){
    radius = sqrt(x_earth*x_earth+y_earth*y_earth);
	// radius = 1;
    x2_earth=x_earth+vx_earth*step;
	y2_earth=y_earth+vy_earth*step;
	vx2_earth=vx_earth-step*(4*pi*pi*x_earth/pow(radius,3));
	vy2_earth=vy_earth-step*(4*pi*pi*y_earth/pow(radius,3));

	x_earth=x2_earth;
	y_earth=y2_earth;
	vx_earth=vx2_earth;
	vy_earth=vy2_earth;
	cout << x_earth << " " << y_earth<< " "<< radius <<endl;
}

}

void Verlet(double step, int numstep, double inx, double iny, double invx, double invy){


double x_earth, vx_earth, y_earth, vy_earth, x2_earth, y2_earth, vx2_earth, vy2_earth, inEnergy, FinalEnergy;
int k=0;
	double pi =3.14159265358979323;
	invy=pi*2;
inEnergy = 0.5*invx*invx + 0.5*invy*invy-(4*pi*pi/Rad(inx,iny));
	x_earth = inx;
	y_earth = iny;
	vx_earth = invx;
	vy_earth = invy;
	
	ofstream Verlet;
	Verlet.open("Verlet_En_SE_NoClass.txt");
	
	for (int i=0; i<numstep; i++){
   k++;
	

    x2_earth=x_earth+vx_earth*step+(step*step/2)*Fun1(x_earth,Rad(x_earth,y_earth));
    y2_earth=y_earth+vy_earth*step+(step*step/2)*Fun1(y_earth,Rad(x_earth,y_earth));
	
	vx2_earth=vx_earth+(step/2)*(Fun1(x2_earth,Rad(x2_earth,y2_earth))+Fun1(x_earth,Rad(x_earth,y_earth)));
	vy2_earth=vy_earth+(step/2)*(Fun1(y2_earth,Rad(x2_earth,y2_earth))+Fun1(y_earth,Rad(x_earth,y_earth)));
	

	x_earth=x2_earth;
	y_earth=y2_earth;
	vx_earth=vx2_earth;
	vy_earth=vy2_earth;
	//cout << x_earth << " " << y_earth<< " "<< Rad(x_earth,y_earth) <<endl;

	FinalEnergy = 0.5*vx_earth*vx_earth + 0.5*vy_earth*vy_earth -(4*pi*pi/Rad(x_earth,y_earth));


	if (k==100){
		Verlet<<fabs(FinalEnergy-inEnergy)<<endl;
		k=0;
	}

}
}

void RunK4(double step, int numstep, double inx, double iny, double invx, double invy){

double x_earth, vx_earth, y_earth, vy_earth, inEnergy, FinalEnergy;
double k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, k1vx, k1vy, k2vx, k2vy, k3vx, k3vy, k4vx, k4vy;
double pi=3.14159265358979323826264;
int k=0;
invy=pi*2;
inEnergy = 0.5*invx*invx + 0.5*invy*invy-(4*pi*pi/Rad(inx,iny));
	x_earth = 1;
	y_earth = 0;
	vx_earth = 0;
	vy_earth = pi*2;
	
	//make text file
	ofstream RK;
	RK.open("RK_En_SE_NoClass.txt");
	
	//Runge Kutta Solve
	for (int i=0; i<numstep; i++){
     k++;


	 k1vx =Fun1(x_earth,Rad(x_earth,y_earth));
     k1vy =Fun1(y_earth,Rad(x_earth,y_earth));
	 k1x = vx_earth;
	 k1y =vy_earth;

	  k2vx =Fun1(x_earth+(0.5*step*k1x),Rad(x_earth+(0.5*step*k1x),y_earth+(0.5*step*k1y)));
     k2vy =Fun1(y_earth+(0.5*step*k1y),Rad(x_earth+(0.5*step*k1x),y_earth+(0.5*step*k1y)));
	 k2x =(vx_earth+(0.5*step*k1vx));
	 k2y =(vy_earth+(0.5*step*k1vy));
	 
	 k3vx =Fun1(x_earth+(0.5*step*k2x),Rad(x_earth+(0.5*step*k2x),y_earth+(0.5*step*k2y)));
     k3vy =Fun1(y_earth+(0.5*step*k2y),Rad(x_earth+(0.5*step*k2x),y_earth+(0.5*step*k2y)));
	 k3x =(vx_earth+(0.5*step*k2vx));
	 k3y =(vy_earth+(0.5*step*k2vy));

	k4vx =Fun1(x_earth+step*k3x,Rad(x_earth+step*k3x,y_earth+step*k3y));
     k4vy =Fun1(y_earth+step*k3y,Rad(x_earth+step*k3x,y_earth+step*k3y));
	 k4x =(vx_earth+step*k3vx);
	 k4y =(vy_earth+step*k3vy);
	
	
//cout << k1vx << " " <<k2vx<<" "<<k3vx<<" "<<k4vx<<endl;

	 x_earth = x_earth + (step/6)*(k1x + (2*k2x) + (2*k3x) + k4x);
	 y_earth = y_earth + (step/6)*(k1y + (2*k2y) + (2*k3y) + k4y);
	 vx_earth = vx_earth + (step/6)*(k1vx + (2*k2vx) + (2*k3vx) + k4vx);
	 vy_earth = vy_earth + (step/6)*(k1vy + (2*k2vy) + (2*k3vy) + k4vy);
	 

	 FinalEnergy = 0.5*vx_earth*vx_earth + 0.5*vy_earth*vy_earth -(4*pi*pi/Rad(x_earth,y_earth));
	 cout << x_earth <<" "<<y_earth<<" "<<Rad(x_earth,y_earth) << endl;
if (k==100){
		RK<<fabs(FinalEnergy-inEnergy)<<endl;
		k=0;
	}


	 //cout <<inEnergy <<" "<<FinalEnergy << endl;


	}


}