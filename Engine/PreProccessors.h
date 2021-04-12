#ifndef _PREPROCESSORS_H_
#define _PREPROCESSORS_H_


//-------------------------------
//Preprocessors for Engine code.
//
//-------------------------------

#define PI 3.141592653589793238462643
#define MU0 4*PI*1e-7 //T.m/A
#define kboltz 1.38064852e-23 //Boltzmann's constant in J/K
#define epsilon0 8.854187817e-12//C^2/N.m^2 
#define q 1.60217662e-19 //elementary charge
#define AMU 1.6726219e-27 //one AMU in kilograms
#define massP (126.9 * AMU)
#define massN (126.9 * AMU)
#define grid(r,z,nz) ((r)*(nz))+(z) //used index a 1D array as if it had 2 dimensions r and z
#define rCell(i,dR) (((i) + .5) * (dR)) //used to find radius of a center cell
#define threadsPerBlock2D 16
#define omega 1.99 //Red-Black iterative solver coefficient
#define maxVoltChange 1e-7
#define CFL 0.2 //time step stability
#define fresh 1 //start program from ignition conditions or some other values in space
#define thrusterTemp 1160.4 // temperature of positive and negative ions directly after exiting t
#define minDensity 1.0e14 //density in space prior to thruster initialization
#define grounded 1 //1 for right end grounded case, 0 for zero gradient



#define Total_Time 10.0e-6 //s
#define numSegR 100
#define numFluxR (numSegR+1)
#define	numSegZ 200
#define numFluxZ (numSegZ+1)
#define	length_R 0.1 //m
#define	length_Z 0.2 //m
#define	thruster_Radius 0.05 //m
#define	plate_Frequency 1.0e6 //Frequency (Hz) of voltage reversals of the accelerating grid
#define plate_Voltage 500.0 //Voltage (V) of the accelerating grid
#define thrusterDensity 1.0e15 //Density (#/m^3) of the exhaust from the thruster grid
#define thrusterVelocity sqrt(2.0 * q * plate_Voltage / massP)

#endif