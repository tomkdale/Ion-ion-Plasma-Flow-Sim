//---------------------------------------------
//calculates fluxes to be used by conservatives.h
//---------------------------------------------

#ifndef _FLUXES_
#define _FLUXES_


__global__ void calc_fluxesR(double* tempP, double* tempN, double* densityP,  double* densityN,  double* velocityRP, double* velocityRN,  double* velocityZP, double* velocityZN, double* eFieldR, double* eFieldZ, double* fluxR1, double* fluxR2, double* fluxR3, double* fluxR4, double* fluxR5, double* fluxR6, double* spectralRadiusR, double dR, double dZ)
{
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index
	if (i >= numFluxR || j >= numSegZ) return;//out of bounds threads do nothing 

	if (i == 0){//if center line value, the only thing should be numerical diffusion correction
		//While densities and velocities are the same across the centerline, the value of 'r' flips sign from -dR/2 to +dR/2. So, the numerical diffusion becomes ... 
		fluxR1[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityP[grid(i, j, numSegZ)];
		fluxR2[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityN[grid(i, j, numSegZ)];
		fluxR3[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)];
		fluxR4[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)];
		fluxR5[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)];
		fluxR6[grid(i, j, numSegZ)] = -.5 * spectralRadiusR[grid(i, j, numSegZ)] * dR * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)];

		//fluxR1[grid(i, j, numSegZ)] = 0.; //not correct, centerline values need numerical diffusion as shown above
		//fluxR2[grid(i, j, numSegZ)] = 0.;
		//fluxR3[grid(i, j, numSegZ)] = 0.;
		//fluxR4[grid(i, j, numSegZ)] = 0.;
		//fluxR5[grid(i, j, numSegZ)] = 0.;
		//fluxR6[grid(i, j, numSegZ)] = 0.;

	}
	else if (i == numFluxR - 1) {//if outer edge value, assume that there is no radial gradient across the boundary
		fluxR1[grid(i, j, numSegZ)] = (rCell(i, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] //above ... the 'i' value = 'i-1' value
										+ rCell(i - 1, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]) / 2.0;//below

		fluxR2[grid(i, j, numSegZ)] = (rCell(i, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] //above
										+ rCell(i - 1, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)])/ 2.0;//below

		fluxR3[grid(i, j, numSegZ)] = (rCell(i, dR)*(densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]
										+ densityP[grid(i - 1, j, numSegZ)] * kboltz * tempP[grid(i - 1, j, numSegZ)] / massP) //above
									+ rCell(i - 1, dR)*(densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]
										+ densityP[grid(i - 1, j, numSegZ)] * kboltz * tempP[grid(i - 1, j, numSegZ)]/massP))/ 2.0; //below

		fluxR4[grid(i, j, numSegZ)] = (rCell(i, dR)*(densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]
										+ densityN[grid(i - 1, j, numSegZ)] * kboltz * tempN[grid(i - 1, j, numSegZ)] / massP)//above
									+ rCell(i - 1, dR)*(densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]
											+ densityN[grid(i - 1, j, numSegZ)] * kboltz * tempN[grid(i - 1, j, numSegZ)] / massP)) / 2.0;//below

		fluxR5[grid(i, j, numSegZ)] = (rCell(i, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)] 
										+ rCell(i - 1, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)]) / 2.0;
		fluxR6[grid(i, j, numSegZ)] = (rCell(i, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]
										+ rCell(i - 1, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]) / 2.0;

		//Numerical diffusion
		fluxR1[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityP[grid(i - 1, j, numSegZ)] - rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)]);
		fluxR2[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityN[grid(i - 1, j, numSegZ)] - rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)]);
		fluxR3[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]
																						- rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]);
		fluxR4[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]
																						- rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]);
		fluxR5[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)]
																						- rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)]);
		fluxR6[grid(i, j, numSegZ)] -= .5 * spectralRadiusR[grid(i - 1, j, numSegZ)] * (rCell(i, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]
																						- rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]);
	}
	else {// if normal interior point
		fluxR1[grid(i, j, numSegZ)] = (rCell(i, dR)*densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] 
									+ rCell(i - 1, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]) / 2.0;

		fluxR2[grid(i, j, numSegZ)] = (rCell(i, dR)*densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] 
									+ rCell(i - 1, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]) / 2.0;

		fluxR3[grid(i, j, numSegZ)] = (rCell(i, dR)*(densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)]
										+ densityP[grid(i, j, numSegZ)] * kboltz * tempP[grid(i, j, numSegZ)] / massP) 
									+ rCell(i - 1, dR)*(densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] 
										+ densityP[grid(i - 1, j, numSegZ)] * kboltz * tempP[grid(i - 1, j, numSegZ)] / massP)) / 2.0;

		fluxR4[grid(i, j, numSegZ)] = (rCell(i, dR)*(densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)]
										+ densityN[grid(i, j, numSegZ)] * kboltz * tempN[grid(i, j, numSegZ)] / massN) 
									+ rCell(i - 1, dR)*(densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] 
										+ densityN[grid(i - 1, j, numSegZ)] * kboltz * tempN[grid(i - 1, j, numSegZ)] / massN)) / 2.0;

		fluxR5[grid(i, j, numSegZ)] = (rCell(i, dR)*densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] 
									+ rCell(i - 1, dR)*densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)]) / 2.0;

		fluxR6[grid(i, j, numSegZ)] = (rCell(i, dR)*densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] 
									+ rCell(i - 1, dR)*densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]) / 2.0;

		//Numerical diffusion
		fluxR1[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] - rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)]);
		fluxR2[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] - rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)]);
		fluxR3[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)]
												- rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityRP[grid(i - 1, j, numSegZ)]);
		fluxR4[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)]
												- rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityRN[grid(i - 1, j, numSegZ)]);
		fluxR5[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
												- rCell(i - 1, dR) * densityP[grid(i - 1, j, numSegZ)] * velocityZP[grid(i - 1, j, numSegZ)]);
		fluxR6[grid(i, j, numSegZ)] -= .5 * .5*(spectralRadiusR[grid(i - 1, j, numSegZ)] + spectralRadiusR[grid(i, j, numSegZ)]) 
											* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
												- rCell(i - 1, dR) * densityN[grid(i - 1, j, numSegZ)] * velocityZN[grid(i - 1, j, numSegZ)]);
	}
}//end of calc_fluxesR function

__global__ void calc_fluxesZ(double* tempP, double* tempN, double* densityP, double* densityN, double* velocityRP, double* velocityRN, double* velocityZP, double* velocityZN, double* eFieldR, double* eFieldZ,  double* fluxZ1, double* fluxZ2, double* fluxZ3, double* fluxZ4, double* fluxZ5, double* fluxZ6, double* spectralRadiusZ, double dR, double dZ, double currentTime)
{
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index
	if (i >= numSegR || j >= numFluxZ) return;//out of bounds threads do nothing

	if (j == 0) {//for left boundary
				 //While densities and velocities are the same across the centerline, the value of 'r' flips sign from -dR/2 to +dR/2. So, the numerical diffusion becomes ... 
		if (i*dR <= thruster_Radius) {//IF LESS THAN THRUSTER RADIUS
				if (fmod(currentTime, (1. / plate_Frequency)) >= (.5 / plate_Frequency)) {//inject positives
					fluxZ1[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] + rCell(i, dR) * thrusterDensity * thrusterVelocity) / 2.0;

					fluxZ2[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ3[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ4[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ5[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
														+ densityP[grid(i, j, numSegZ)] * kboltz*tempP[grid(i, j, numSegZ)] / massP)
												+ rCell(i, dR)*((thrusterDensity *thrusterVelocity * thrusterVelocity)
														+ (thrusterDensity * kboltz * thrusterTemp / massP))) / 2.0;

					fluxZ6[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
																	+ densityN[grid(i, j, numSegZ)] * kboltz*tempN[grid(i, j, numSegZ)] / massN)
													+ rCell(i, dR)*(0.0 + (minDensity * kboltz * thrusterTemp / massN))) / 2.0;

					//Numerical Diffusion
					fluxZ1[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityP[grid(i, j, numSegZ)] - rCell(i, dR) * thrusterDensity);
					fluxZ2[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityN[grid(i, j, numSegZ)] - rCell(i, dR) * minDensity);
					fluxZ3[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] - 0.0);
					fluxZ4[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] - 0.0);
					fluxZ5[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)])	* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
																									- rCell(i, dR) * thrusterDensity * thrusterVelocity);
					fluxZ6[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] - 0.0);
				}
				else {//inject negatives
					fluxZ1[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ2[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] + rCell(i, dR) * thrusterDensity * thrusterVelocity) / 2.0;

					fluxZ3[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ4[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] + 0.0) / 2.0;

					fluxZ5[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
																	+ densityP[grid(i, j, numSegZ)] * kboltz*tempP[grid(i, j, numSegZ)] / massP)
													+ rCell(i, dR)*(0.0 + (minDensity * kboltz * thrusterTemp / massP))) / 2.0;

					fluxZ6[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
																	+ densityN[grid(i, j, numSegZ)] * kboltz*tempN[grid(i, j, numSegZ)] / massN)
													+ rCell(i, dR)*((thrusterDensity * thrusterVelocity * thrusterVelocity)
																	+ (thrusterDensity * kboltz * thrusterTemp / massN))) / 2.0;
					
					//Numerical Diffusion
					fluxZ1[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityP[grid(i, j, numSegZ)] - rCell(i, dR) * minDensity);
					fluxZ2[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityN[grid(i, j, numSegZ)] - rCell(i, dR) * thrusterDensity);
					fluxZ3[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] - 0.0);
					fluxZ4[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)]) * (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] - 0.0);
					fluxZ5[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)])	* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] - 0.0);
					fluxZ6[grid(i, j, numFluxZ)] -= .5 * (spectralRadiusZ[grid(i, j, numSegZ)])	* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
																									- rCell(i, dR) * thrusterDensity * thrusterVelocity);
				}		
		}//end of thruster radius
		else {//Outside thruster radius at j=0
			fluxZ1[grid(i, j, numFluxZ)] = rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)];
			fluxZ2[grid(i, j, numFluxZ)] = rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)];
			fluxZ3[grid(i, j, numFluxZ)] = rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)];
			fluxZ4[grid(i, j, numFluxZ)] = rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)];
			fluxZ5[grid(i, j, numFluxZ)] = rCell(i, dR) * (densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
														+ densityP[grid(i, j, numSegZ)] * kboltz*tempP[grid(i, j, numSegZ)] / massP);
			fluxZ6[grid(i, j, numFluxZ)] = rCell(i, dR) * (densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
														+ densityN[grid(i, j, numSegZ)] * kboltz*tempN[grid(i, j, numSegZ)] / massN);

		//Numerical Diffusion
			fluxZ1[grid(i, j, numFluxZ)] -= 0.0;
			fluxZ2[grid(i, j, numFluxZ)] -= 0.0;
			fluxZ3[grid(i, j, numFluxZ)] -= 0.0;
			fluxZ4[grid(i, j, numFluxZ)] -= 0.0;
			fluxZ5[grid(i, j, numFluxZ)] -= 0.0;
			fluxZ6[grid(i, j, numFluxZ)] -= 0.0;
		}	
	}//end of if j==0
	else if (j == numFluxZ - 1) {//for right boundary: with zero gradient, the outside point 'j' has the same value as the inside point 'j-1'
		fluxZ1[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)] 
										+ rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]) / 2.0;
		fluxZ2[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)] 
										+ rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ3[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityRP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]
										+ rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityRP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ4[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityRN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]
										+ rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityRN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ5[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]
										+ densityP[grid(i, j - 1, numSegZ)] * kboltz*tempP[grid(i, j - 1, numSegZ)] / massP)
										+ rCell(i, dR) * (densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]
										+ densityP[grid(i, j - 1, numSegZ)] * kboltz*tempP[grid(i, j - 1, numSegZ)] / massP)) / 2.0;

		fluxZ6[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]
										+ densityN[grid(i, j - 1, numSegZ)] * kboltz*tempN[grid(i, j - 1, numSegZ)] / massN)
										+ rCell(i, dR) * (densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]
										+ densityN[grid(i, j - 1, numSegZ)] * kboltz*tempN[grid(i, j - 1, numSegZ)] / massN)) / 2.0;

		//Numerical Diffusion
		fluxZ1[grid(i, j, numFluxZ)] -= 0.0;
		fluxZ2[grid(i, j, numFluxZ)] -= 0.0;
		fluxZ3[grid(i, j, numFluxZ)] -= 0.0;
		fluxZ4[grid(i, j, numFluxZ)] -= 0.0;
		fluxZ5[grid(i, j, numFluxZ)] -= 0.0;
		fluxZ6[grid(i, j, numFluxZ)] -= 0.0;
	}//end of j==numFluxZ-1
	else {//for interior point
		fluxZ1[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] + rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]) / 2.0;
		fluxZ2[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] + rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ3[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
										+ rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityRP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ4[grid(i, j, numFluxZ)] = (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
										+ rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityRN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]) / 2.0;

		fluxZ5[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)] 
											+ densityP[grid(i, j, numSegZ)] * kboltz*tempP[grid(i, j, numSegZ)] / massP) 
										+ rCell(i, dR) * (densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)] 
											+ densityP[grid(i, j - 1, numSegZ)] * kboltz*tempP[grid(i, j - 1, numSegZ)] / massP)) / 2.0;

		fluxZ6[grid(i, j, numFluxZ)] = (rCell(i, dR) * (densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
											+ densityN[grid(i, j, numSegZ)] * kboltz*tempN[grid(i, j, numSegZ)] / massN)
										+ rCell(i, dR) * (densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]
											+ densityN[grid(i, j - 1, numSegZ)] * kboltz*tempN[grid(i, j - 1, numSegZ)] / massN)) / 2.0;


		//Numerical Diffusion
		fluxZ1[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)]) 
										* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] - rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)]);
		fluxZ2[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)])
										* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] - rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)]);

		fluxZ3[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)])
										* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityRP[grid(i, j, numSegZ)]
										- rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityRP[grid(i, j - 1, numSegZ)]);

		fluxZ4[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)])
										* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityRN[grid(i, j, numSegZ)]
										- rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityRN[grid(i, j - 1, numSegZ)]);
		fluxZ5[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)])
										* (rCell(i, dR) * densityP[grid(i, j, numSegZ)] * velocityZP[grid(i, j, numSegZ)]
										- rCell(i, dR) * densityP[grid(i, j - 1, numSegZ)] * velocityZP[grid(i, j - 1, numSegZ)]);
		fluxZ6[grid(i, j, numFluxZ)] -= .5 * .5 * (spectralRadiusZ[grid(i, j, numSegZ)] + spectralRadiusZ[grid(i, j - 1, numSegZ)])
										* (rCell(i, dR) * densityN[grid(i, j, numSegZ)] * velocityZN[grid(i, j, numSegZ)]
										- rCell(i, dR) * densityN[grid(i, j - 1, numSegZ)] * velocityZN[grid(i, j - 1, numSegZ)]);
	}//end of interior points
}//end of calc_fluxesZ function

#endif