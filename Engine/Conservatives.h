//---------------------------------------------
//calculates conservative equations for zero gradient edge values
//---------------------------------------------

#ifndef _CONSERVATIVES_
#define _CONSERVATIVES_

__global__ void calc_conservatives(double* fluxR1, double* fluxR2, double* fluxR3, double* fluxR4, double* fluxR5, double* fluxR6, double* fluxZ1, double* fluxZ2, double* fluxZ3, double* fluxZ4, double* fluxZ5, double* fluxZ6, double* tempP, double* tempN, double* newDensityP, double* oldDensityP, double* newDensityN, double* oldDensityN, double* newVelocityRP, double* oldVelocityRP, double* newVelocityRN, double* oldVelocityRN, double* newVelocityZP, double* oldVelocityZP, double* newVelocityZN, double* oldVelocityZN, double* eFieldR, double* eFieldZ, double* collisionFreq, double dR, double dZ, double dT) {

	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index

	if (i >= numSegR || j >= numSegZ) return;//out of bounds threads do nothing

	//1 = Conservation of mass (positives)
	//2 = Conservation of mass (negatives)
	//3 = Conservation of momentum in R direction (positives)
	//4 = Conservation of momentum in R direction (negatives)
	//5 = Conservation of momentum in Z direction (positives)
	//6 = Conservation of momentum in Z direction (negatives)


	//1
	newDensityP[grid(i, j, numSegZ)] =  oldDensityP[grid(i, j, numSegZ)]
										- dT*((fluxR1[grid(i+1 ,j , numSegZ)] - fluxR1[grid(i,j,numSegZ)]) / dR	+ (fluxZ1[grid(i, j + 1, numFluxZ)] - fluxZ1[grid(i, j, numFluxZ)]) / dZ)
											/ rCell(i, dR);
	//2
	newDensityN[grid(i, j, numSegZ)] =  oldDensityN[grid(i, j, numSegZ)]
										- dT*((fluxR2[grid(i + 1, j, numSegZ)] - fluxR2[grid(i, j, numSegZ)]) / dR + (fluxZ2[grid(i, j + 1, numFluxZ)] - fluxZ2[grid(i, j, numFluxZ)]) / dZ)
											/ rCell(i, dR);
	//3
	newVelocityRP[grid(i, j, numSegZ)] = oldDensityP[grid(i, j, numSegZ)]*oldVelocityRP[grid(i, j, numSegZ)] / newDensityP[grid(i, j, numSegZ)]
											- dT*((fluxR3[grid(i + 1, j, numSegZ)] - fluxR3[grid(i, j, numSegZ)]) / dR + (fluxZ3[grid(i, j + 1, numFluxZ)] - fluxZ3[grid(i, j, numFluxZ)]) / dZ
												- (oldDensityP[grid(i, j, numSegZ)]*kboltz*tempP[grid(i, j, numSegZ)] / massP)
												- (rCell(i,dR)*oldDensityP[grid(i, j, numSegZ)]
													*(q*eFieldR[grid(i, j, numSegZ)] / massP
														- collisionFreq[grid(i, j, numSegZ)]*oldVelocityRP[grid(i, j, numSegZ)])))
												/(rCell(i, dR)*newDensityP[grid(i, j, numSegZ)]);
	//4
	newVelocityRN[grid(i, j, numSegZ)] = oldDensityN[grid(i, j, numSegZ)] * oldVelocityRN[grid(i, j, numSegZ)] / newDensityN[grid(i, j, numSegZ)]
											- dT*((fluxR4[grid(i + 1, j, numSegZ)] - fluxR4[grid(i, j, numSegZ)]) / dR	+ (fluxZ4[grid(i, j + 1, numFluxZ)] - fluxZ4[grid(i, j, numFluxZ)]) / dZ
												- (oldDensityN[grid(i, j, numSegZ)]*kboltz*tempN[grid(i, j, numSegZ)] / massN)
												- (rCell(i, dR)*oldDensityN[grid(i, j, numSegZ)]
													*(-q*eFieldR[grid(i, j, numSegZ)] / massN
														+ collisionFreq[grid(i, j, numSegZ)] * oldVelocityRN[grid(i, j, numSegZ)])))
											/ (rCell(i, dR)*newDensityN[grid(i, j, numSegZ)]);
	//5
	newVelocityZP[grid(i, j, numSegZ)] = oldDensityP[grid(i, j, numSegZ)] * oldVelocityZP[grid(i, j, numSegZ)] / newDensityP[grid(i, j, numSegZ)]
											- dT*((fluxR5[grid(i + 1, j, numSegZ)] - fluxR5[grid(i, j, numSegZ)]) / dR + (fluxZ5[grid(i, j + 1, numFluxZ)] - fluxZ5[grid(i, j, numFluxZ)]) / dZ
												- (rCell(i, dR)*oldDensityP[grid(i, j, numSegZ)]
													*(q*eFieldZ[grid(i, j, numSegZ)] / massP
														- collisionFreq[grid(i, j, numSegZ)] * oldVelocityZP[grid(i, j, numSegZ)])))
											/ (rCell(i, dR)* newDensityP[grid(i, j, numSegZ)]);
	//6
	newVelocityZN[grid(i, j, numSegZ)] = oldDensityN[grid(i, j, numSegZ)] * oldVelocityZN[grid(i, j, numSegZ)] / newDensityN[grid(i, j, numSegZ)]
											- dT*((fluxR6[grid(i + 1, j, numSegZ)] - fluxR6[grid(i, j, numSegZ)]) / dR + (fluxZ6[grid(i, j + 1, numFluxZ)] - fluxZ6[grid(i, j, numFluxZ)]) / dZ
												- (rCell(i, dR)*oldDensityN[grid(i, j, numSegZ)]
													*(-q*eFieldZ[grid(i, j, numSegZ)] / massN
														+ collisionFreq[grid(i, j, numSegZ)] * oldVelocityZN[grid(i, j, numSegZ)])))
											/ (rCell(i, dR)* newDensityN[grid(i, j, numSegZ)]);
		
}//end of function
#endif