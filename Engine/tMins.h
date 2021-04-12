//----------------------------------
//Changes tMins[], spectral radiusR and spectral RadiusZ
//----------------------------------

#ifndef _tMINS_
#define _tMINS_
__global__ void calc_tmins(double *oldTempP, double *oldTempN, double *oldVelocityRP, double *oldVelocityRN, double *oldVelocityZP, double *oldVelocityZN, double* tMins,double* spectralRadiusR,double* spectralRadiusZ,double dR,double dZ) {
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index
	if (i >= numSegR || j >= numSegZ) return; //return if thread out of bounds

	double gamma = 1.3;//TODO, get actual value this is an estimate


	//calculate speed of sound for positive and negative charges
	double speedSoundP = sqrt(gamma * kboltz * oldTempP[grid(i, j, numSegZ)] / massP);
	double speedSoundN = sqrt(gamma * kboltz * oldTempN[grid(i, j, numSegZ)] / massN);

	//calculate maximum spectralRadiusR and spectralRadiusZ
	if(fabs(oldVelocityRP[grid(i, j, numSegZ)]) > fabs(oldVelocityRN[grid(i, j, numSegZ)]))
		spectralRadiusR[grid(i, j, numSegZ)] = fabs(oldVelocityRP[grid(i,j,numSegZ)]) + speedSoundP;
	else
		spectralRadiusR[grid(i, j, numSegZ)] = fabs(oldVelocityRN[grid(i, j, numSegZ)]) + speedSoundN;
	if (fabs(oldVelocityZP[grid(i, j, numSegZ)]) > fabs(oldVelocityZN[grid(i, j, numSegZ)]))
		spectralRadiusZ[grid(i, j, numSegZ)] = fabs(oldVelocityZP[grid(i, j, numSegZ)]) + speedSoundP;
	else
		spectralRadiusZ[grid(i, j, numSegZ)] = fabs(oldVelocityZN[grid(i, j, numSegZ)]) + speedSoundN;

	//calculate needed time for R and Z
	double tMinR = CFL * dR / spectralRadiusR[grid(i, j, numSegZ)];
	double tMinZ = CFL * dZ / spectralRadiusZ[grid(i, j, numSegZ)];

	//save the smallest tMin time
	if (tMinR < tMinZ)
		tMins[grid(i, j, numSegZ)] = tMinR;
	else
		tMins[grid(i, j, numSegZ)] = tMinZ;

	if (tMins[grid(i, j, numSegZ)] > 0.5 / plate_Frequency)
		tMins[grid(i, j, numSegZ)] = 0.5 / plate_Frequency;
}
#endif