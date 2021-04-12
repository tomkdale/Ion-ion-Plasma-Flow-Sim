//-----------------------------------
//Function calculates new temperatures from old temperatures in each location on center grid
//TODO:FOR NOW PROGRAM DOES NOT CHANGE TEMPERATURE EVER
//-----------------------------------

#ifndef _TEMPERATURES_
#define _TEMPERATURES_

__global__ void calc_temperatures(double* newTempP, double* oldTempP, double* newTempN, double* oldTempN) {
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index

	if (i >= numSegR || j >= numSegZ) return;//out of bounds threads do nothing

	//just set new values to old
	newTempP[grid(i,j,numSegZ)] = oldTempP[grid(i,j,numSegZ)];
	newTempN[grid(i,j,numSegZ)] = oldTempN[grid(i,j,numSegZ)];

}
#endif
