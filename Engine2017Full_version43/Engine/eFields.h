//---------------------------------
//Function to calculate efields from surrounding voltages
//
//---------------------------------

#ifndef __EFIELDS_H
#define __EFIELDS_H

__global__ void calc_eFields(double* eFieldR, double* eFieldZ, double* volt, double dR, double dZ) {
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index
	

	if (i >= numSegR || j >= numSegZ) {
		return;
	}
	int numEdgesZ = numSegZ + 1;

	double vIn = (volt[grid(i, j, numEdgesZ)]+volt[grid(i,j+1,numEdgesZ)])/2;
	double vOut = (volt[grid(i + 1, j, numEdgesZ)] + volt[grid(i+1, j+1, numEdgesZ)]) / 2;
	double vLeft = (volt[grid(i+1, j, numEdgesZ)] + volt[grid(i , j, numEdgesZ)]) / 2;
	double vRight = (volt[grid(i + 1, j+1, numEdgesZ)] + volt[grid(i, j+1, numEdgesZ)]) / 2;

	//calculate eFields
	eFieldR[grid(i, j, numSegZ)] = -(vOut - vIn) / dR;
	eFieldZ[grid(i, j, numSegZ)] = -(vRight - vLeft) / dZ;
}

#endif // !__EFIELDS_H
