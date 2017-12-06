//--------------------------------------------------
//Function calculates new voltages and new E fields in the R and Z direction for every point on edge grid for zero gradient case
//
//-------------------------------------------------

#ifndef _VOLTAGESZERO_
#define _VOLTAGESZERO_

//CALCULATES RED VALUES FOR ZERO GRADIENT CASE
__global__ void calc_voltagesRedZeroGradient(double*  newVolt, double* oldVolt, double* densityP, double* densityN, double dR, double dZ, bool* converged) {
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index

	i = (2 * i) + (j % 2); // edit i indexing term to switch from compressed red thread i j grid locations to actual grid locations

	int numEdgeR = numSegR + 1;
	int numEdgeZ = numSegZ + 1;
	if (i >= numEdgeR || j >= numEdgeZ) {
		return;
	}

	double densityP_avg, densityN_avg;
	double deltaDensity;

	if (i == 0 && j == 0) {//left inner corner
		newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded

	}
	else if (i == numEdgeR - 1 && j == 0) {//left outer corner
		if (i*dR <= thruster_Radius) {//if on thruster edge ground side
			newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded

		}
		else {//float values not on engine
			densityP_avg = densityP[grid(i - 1, j, numSegZ)];//density weighted average
			densityN_avg = densityN[grid(i - 1, j, numSegZ)];//density weighted average
			deltaDensity = densityP_avg - densityN_avg;

			newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
				(((dR *dR * q *deltaDensity) / epsilon0) +
				((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
					((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
					((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j + 1, numEdgeZ)] + //point right
						oldVolt[grid(i, j, numEdgeZ)]))); // point left

		}
	}
	else if (i == numEdgeR - 1 && j == numEdgeZ - 1) {//right outer corner
		densityP_avg = densityP[grid(i - 1, j - 1, numSegZ)];//density weighted average
		densityN_avg = densityN[grid(i - 1, j - 1, numSegZ)];//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j, numEdgeZ)] + //point right
					oldVolt[grid(i, j - 1, numEdgeZ)]))); // point left

	}
	else if (i == 0) {//centerline values
		newVolt[grid(i, j, numEdgeZ)] = oldVolt[grid(i + 1, j, numEdgeZ)];

	}

	else if (j == 0) {//left edge value
		if (i*dR <= thruster_Radius) {//if on thruster edge ground side
			newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded

		}
		else {//float values not on engine
			densityP_avg = (densityP[grid(i - 1, j, numSegZ)] + densityP[grid(i, j, numSegZ)]) / 2;//density weighted average
			densityN_avg = (densityN[grid(i - 1, j, numSegZ)] + densityN[grid(i, j, numSegZ)]) / 2;//density weighted average
			deltaDensity = densityP_avg - densityN_avg;

			newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
				(((dR *dR * q *deltaDensity) / epsilon0) +
				((1 + (1 / (2 * (double)i)))*oldVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
					((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
					((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j + 1, numEdgeZ)] + //point right
						oldVolt[grid(i, j, numEdgeZ)]))); // point left

		}
	}
	else if (i == numEdgeR - 1) {//outside edge value 
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i - 1, j, numSegZ)]) / 2;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i - 1, j, numSegZ)]) / 2;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j + 1, numEdgeZ)] + //point right
					oldVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}
	else if (j == numEdgeZ - 1) {//right edge value
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i, j - 1, numSegZ)]) / 2;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i, j - 1, numSegZ)]) / 2;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j, numEdgeZ)] + //point right
					oldVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}
	else {
		//if interior point
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i, j - 1, numSegZ)] + densityP[grid(i - 1, j, numSegZ)] + densityP[grid(i, j, numSegZ)]) / 4;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i, j - 1, numSegZ)] + densityN[grid(i - 1, j, numSegZ)] + densityN[grid(i, j, numSegZ)]) / 4;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* oldVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j + 1, numEdgeZ)] + //point right
					oldVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}

	if (fabs(newVolt[grid(i, j, numEdgeZ)] - oldVolt[grid(i, j, numEdgeZ)]) < maxVoltChange)
		converged[grid(i, j, numEdgeZ)] = true;
	else
		converged[grid(i, j, numEdgeZ)] = false;
}




//CALCULATES BLACK VALUES FOR ZERO GRADIENT CASE
__global__ void calc_voltagesBlackZeroGradient(double*  newVolt, double* oldVolt, double* densityP, double* densityN, double dR, double dZ, bool* converged) {
	int i = blockIdx.y * blockDim.y + threadIdx.y;//r position index
	int j = blockIdx.x * blockDim.x + threadIdx.x;//z position index

	i = (2 * i) + ((j + 1) % 2); // edit i indexing term to switch from compressed black thread i j grid locations to actual grid locations


	int numEdgeR = numSegR + 1;
	int numEdgeZ = numSegZ + 1;
	if (i >= numEdgeR || j >= numEdgeZ) {
		return;
	}

	double densityP_avg, densityN_avg;
	double deltaDensity;


	if (i == 0 && j == 0) {//left inner corner
		newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded
	}
	else if (i == numEdgeR - 1 && j == 0) {//left outer corner
		if (i*dR <= thruster_Radius) {//if on thruster edge ground side
			newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded
		}
		else {
			densityP_avg = densityP[grid(i - 1, j, numSegZ)];//density weighted average
			densityN_avg = densityN[grid(i - 1, j, numSegZ)];//density weighted average
			deltaDensity = densityP_avg - densityN_avg;

			newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
				(((dR *dR * q *deltaDensity) / epsilon0) +
				((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
					((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
					((dR * dR) / (dZ * dZ) * (newVolt[grid(i, j + 1, numEdgeZ)] + //point right
						oldVolt[grid(i, j, numEdgeZ)]))); // point left
		}
	}
	else if (i == numEdgeR - 1 && j == numEdgeZ - 1) {//right outer corner
		densityP_avg = densityP[grid(i - 1, j - 1, numSegZ)];//density weighted average
		densityN_avg = densityN[grid(i - 1, j - 1, numSegZ)];//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j, numEdgeZ)] + //point right
					newVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}
	else if (i == 0) {//if center value
		newVolt[grid(i, j, numEdgeZ)] = newVolt[grid(i + 1, j, numEdgeZ)];
	}

	else if (j == 0) {//left edge value
		if (i*dR <= thruster_Radius) {//if on thruster edge ground side
			newVolt[grid(i, j, numEdgeZ)] = 0.0;//Exit plate of the thruster is grounded
		}
		else {//float values not on engine
			densityP_avg = (densityP[grid(i - 1, j, numSegZ)] + densityP[grid(i, j, numSegZ)]) / 2;//density weighted average
			densityN_avg = (densityN[grid(i - 1, j, numSegZ)] + densityN[grid(i, j, numSegZ)]) / 2;//density weighted average
			deltaDensity = densityP_avg - densityN_avg;

			newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
				(((dR *dR * q *deltaDensity) / epsilon0) +
				((1 + (1 / (2 * (double)i)))*newVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
					((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
					((dR * dR) / (dZ * dZ) * (newVolt[grid(i, j + 1, numEdgeZ)] + //point right
						oldVolt[grid(i, j, numEdgeZ)]))); // point left
		}
	}
	else if (i == numEdgeR - 1) {//outside edge value
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i - 1, j, numSegZ)]) / 2;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i - 1, j, numSegZ)]) / 2;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*oldVolt[grid(i, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (newVolt[grid(i, j + 1, numEdgeZ)] + //point right
					newVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}
	else if (j == numEdgeZ - 1) {//right edge value
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i, j - 1, numSegZ)]) / 2;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i, j - 1, numSegZ)]) / 2;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*newVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (oldVolt[grid(i, j, numEdgeZ)] + //point right
					newVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}
	else {
		//if interior point
		densityP_avg = (densityP[grid(i - 1, j - 1, numSegZ)] + densityP[grid(i, j - 1, numSegZ)] + densityP[grid(i - 1, j, numSegZ)] + densityP[grid(i, j, numSegZ)]) / 4;//density weighted average
		densityN_avg = (densityN[grid(i - 1, j - 1, numSegZ)] + densityN[grid(i, j - 1, numSegZ)] + densityN[grid(i - 1, j, numSegZ)] + densityN[grid(i, j, numSegZ)]) / 4;//density weighted average
		deltaDensity = densityP_avg - densityN_avg;

		newVolt[grid(i, j, numEdgeZ)] = (1 - omega)* oldVolt[grid(i, j, numEdgeZ)] + omega *(1.0 / (2.0 * (1.0 + ((dR * dR) / (dZ * dZ))))) *
			(((dR *dR * q *deltaDensity) / epsilon0) +
			((1 + (1 / (2 * (double)i)))*newVolt[grid(i + 1, j, numEdgeZ)]) + //point outside
				((1 - (1 / (2 * (double)i)))* newVolt[grid(i - 1, j, numEdgeZ)]) + //point inside
				((dR * dR) / (dZ * dZ) * (newVolt[grid(i, j + 1, numEdgeZ)] + //point right
					newVolt[grid(i, j - 1, numEdgeZ)]))); // point left
	}

	if (fabs(newVolt[grid(i, j, numEdgeZ)] - oldVolt[grid(i, j, numEdgeZ)]) < maxVoltChange)
		converged[grid(i, j, numEdgeZ)] = true;
	else
		converged[grid(i, j, numEdgeZ)] = false;
}
#endif
