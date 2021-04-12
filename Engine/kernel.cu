//-----------------------------------------------------
//Engine code analysing ion-ion flow
//Written by Dr. Kamesh Sankaran, Thomas Dale and Travis Widmer
//Summer 2017
//-----------------------------------------------------

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <time.h>
#include "PreProccessors.h"
#include "tMins.h"
#include "VoltagesGrounded.h"
#include "VoltagesZeroGradient.h"
#include "Conservatives.h"
#include "Temperatures.h"
#include "Fluxes.h"
#include "eFields.h"
#include <random>


int main() {
	//-----------------------------------------
	clock_t start = clock(); // real run time
	double t = 0.0; //current Engine time
	//-----------------------------------------
	//Do necessary calculations for code setup
	int numCenters = numSegR *numSegZ;
	int numEdges = (numSegR + 1)*(numSegZ + 1);
	int numRFluxes = (numSegR + 1) * numSegZ;
	int numZFluxes = numSegR * (numSegZ + 1);
	double dR = length_R / numSegR;
	double dZ = length_Z / numSegZ;
	size_t sizeCenters = numCenters * sizeof(double);
	size_t sizeEdges = numEdges * sizeof(double);
	size_t sizeBoolEdges = numEdges * sizeof(bool);
	size_t sizeFluxR = numRFluxes * sizeof(double);
	size_t sizeFluxZ = numZFluxes * sizeof(double);
	dim3 blocksPerGridCenters((numSegZ / threadsPerBlock2D) + 1, (numSegR / threadsPerBlock2D) + 1);
	dim3 blocksPerGridEdges(((numSegZ +1) / threadsPerBlock2D) + 1, ((numSegR+1)/ threadsPerBlock2D) + 1);
	dim3 blocksPerGridFluxRs((numSegZ / threadsPerBlock2D) + 1, ((numFluxR) / threadsPerBlock2D) + 1);
	dim3 blocksPerGridFluxZs(((numFluxZ) / threadsPerBlock2D) + 1, (numSegR / threadsPerBlock2D) + 1);
	dim3 threadsPerBlockDim3(threadsPerBlock2D, threadsPerBlock2D);
	dim3 blocksPerGridCompressed((numSegZ / threadsPerBlock2D) + 1, ((numSegR / 2) / threadsPerBlock2D) + 1);
	//-----------------------------------------

	//-----------------------------------------
	//create array variables and allocate memory
	double *tempP, *tempN, *volt, *densityP, *densityN, *velocityRP, *velocityRN, *velocityZP, *velocityZN, *tMins, *eFieldR, *eFieldZ, *spectralRadiusR, *spectralRadiusZ, *collisionFreq;
	bool *redBlackConvergence;
	tempP = (double*)malloc(sizeCenters);
	tempN = (double*)malloc(sizeCenters);
	volt = (double*)malloc(sizeEdges);
	densityN = (double*)malloc(sizeCenters);
	densityP = (double*)malloc(sizeCenters);
	velocityRP = (double*)malloc(sizeCenters);
	velocityRN = (double*)malloc(sizeCenters);
	velocityZP = (double*)malloc(sizeCenters);
	velocityZN = (double*)malloc(sizeCenters);
	tMins = (double*)malloc(sizeCenters);
	eFieldR = (double*)malloc(sizeCenters);
	eFieldZ = (double*)malloc(sizeCenters);
	spectralRadiusR = (double*)malloc(sizeCenters);
	spectralRadiusZ = (double*)malloc(sizeCenters);
	collisionFreq = (double*)malloc(sizeCenters);
	redBlackConvergence = (bool*)malloc(sizeBoolEdges);


	if (fresh == 1) {// Create fresh start initial variables
		for (int it = 0; it < numEdges; it++)
			volt[it] = 0.;
		for (int it = 0; it < numCenters; it++)
			tempP[it] = thrusterTemp;
		for (int it = 0; it < numCenters; it++)
			tempN[it] = thrusterTemp;
		for (int it = 0; it < numCenters; it++)
			densityP[it] = minDensity;
		for (int it = 0; it < numCenters; it++)
			densityN[it] = minDensity;
		for (int it = 0; it < numCenters; it++)//velocity of positive ions in r direction
			velocityRP[it] = 0;
		for (int it = 0; it < numCenters; it++)// velocity of negative in r direction
			velocityRN[it] = 0;
		for (int it = 0; it < numCenters; it++)// velocity of positive ions in z direction
			velocityZP[it] = 0;
		for (int it = 0; it < numCenters; it++)// velocity of negative ions in z direction
			velocityZN[it] = 0;
	}
	else {//assign initial variables from initialValues.txt
		FILE* arrays;
		arrays = fopen("initialValues.txt", "r");
		fscanf(arrays, "%*s");
		for (int it = 0; it < numEdges; it++)
			fscanf(arrays, "%le", &volt[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &tempP[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &tempN[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &densityP[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &densityN[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &velocityRP[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &velocityRN[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &velocityZP[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &velocityZN[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &eFieldR[it]);
		fscanf(arrays, "%*s");
		for (int it = 0; it < numCenters; it++)
			fscanf(arrays, "%le", &eFieldZ[it]);
		fscanf(arrays, "%*s");
		fscanf(arrays, "%le", &t);
		fclose(arrays);
	}

	//initialize collision frequency
	for (int it = 0; it < numCenters; it++)
	{
		//collisionFreq[it] = 1e5; //TODO solve problem using nonzero value for collisionfreq once test case is verified
		collisionFreq[it] = 0.;
	}
	//--------------------------------------------
	

	//-----------------------------------------
	//intialize Device side variables
	double *d_oldTempP, *d_oldTempN, *d_oldVolt, *d_oldDensityP, *d_oldDensityN, *d_oldVelocityRP, *d_oldVelocityRN, *d_oldVelocityZP, *d_oldVelocityZN, *d_tMins, *d_eFieldR, *d_eFieldZ,*d_spectralRadiusR, *d_spectralRadiusZ, *d_collisionFreq, *d_newTempP, *d_newTempN, *d_newVolt, *d_newDensityP, *d_newDensityN, *d_newVelocityRP, *d_newVelocityRN, *d_newVelocityZP, *d_newVelocityZN, *fluxR1,*fluxR2,*fluxR3,*fluxR4,*fluxR5,*fluxR6,*fluxZ1,*fluxZ2,*fluxZ3,*fluxZ4,*fluxZ5,*fluxZ6;
	bool *d_redBlackConvergence;

	cudaMalloc(&d_tMins, sizeCenters);
	cudaMalloc(&d_oldTempP, sizeCenters);
	cudaMalloc(&d_oldTempN, sizeCenters);
	cudaMalloc(&d_oldVolt, sizeEdges);
	cudaMalloc(&d_oldDensityP, sizeCenters);
	cudaMalloc(&d_oldDensityN, sizeCenters);
	cudaMalloc(&d_oldVelocityRP, sizeCenters);
	cudaMalloc(&d_oldVelocityRN, sizeCenters);
	cudaMalloc(&d_oldVelocityZP, sizeCenters);
	cudaMalloc(&d_oldVelocityZN, sizeCenters);
	cudaMalloc(&d_eFieldR, sizeCenters);
	cudaMalloc(&d_eFieldZ, sizeCenters);
	cudaMalloc(&d_spectralRadiusR, sizeCenters);
	cudaMalloc(&d_spectralRadiusZ, sizeCenters);
	cudaMalloc(&d_collisionFreq, sizeCenters);

	cudaMalloc(&d_newTempP, sizeCenters);
	cudaMalloc(&d_newTempN, sizeCenters);
	cudaMalloc(&d_newVolt , sizeEdges);
	cudaMalloc(&d_newDensityP, sizeCenters);
	cudaMalloc(&d_newDensityN, sizeCenters);
	cudaMalloc(&d_newVelocityRP, sizeCenters);
	cudaMalloc(&d_newVelocityRN, sizeCenters);
	cudaMalloc(&d_newVelocityZP, sizeCenters);
	cudaMalloc(&d_newVelocityZN, sizeCenters);
	cudaMalloc(&d_redBlackConvergence,sizeBoolEdges);

	cudaMalloc(&fluxR1, sizeFluxR);
	cudaMalloc(&fluxR2, sizeFluxR);
	cudaMalloc(&fluxR3, sizeFluxR);
	cudaMalloc(&fluxR4, sizeFluxR);
	cudaMalloc(&fluxR5, sizeFluxR);
	cudaMalloc(&fluxR6, sizeFluxR);

	cudaMalloc(&fluxZ1, sizeFluxZ);
	cudaMalloc(&fluxZ2, sizeFluxZ);
	cudaMalloc(&fluxZ3, sizeFluxZ);
	cudaMalloc(&fluxZ4, sizeFluxZ);
	cudaMalloc(&fluxZ5, sizeFluxZ);
	cudaMalloc(&fluxZ6, sizeFluxZ);


	cudaMemcpy(d_oldVolt, volt, sizeEdges, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldTempP, tempP, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldTempN, tempN, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldDensityP, densityP, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldDensityN, densityN, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldVelocityRP, velocityRP, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldVelocityRN, velocityRN, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldVelocityZP, velocityZP, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_oldVelocityZN,velocityZN, sizeCenters, cudaMemcpyHostToDevice);
	cudaMemcpy(d_collisionFreq, collisionFreq, sizeCenters, cudaMemcpyHostToDevice);

	//if starting from initialValues.txt results copy data to device
	if (fresh==0)
	{
		cudaMemcpy(d_eFieldR, eFieldR, sizeCenters, cudaMemcpyHostToDevice);
		cudaMemcpy(d_eFieldZ, eFieldZ, sizeCenters, cudaMemcpyHostToDevice);
		printf("Starting from previous results.\nCopied electric fields to the GPU.\n");
	}
	//-----------------------------------------


	//-----------------------------------------
	//Start Actual Calculation portion
	unsigned int counter = 0;//time loop counter
	bool converged = false;
	if (fresh == 1) {
		do {//first Red-Black method loop to solve for initial voltages
			calc_voltagesRedZeroGradient << < blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_oldDensityP, d_oldDensityN, dR, dZ, d_redBlackConvergence);
			calc_voltagesBlackZeroGradient << <blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_oldDensityP, d_oldDensityN, dR, dZ, d_redBlackConvergence);
			cudaMemcpy(redBlackConvergence, d_redBlackConvergence, sizeBoolEdges, cudaMemcpyDeviceToHost);
			cudaMemcpy(d_oldVolt, d_newVolt, sizeEdges, cudaMemcpyDeviceToDevice);
			for (int it = 0; it < numEdges; it++) {
				if (!redBlackConvergence[it]) break;
				if (it == numEdges - 1) converged = true;
			}
			counter++;
		} while (!converged);
		printf("First calc_voltage call took %d steps.\n", counter);
	}

	else
	{
		cudaMemcpy(d_newVolt, d_oldVolt, sizeEdges, cudaMemcpyDeviceToDevice);
	}


	counter = 0; //Reset the counter

	while (t < Total_Time) { //UNCOMMENT TO RUN TO TOTAL TIME
	//while (counter < 10000) {  //UNCOMMENT TO RUN BY COUNTER. NOTE: REMEMBER TO SET TOTAL TIME much higher than count would every reach to avoid crash
		//Calculate the electric field from voltage
		calc_eFields <<<blocksPerGridCenters, threadsPerBlockDim3 >>> (d_eFieldR, d_eFieldZ, d_newVolt, dR, dZ);

		//Find the smallest time step needed at any location in the domain
		calc_tmins <<<blocksPerGridCenters, threadsPerBlockDim3 >>> (d_oldTempP, d_oldTempN, d_oldVelocityRP, d_oldVelocityRN, d_oldVelocityZP, d_oldVelocityZN, d_tMins, d_spectralRadiusR, d_spectralRadiusZ, dR, dZ);
		cudaMemcpy(tMins, d_tMins, sizeCenters, cudaMemcpyDeviceToHost);
		for (int it = 0; it < numCenters; it++) {//move smallest tMins to tMins[0]
			if (tMins[0] > tMins[it]) tMins[0] = tMins[it];
		}
		double timeStep = tMins[0];

		//if on last step, shorten timestep to end exactly on total_time
		if (timeStep + t > Total_Time) {
			timeStep = Total_Time - t;
		}
		//calc fluxes
		calc_fluxesR << <blocksPerGridFluxRs, threadsPerBlockDim3 >> > (d_oldTempP, d_oldTempN, d_oldDensityP, d_oldDensityN, d_oldVelocityRP, d_oldVelocityRN, d_oldVelocityZP, d_oldVelocityZN, d_eFieldR, d_eFieldZ, fluxR1, fluxR2, fluxR3, fluxR4, fluxR5, fluxR6, d_spectralRadiusR, dR, dZ);
		calc_fluxesZ << <blocksPerGridFluxZs, threadsPerBlockDim3 >> > (d_oldTempP, d_oldTempN, d_oldDensityP, d_oldDensityN, d_oldVelocityRP, d_oldVelocityRN, d_oldVelocityZP, d_oldVelocityZN, d_eFieldR, d_eFieldZ, fluxZ1, fluxZ2, fluxZ3, fluxZ4, fluxZ5, fluxZ6, d_spectralRadiusZ, dR, dZ, t);


		//Calculate temperatures ... for now, this is redundant since they are held constant
		calc_temperatures <<<blocksPerGridCenters, threadsPerBlockDim3 >>> (d_newTempP, d_oldTempP, d_newTempN, d_oldTempN);

		//Calculate the conservation of mass and momentum
		calc_conservatives<<<blocksPerGridCenters, threadsPerBlockDim3 >>> (fluxR1, fluxR2, fluxR3, fluxR4, fluxR5, fluxR6, fluxZ1, fluxZ2, fluxZ3, fluxZ4, fluxZ5, fluxZ6, d_newTempP, d_newTempN, d_newDensityP, d_oldDensityP, d_newDensityN, d_oldDensityN, d_newVelocityRP, d_oldVelocityRP, d_newVelocityRN, d_oldVelocityRN, d_newVelocityZP, d_oldVelocityZP, d_newVelocityZN, d_oldVelocityZN, d_eFieldR, d_eFieldZ, d_collisionFreq, dR, dZ, timeStep);

		//Use the Red-Black iterative method to solve for voltages based on the new charge densities
		converged = false;
		do {
			if (grounded == 0) {
				calc_voltagesRedZeroGradient << < blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_newDensityP, d_newDensityN, dR, dZ, d_redBlackConvergence);
				calc_voltagesBlackZeroGradient << <blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_newDensityP, d_newDensityN, dR, dZ, d_redBlackConvergence);
			}
			else {
				calc_voltagesRedGrounded << < blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_newDensityP, d_newDensityN, dR, dZ, d_redBlackConvergence);
				calc_voltagesBlackGrounded << <blocksPerGridCompressed, threadsPerBlockDim3 >> > (d_newVolt, d_oldVolt, d_newDensityP, d_newDensityN, dR, dZ, d_redBlackConvergence);
			}
			cudaMemcpy(redBlackConvergence, d_redBlackConvergence, sizeBoolEdges, cudaMemcpyDeviceToHost);
			cudaMemcpy(d_oldVolt, d_newVolt, sizeEdges, cudaMemcpyDeviceToDevice);
			for (int it = 0; it < numEdges; it++) {
				if (!redBlackConvergence[it]) break;
				if (it == numEdges - 1) converged = true;
			}
		} while (!converged);

		//update olds from new
		cudaMemcpy(d_oldTempP, d_newTempP, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldTempN, d_newTempN, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldDensityP, d_newDensityP, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldDensityN, d_newDensityN, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldVelocityRP, d_newVelocityRP, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldVelocityRN, d_newVelocityRN, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldVelocityZP, d_newVelocityZP, sizeCenters, cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_oldVelocityZN, d_newVelocityZN, sizeCenters, cudaMemcpyDeviceToDevice);

		//periodically (e.g., every 1000 steps) display a screen message
		counter++;
		t += timeStep;
		if (counter % 10 == 0) {
			printf("Excecuted %d steps using timeStep=%le to get to t=%le.\n", counter, timeStep, t);
		}
	}//end of calculation while loop

	printf("Timestep loop ran for %d steps.\n", counter);
	//-----------------------------------------
	
	//calculate eFields one more time so they update to the most recent voltage
	calc_eFields <<<blocksPerGridCenters, threadsPerBlockDim3 >>> (d_eFieldR, d_eFieldZ, d_newVolt, dR, dZ);
	//-----------------------------------------
	//copy final GPU data back to CPU
	cudaMemcpy(volt,d_newVolt, sizeEdges, cudaMemcpyDeviceToHost);
	cudaMemcpy(tempP,d_newTempP, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(tempN,d_newTempN, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(densityP,d_newDensityP, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(densityN,d_newDensityN, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocityRP,d_newVelocityRP, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocityRN,d_newVelocityRN, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocityZP,d_newVelocityZP, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(velocityZN,d_newVelocityZN, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(eFieldR, d_eFieldR, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(eFieldZ, d_eFieldZ, sizeCenters, cudaMemcpyDeviceToHost);

	cudaMemcpy(spectralRadiusR, d_spectralRadiusR, sizeCenters, cudaMemcpyDeviceToHost);
	cudaMemcpy(spectralRadiusZ, d_spectralRadiusZ, sizeCenters, cudaMemcpyDeviceToHost);

	//-----------------------------------------
	
	//-----------------------------------------
	//print results to results.txt
	FILE* finals;
	finals = fopen("results.txt", "w");
	fprintf(finals,"Voltage ");
	for (int it = 0; it < numEdges; it++)
		fprintf(finals, "%le ",volt[it]);
	fprintf(finals, "\nTemperature_Positive ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", tempP[it]);
	fprintf(finals, "\nTemperature_Negative ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", tempN[it]);
	fprintf(finals, "\nDensity_Positive ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", densityP[it]);
	fprintf(finals, "\nDensity_Negative ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", densityN[it]);
	fprintf(finals, "\nVelocity_Positive_Radial ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", velocityRP[it]);
	fprintf(finals, "\nVelocity_Negative_Radial ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", velocityRN[it]);
	fprintf(finals, "\nVelocity_Positive_ZDirection ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", velocityZP[it]);
	fprintf(finals, "\nVelocity_Negative_ZDirection ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", velocityZN[it]);
	fprintf(finals, "\nEField_Radial ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", eFieldR[it]);
	fprintf(finals, "\nEfield_Zdirection ");
	for (int it = 0; it < numCenters; it++)
		fprintf(finals, "%le ", eFieldZ[it]);
	fprintf(finals, "\nTime_Reached: %le", t);
	fclose(finals);
	clock_t end = clock();
	printf("Program took %le to execute.\n", (float)(end - start) /CLOCKS_PER_SEC);
	//-----------------------------------------

	//As your mom taught you, clean up after yourself
	cudaFree(d_oldTempP);
	cudaFree(d_oldTempN);
	cudaFree(d_oldVolt);
	cudaFree(d_oldDensityP);
	cudaFree(d_oldDensityN);
	cudaFree(d_oldVelocityRP);
	cudaFree(d_oldVelocityRN);
	cudaFree(d_oldVelocityZP);
	cudaFree(d_oldVelocityZN);
	cudaFree(d_tMins);
	cudaFree(d_eFieldR);
	cudaFree(d_eFieldZ);
	cudaFree(d_spectralRadiusR);
	cudaFree(d_spectralRadiusZ);
	cudaFree(d_collisionFreq);
	cudaFree(d_newTempP);
	cudaFree(d_newTempN);
	cudaFree(d_newVolt);
	cudaFree(d_newDensityP);
	cudaFree(d_newDensityN);
	cudaFree(d_newVelocityRP);
	cudaFree(d_newVelocityRN);
	cudaFree(d_newVelocityZP);
	cudaFree(d_newVelocityZN);
	cudaFree(d_eFieldR);
	cudaFree(d_eFieldZ);
	cudaFree(d_redBlackConvergence);
	cudaFree(fluxR1);
	cudaFree(fluxR2);
	cudaFree(fluxR3);
	cudaFree(fluxR4);
	cudaFree(fluxR5); 
	cudaFree(fluxR6);
	cudaFree(fluxZ1);
	cudaFree(fluxZ2);
	cudaFree(fluxZ3);
	cudaFree(fluxZ4);
	cudaFree(fluxZ5);
	cudaFree(fluxZ6);

	return 0;
}//DONE
