Explination of Engine Code
By: Travis Widmer and Tom Dale

---------------------------------------------------------------------------------------------------------------------------
Variable names and explinations:


PREPROCESSORS:
	PI 3.14...
	omega 1.5 -- used for Red-Black solver
	kboltz -- Boltzmann's Constant
	epsilon0 --  vacuum permittivity (electric constant)
	q -- elementary charge
	AMU -- one AMU in Kilograms
	massP -- mass of positive ion
	massN -- mass of negative ion
	grid(r,z,nz) (r)*nz+(z) -- conversion from 2D array to 1D array for calculations in GPU
	rCell(i,dR) ------ used to find effective radius of a center value cell
	threadsPerBlock2D 32 -- single dimension for a block creating a 32x32 2D block (1024 threads total)
	omega ---------Red-Black iterative solver coefficient. Must be between 1 and 2. Look for our excel sheet that recommends optimal values
	maxVoltChange -------error tolerance for voltage convergence should be equal to or less than 1e-6.
	CFL ----------time step stabliliy. Should be less than 1. Usually around .75 or .85
	fresh ----- 0 means start program by reading in initialValues.txt file for starting data----1 means start from scratch
	constLeftVolt, constRightVolt, constOuterVolt ------ constant held sides for when code uses voltagesGrounded. Should not be 0
	thrusterTemp -----temerature of the thurster inlet
	minDensity ----- starting density of ion in space
	Total_Time total time to run the engine up to. (check Kernel to make sure CPU loop is running of total_time and not counter amount
	numSegZ,numSegR ------- number of segmentations to break r and z into. This value is the amount of centers, not edges or faces
	lengthZ,lengthR ------- length of map in r and z direction
	plateFrequency ------ the frequency of the alternating charged plates
	thrusterDensity ------ density of exhaust right as it leaves the plate
	plate_Voltage ------ voltage of thruster plate
	thruster_Radius ------- radius for thruster, should be equal or less than lengthR
	numFluxR, numFluxZ ------ used to create rectangular grid for R and Z fluxes. Do not edit
	


counter -------- number of times cpu has gone through main loop and increased totTime
dR,dZ ------length of segments in r and z directions
timeStep ----- calcualted time segmentation
timer ---- time the duration of the simulation

fluxRP* ----flux of positive ions in r direction(face)
fluxRN* ---flux of negative ions in r direction(face)
fluxZP* ----flux of positive ions in z direction(face)
fluxZN* ----flux of negative ions in z direction(face)

oldTempP* ----- old positive ion temp (center)
oldTempN* ------old negative ion temp (center)
newTempP* ----- new positive ion temp (center)
newTempN* ------ new negative ion temp (center)

oldVolt* ------old voltage (corner)
newVolt* ------new voltage (corner)

eFieldR* -- Electric field in the r direction (center)
eFieldZ* -- Electric field in the z direction (center)

oldDensityP*----old positive ion density (center)
oldDenstiyN*----old negative ion denstiy (center)
newDensityP*----new positive ion density (center)
newDenstiyN*----new negative ion denstiy (center)

oldVelocityRP*--old positive velocity in r direction (center)
oldVelocityRN*--old negative velocity in r direction (center)
oldVelocityZP*--old positive velocity in z direction (center)
oldVelocityZN*--old negative velocity in z direction (center)
newVelocityRP*--new positive velocity in r direction (center)
newVelocityRN*--new negative velocity in r direction (center)
newVelocityZP*--new positive velocity in z direction (center)
newVelocityZN*--new negative velocity in z direction (center)

sizeCenter -- number of elements (in bytes) of center values
sizeEdge -- number of elements (in bytes) of edge values (i.e. corner and flux)

dim3 blocksPerGridCenters((numSegZ/threadsPerBlock)+1, (numSegR/threadsPerBlock)+1) 
		-- number of blocks based on number of segments in r and z direction while adding a row and a column to prevent exceeding Grid size

dim3 blocksPerGridEdges(((numSegZ+1)/threadsPerBlock)+1, ((numSegR+1)/threadsPerBlock)+1) 
		-- number of blocks based on number of segments in r and z direction while adding a row and a column to prevent exceeding Grid size

dim3 threadsPerBlock(32, 32) -- creates a 32x32 block using all 1024 threads per block

tMins* -- minimum time for each center
spectralRadiusR* -- "net propultion" in r direction ('velocity' - 'speed of sound') (center)
spectralRadiusR* -- "net propultion" in z direction ('velocity' - 'speed of sound') (center)
collisionFreq* -- calculated using densities and temperature (center)


NOTATION:
	i -- position in r on Map (r = deltaR*i)
	j -- position in z on Map (z = deltaZ*j) 

MAP:
  o--> z
  |    "INNER"
  v       |corner(i,j)            | corner(i,j+1)
  r     --O----------FR-----------O--
	  |       fluxR(i,j)      |
	  |                       |
	  |                       |
	  |       center(i,j)     |
	 FZ           X           FZ  fluxZ(i, j+1)
	  |fluxZ(i,j)             |
	  |                       |
	  |                       |
	  |       fluxR(i+1,j)    |
 "LEFT" --O----------FR-----------O-- "RIGHT"
	  |corner(i+1,j)          | corner(i+1,j+1)
       "OUTER"

FUNCTIONS:
		
	calc_voltage<<< blocksPerGridEdges, threadsPerBlock >>>(oldDensityP*, oldDenistyN*, numSegR, numSegZ, dR, dZ)
		-- calculate Voltage at the corners based on the ion densities using the Laplace equation
		-- Uses: oldDensityP*, oldDenistyN*, elementary charge, number of segments, dR, dZ
		-- Updates: new voltage at each corner
		
	calc_conservation<<< blocksPerGridCenters, threadsPerBlock >>>()
		-- calculate every conservation equation (6 total)
			Conservation of Mass (2) - calculate positive and negative ion densities, seperately, using radius position, oldDensityP, oldDensityN, and all old velocities (in all directions) in a hyperbolic equation
				- Uses: old values of density, old values of velocity in R and Z, 
				- Updates: positive and negative ion denisties into newDensityP* and newDensityN*
			Conservatoin of Momentum (4) - calculate positive and negative ion velocities in r and z directions (all seperately) using oldDensites, oldVelocities, massP, massN, elementary charge, boltzmanns constant, collision freq and Temperature in a hyperbolic equation
				-Updates: positive and negative ion velocities in both R and Z directions into newVelocityRN*, newVelocityRP*, newVelocityZN*, and newVelocityZP*
			
	calc_tMins <<<blocksPerGridCenters, threadsPerBlock >>>()
			--calculate smallest deltaT (to ensure stability) for each block which will then be compared across blocks where the smallest will be used as the next timeStep

	calc_temperatures<<<blocksPerGridCenters, threadPerBlock>>>(oldTempP*,oldTempN*,newTempP*,newTempN*)
			--updates new temperatures, for initial version of code keep constant
	
	calc_eFields<<<blocksPerGridCenters,threadsPerBlock>>>(eFieldR*, eFieldZ*, volt*, dR,dZ,numSegR,numSegZ)
			-- E-Fields are the negative gradient of the Voltage
	
	calc_Fluxes<<<>>>(...)
			-- Global variables for fluxes through edges of cells are calculated here and used in conservatives

	
			
			
			
			
			