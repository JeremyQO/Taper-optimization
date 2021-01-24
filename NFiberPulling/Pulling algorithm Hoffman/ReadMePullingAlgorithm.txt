
BRIEF DESCRIPTION:


Pulling_algorithm_DRUM is a MATLAB script with attached functions which calculates the necessary motor trajectories to achieve a user specified profile for an optical fiber taper design.   The program is designed to work with a heat and pull method of tapering requiring a specified hot zone.  The script also includes a simulation of the expected step-by-step profile of the fiber during the tapering process.  


The script is composed of three sections:  calculation of the pull profile, the simulation of the pull profile , and saving of the pull data.



Section 1:  Calculation of the pull profile



Here the script recursively determines all of the necessary parameters to fully define the fiber profile and motor trajectories during a pull.
  
SubSections A and B
 		
Take the user defined inputs (see HOW TO USE section) to describe the desired profile.

SubSections C, D, E , F and G
		
Here the script initializes, SubSec. C, and then recursively determines the values for all the parameters of interest, SubSec. D. exponential profile, Subsec. E.  linear angle 	1 profile, SubSec. F. linear angle 2 profile, Subsection G ends the loop.

Parameters:

cn			radius ratios between two consecutive steps	
vf			pulling velocity (mm/s)
rwaist		radius of the waist
Lwn			length of the waist
Ln			traveling length of the flame
t0			heating time of a point totally swept by the flame
zt0			traveling length of a point totally swept by the flame
zan			position of the last point totally swept by the flame
zend			position of the end of the fiber
zta 			traveling length of a point that ends in the flame
Lue			difference between extension length of the taper and L0

Each of these elements are vectors of length n, where n is the total number of steps in the 		pull, and the nth element  corresponds to the value during the nth step.    By step we 			mean one sweep of the flame along the fiber.  The flame is swept so that we can control 			the hot zone in each step and created linear taper profiles, and waists of specified length 		Lw.



Section 2:  The simulation



We produce a simulation of the pull that using the parameters calculated in Section 1 to describe the profile of the fiber (using the N = 10000 from Section 1, Subsection B) at the end of each step.  The vectors znew and rnew, have length N, and are the position and radius vectors for the N markers that we break the fiber into. Inside the for loop we update znew and rnew at the end of each step and display it.  It is worth noting, that to simplify the simulation we redefine the zero of the z-axis at the end of every step.  It is defined as the outer edge of the flame, this is achieve by the parameter zflip in the code.


We use XML 210 motors controlled by a Newport XPS controller. The controller provides a means of sending a PVT file to the motors.  The PVT file allows us to specify the the position and velocity of a given motor to be achieved within a specified time.  This information is specified below in an array, the PVTM array:

PVTM 	             a vector [T x1 v1 x2 v2]
T 				duration of a step
x1				relative position motor1 has to reach at T
x2 				relative position motor2 has to reach at T
v1 				velocity motor1 has to reach at T
v2				velocity motor2 has to reach at T

Where each row of PVTM corresponds the nth step of the pull.  In an ideal world the profile of the velocities from step-to-step would be a square wave. Unfortunately, this is unachievable in reality.  To compensate for this we break each 'step' into 5 parts.  Where the first two parts of step n correspond to a ramp 'up' to the desired velocity, part 3 corresponds to the actual pull step and finally parts 4 and 5 correspond to a ramp down to zero velocity.  The ramp 'up' (or 'down') is broken down into two steps to allow for the motors to achieve their different individual velocities in a jerkless manner.  We use the max acceleration allowed by the motors, specified earlier in Section 1, Subsection B, to allow for the smallest deviation from the ideal velocity profile. The PVTM matrix is generated inside the for loop. 

In our pulling apparatus we do not sweep the flame, instead each motor moves with the velocity of the burner (or flame), vb, in each step.  That is, we transform into the rest frame of the flame.  This reduces fluctuations on the flame from the flame moving.  It is for that reason that we add the burner velocity to the velocity of the motors inside the for loop under gfeneration of the PVTM matrix.



Section 3: Saving the Files



The PVTM matrix is saved and can be sent to a controller to define the trajectories the motors need to produce the profile specified in Section 1 Subsection A.  We also save a vector corresponding to the radius of the waist at the end of each step and a vector with the corresponding times at the end of each step.  The final taper profile is saved into the vectors znew and rnew. 

Further Details: 

Simple modifications can be made to the existing sections to add more or fewer linear sections.  All of the vectors generated should be usable as inputs for any set of computer controlled translations stages, although this code is designed to work with an XPS controller, the modifications should be minor.  If one desires to know the profile at the end of each step it would be straightforward to save the profile (znew_n, rnew_n ) for each step of the pull.




HOW TO USE:




Sections 1 and 2 take in all of the necessary user inputs.  Once that is entered simply execute the MATLAB script and the simulation will run a movie displaying the expected profile of the fiber at the end of each step.  Finally the program will save a PVT file (described above) and vectors corresponding to the radius of the waist at the end of each step, the time at the end of each step, and the final taper profile. 


Section 1:

Subsection A: Definition of the target profile

In the section the user sets the external parameters (fiber radius, flame size, and burner velocity)  and chooses the desired inputs (waist length,  waist radius, taper angle, radius to start exponential section, and radius to change from angle 1 to angle 2) .

Here it is necessary to specify:

r0   		fiber-dependent input: the initial radius of the fiber to be tapered in mm 

L0   		nozzle and flow rate dependent input:the size of the flame in mm (this must be 		experimentally found for a fiber)

vb 		nozzle and flow rate dependent inputthe velocity of the burner or flame in mm/s this 		will need to be experimentally 

Lw    	  	user specified input: desired length of the fiber waist

rw          	user specified input: desired final radius of fiber waist

omega   	user specified input:  a vector [angle 1, angle 2]  specified in radians typically 		0.004-0.01 for pulls.  The code could easily be modified to handle more angles.  This 		is the angle from the horizontal that specifies the taper geometry.

r_exp  	the radius in mm at which the  two angle pulls change to an exponential.  The 		exponential section starts with an initial angle matched to the input in omega.  		Typically for sub-micrometer pulls r_exp is set at 0.006.

r1 	  	The radius that angle1 changes to angle 2 in mm

At the moment there are only two angle inputs.  The code is easily modified to handle n angle inputs.

SubSection B:  Algorithm Parameters

The section contains two parts: the definition of the number of points or markers for the simulation, N, and the specification for the max acceleration of the motors.  The values specified for amax  correspond to Newport XML 210 stages.  If using other motor stages it is critical that amax be specified properly for those motors.

We have found that N  =  10000 markers runs quickly and provides enough resolution to have a smooth taper profile for sub-micrometer pulls.  The markers are representative points of the fiber that we track the stretch (position in z) and radius reduction of as  function of step.  

