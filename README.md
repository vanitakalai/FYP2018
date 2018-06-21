# FYP2018
## Final Year Project: Fault Tolerant Control
## Supervisor: Dr Imad Jaimoukha


This repository contains all the code used in this project along with the final report. Please refer to the report for more information. 

### File Clarifications:

* `sysgen` -> creates matrix parameters for observer e.g. A, B, C
* `faultsgen` -> function that creates fault matrix, Delta
* `testScript` -> used with FaultTolerantLS and FaultTolerantRPP for simulation testing 
* `testLoop` -> used to implement loop testing (e.g. to run code 100 times and collect data)
* plotRPP -> function that creates RPP plot region using r, alpha and theta

* FaultTolerantLoop -> iterative algorithm with LS
* FaultTolerantLS -> NFTC, FTC and FFTC design for H infinity with LS method
* FaultTolerantRPP -> NFTC, FTC and FFTC design for H infinity with RPP method
* binarySearchTheta -> NFTC, FTC and FFTC design for binary search for theta with RPP method

All files designed for loop testing
..* pConstrainedLS -> NFTC, FTC and FFTC design for H infinity with LS method with p constraints
..* pConstrainedRPP -> NFTC, FTC and FFTC design for H infinity with RPP method with p constraints
..* binarySearchThetawithP -> NFTC, FTC and FFTC design for binary search for theta with RPP method with p constraints

..* linkedSensorsFTC -> FTC design for linked sensors implementation for all design methods
..* linkedSensorsFTC -> FFTC design for linked sensors implementation for all design methods

..* HFTC -> FTC and FFTC design for hybrid design method

### Models:

..* Model2016a -> switching state observer model from no faults to faults in Simulink version 2016a
..* modelext -> switching state observer model from no faults to faults in Simulink version 2017a (for hybrid FTC but can also be used for all FTC)
