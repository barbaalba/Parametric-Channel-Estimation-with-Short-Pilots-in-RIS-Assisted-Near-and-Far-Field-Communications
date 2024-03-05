# NearField_RIS

# Tutorial Functions
- BeamWidth.m: It plots a figure where we understand that moving toward edges the beam gets wider and stretched in near field region. 
- EigenDecompositionNearField.m: It plots the Eigenvalues in order.
- phaseOverAntennas.m: It demonstrate how near field project phase over the array surface compared to far field. 
- stretchedBeam.m: It plots the orthogonal crosses to show that increasing elevation leads to points going out of the border.
- Tightness.m: To understand far-field and first order approximation of the channel. 

# Improtant Functions
- widebeam.m: To generate Figures 6 and 7 in the paper.
- wideAlgTest.m: To generate Figures 8 and 9.
- AlgTestOptimized.m: To generate Figures 10 and 11. It needs  to be configured to correspond to either far-field or near-field scenarios.
- randomwalkRun.m: To generate Figure 12.
- AlgTest.m: The main code to test the efficiency and effectiveness of the generic channel estimation algorithm. It generates the figures [] and []
- AzElGraph.m: It plots set of angles that leads to orthogonal beams. Figure 3 in the paper. 
- BestCodeBook.m: It plots 3D figure where we understand what is the best initialization angle to generate the code book. Figure [] in the paper.
- MLE.m: It is the function that estimate the channel by assuming that it is approximately far-field. 
- MLE3D.m: It estimate the channel by estimating distance and azimuth-elevation pairs. It is work the best for both near field and far-field region.
- nearFieldChan.m: The function that generate realistic channel based on the propagation distance between user and RIS elements. 
- UPA_BasisElupnew.m: It find all the angle pairs that results in orthogonal beams.
- UPA_Codebook.m: Generate array responses of the orthogonal angle pairs.
- ChanParGen.m: 
- distanceCorr.m: 
  
# Random Walk functions
- plotTrajectory.m: plot the trajectory of the random walk in the room 
- randomwalk.m: It generate a random walk scenario in a confined room with given size 
- realChanGen.m: 
- Run.m: 
 
