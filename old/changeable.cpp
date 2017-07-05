#define EPS 1.0e-05
#include "cluster.hpp"

#define NONE -1
#define SQRT -2

/*****************************************************************
generate_CollocationPts - For generating Collocation Points given the Start, End and the number of points needed.
dStart		-	Double. It is a temporaray variable denoting the start of the collocation points.
dEnd			-	Double. It is a  temporaray variable denoting the end of the collocation points.
iNumberofPts	-	Integer. It is a temporary variable denoting the number of points in the range provided by dStart and dEnd variables.
pvecdCollocationPts	-	Pointer to Vector Double. It is a pointer to a Vector Double storing the Collocation points in a single dimenstion.
iOperation			-	Integer. It is as an operation to be performed on the collocation points.
					The collocation point are initially created as equally distributed points between given range by iStart and iEnd variables.
					Then these collocations points can transformed to some other points by this iOperation variable.
					Currently only one such function has been implemented. 
					SQRT is one of such function. Its is defined to be -2 in cluster.h
					In the below a function switch is used to test for value of iOperation. New cases can be introduced for using new transformation on the equally distributed points in the given range.
*****************************************************************/	
	void clsClusteringApproximation::generate_CollocationPts(double dStart, double dEnd, int iNumberofPts, vector<double>* pvecdCollocationPts, int iOperation){
		int iLoopVari;
		double dStep;
		
		printf("Creating Colloc Pts.......");
		dStep = (dEnd - dStart)/(iNumberofPts - 1);
		pvecdCollocationPts->resize(iNumberofPts);
		for (iLoopVari = 0; iLoopVari < iNumberofPts; iLoopVari++){
			switch(iOperation){
				case NONE :
					(*pvecdCollocationPts)[iLoopVari] = dStart + (iLoopVari*dStep);
					break;
				case SQRT :
					(*pvecdCollocationPts)[iLoopVari] = sqrt(dStart + (iLoopVari*dStep));
					break;
			}
		}
		printf("done\n");
	}

/*****************************************************************
admissible_ClusterPair - This function is used to define the criteria for admissible clusters.
				The function should return True in case the given Pair is admissible. otherwise False.
iIndex1	-	Integer. It is used to denote the Cluster Linear forming the Cluster Pair.  The iIndex1 provides the loaction of the Cluster Linear in the vector containing all the Cluster Linear of 1st Dimesnion.
iIndex2	-	Integer. It is used to denote the Cluster Linear forming the Cluster Pair.  The iIndex2 provides the loaction of the Cluster Linear in the vector containing all the Cluster Linear of 2nd Dimesnion.
******************************************************************
Return Type	:	Boolean
Returns		:	True if Admissible Pair else False.
*****************************************************************/	
	bool clsClusteringApproximation::admissible_ClusterPair(int iIndex1, int iIndex2){
		/// ///////////////////////////////////////// Put admissibility condition here ///////////////////////////////////////// /// 
		double dA,dB, dC, dD, dDist, dTemp1, dTemp2;
		double dLengthTemp, dBreadthTemp;

		dA = vec2dCollocationPts[0][vec2pclsClusterLinear[0][iIndex1]->get_iStart()];
		dB = vec2dCollocationPts[0][vec2pclsClusterLinear[0][iIndex1]->get_iEnd()];
		dC = vec2dCollocationPts[1][vec2pclsClusterLinear[1][iIndex2]->get_iStart()];
		dD = vec2dCollocationPts[1][vec2pclsClusterLinear[1][iIndex2]->get_iEnd()];
		dLengthTemp = dB - dA;
		dBreadthTemp = dD - dC;
		dTemp1 = dC - dB;
		dTemp2 = dA - dD;
		dDist = max(0.0, max(dTemp1, dTemp2));
		if (dAdmissibilityCoefficient*dDist >= max(dLengthTemp,dBreadthTemp) ){
			return (true);
            }
        else{
			return (false);
            }

    }

/*****************************************************************
kernel_function - This function is used to define the kernel function.
dX	-	Double. It is used to denote the co-ordiante in the first Dimension at which the kernel function's value is to be calculated.
dY	-	Double. It is used to denote the co-ordiante in the second Dimension at which the kernel function's value is to be calculated.
******************************************************************
Return Type	:	Double
Returns		:	The value of the function at the given co-ordinates.
*****************************************************************/	
	double clsClusteringApproximation::kernel_function(double dX, double dY){
		/// ///////////////////////////////////////// Write the kernel function here ///////////////////////////////////////// /// 
		if (fabs(dX - dY) < EPS)
			return(0);
		else
			return (1.0/(fabs(dX - dY)));

       
	}

/*****************************************************************
get_iOperation - It is the function to determine the type of Operation to be used for generating collocation points.
Agruments:
iDimension 	-	Integer. This variable determines Dimension of Cluster Linear. 
******************************************************************
Return Type	:	Integer
Returns		:	The type of Operation. This has to be defined in the start of this file as an Integer.
*****************************************************************/
	int clsClusteringApproximation::get_iOperation(int iDimension){
		if (iDimension == 0)
			return (SQRT);
		else
			return (SQRT);
	}

/*****************************************************************
get_iNumberofTimes - It is the function to determine the number of times the process is to be repeated for calculating average time complexity.
******************************************************************
Return Type	:	Integer
Returns		:	The Number of times the process has to be repeated.
*****************************************************************/
	int clsClusteringApproximation::get_iNumberofTimes(){
		return(1);
	}
