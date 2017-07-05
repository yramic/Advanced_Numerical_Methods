#define EPS2 1.0e-16
#define PI 3.14159
#include "cluster.hpp"
#include <cstdlib>
#include <cstring>

/*****************************************************************
Functions of Class Cluster Linear
*****************************************************************/

/*****************************************************************
clsClusterLinear - It is the constructor for class ClusterLinear.
Arguments:
iTempStart	-	Integer. It deontes the start of the Linear Cluster. This value is stored in iStart variable of class Cluster Linear.
iTempEnd	-	Integer. It deontes the end of the Linear Cluster. This value is stored in iEnd variable of class Cluster Linear.
iTempLeft	-	Integer. It deontes the Left Child of the Linear Cluster. This value is stored in iLeftChild variable of class Cluster Linear.
			This value indicates the left child position of the this Cluster Linear in a vector of Cluster Linear.
iTempRight	-	Integer. It deontes the Left Child of the Linear Cluster. This value is stored in iLeftChild variable of class Cluster Linear.
			This value indicates the left child position of the this Cluster Linear in a vector of Cluster Linear.
iDegree	-	This value denotes the degree. It is used for reserving space for vecdW variable of Cluster Linear.
*****************************************************************/
	clsClusterLinear::clsClusterLinear(int iTempStart, int iTempEnd, int iTempLeft, int iTempRight, int iDegree){
        iStart = iTempStart;
		iEnd = iTempEnd;
		iLeft = iTempLeft;
		iRight = iTempRight;
		vecdW.reserve(iDegree);
	}

/*****************************************************************
~clsClusterLinear - It is the destructor for class ClusterLinear.
*****************************************************************/	
	clsClusterLinear::~clsClusterLinear() {
	}

/*****************************************************************
increase_size_LagValAtCollocPtsInCluster - It is the function to increase size of LagValAtCollocPtsInCluster variable of class ClusterLinear.
Arguments:
iDegree	-	This value denotes the degree. It is used for sizing the vector LagValAtCollocPtsInCluster.
*****************************************************************/	
	void clsClusterLinear::increase_size_LagValAtCollocPtsInCluster(int iDegree){
        int iLoopVari; 
		
		vec2dLagValAtCollocPtsInCluster.resize((iEnd - iStart) + 1);
		for (iLoopVari = 0; iLoopVari < vec2dLagValAtCollocPtsInCluster.size(); iLoopVari++){
			vec2dLagValAtCollocPtsInCluster[iLoopVari].resize(iDegree + 1);		
		}
	}

/*****************************************************************
get_iStart - It is the function to return the iStart variable of the class ClusterLinear.
Return Type	- 	int.
Returns		-	iStart variable value
*****************************************************************/		
	int clsClusterLinear::get_iStart(){
		return iStart;
	}

/*****************************************************************
get_iEnd - It is the function to return the iEnd variable of the class ClusterLinear.
Return Type - 	int.
Returns 	- 	iEnd variable value
*****************************************************************/	
	int clsClusterLinear::get_iEnd(){
		return iEnd;
	}

/*****************************************************************
get_iLeft - It is the function to return the iLeft variable of the class ClusterLinear.
Return Type	-	int.
Returns		-	iLeft variable value
*****************************************************************/
	int clsClusterLinear::get_iLeft(){
		return iLeft;
	}	

/*****************************************************************
get_iRight - It is the function to return the iRight variable of the class ClusterLinear.
Return Type	-	int.
Returns		-	iRight variable value
*****************************************************************/	
	int	clsClusterLinear::get_iRight(){
		return iRight;
	}

/*****************************************************************
get_ptr_LagValAtCollocPtsInCluster - It is the function to return the pointer to LagValAtCollocPtsInCluster variable of the class ClusterLinear.
Return Type	-	Pointer to Vector of Vector Double
Returns		-	Pointer to LagValAtCollocPtsInCluster variable
*****************************************************************/	
	vector< vector <double> >* clsClusterLinear::get_ptr_LagValAtCollocPtsInCluster(){
		return &vec2dLagValAtCollocPtsInCluster;
	}

/*****************************************************************
get_ptr_AppearsIn - It is the function to return the pointer to AppearsIn variable of the class ClusterLinear.
Return Type	-	Pointer to Vector Integer
Returns		-	 Pointer to AppearsIn variable
*****************************************************************/	
    vector<int>* clsClusterLinear::get_ptr_AppearsIn(){
		return &veciAppearsIn;
    }

/*****************************************************************
get_ptr_vecdW - It is the function to return the pointer to vecdW variable of the class ClusterLinear.
Return Type - 	Pointer to Vector Double
Returns 	- 	Pointer to vecdW variable
*****************************************************************/    
	vector<double>* clsClusterLinear::get_ptr_vecdW(){
		return &vecdW;
	}

/*****************************************************************
print  - It is the function to prnt the class ClusterLinear.
Arguments:
iIndexNumber	-	Integer. This variable is used to indicate the position of this instance of Cluster Linear in the vector of Cluster Linear.
*****************************************************************/  	
	void clsClusterLinear::print(int iIndexNumber){
		int iLoopVari,iLoopVarj;
		
		printf("\nLinear Cluster Number = %d\n", iIndexNumber);
		printf("\tStart = %d\n", iStart);
		printf("\tEnd = %d\n", iEnd);
		printf("\tLeft = %d\n", iLeft);
		printf("\tRight = %d\n", iRight);
		printf("\tLagrange Values at Collocation Pts in Cluster:\n\t\tIndex of Collocation Pt \tLagrange Value");
		for (iLoopVari = 0; iLoopVari < vec2dLagValAtCollocPtsInCluster.size(); iLoopVari++){
			printf("\n\t\t       %d",iStart + iLoopVari);
			for (iLoopVarj = 0; iLoopVarj < vec2dLagValAtCollocPtsInCluster[iLoopVari].size(); iLoopVarj++){
				printf("\t%4.3lf", vec2dLagValAtCollocPtsInCluster[iLoopVari][iLoopVarj]);
			}
		}
		printf("\n\tAppears in: ");
		for (iLoopVari = 0; iLoopVari < veciAppearsIn.size(); iLoopVari++){
			if (veciAppearsIn[iLoopVari] < 0)
				printf("Near");
			else
				printf("Far");
			printf("%d\t", abs(veciAppearsIn[iLoopVari]) - 1);
		}
		printf("\n\tW Values: ");
		for (iLoopVari = 0; iLoopVari < vecdW.size(); iLoopVari++){
			printf("%4.3lf\t", vecdW[iLoopVari]);
		}
		printf("\n");
	}


	
/*****************************************************************
Functions of Class Cluster Pair
*****************************************************************/

/*****************************************************************
clsClusterPair - It is the constructor for class ClusterLinear.
Arguments:
iDimension			-	Integer. It deontes number of Dimensions.
pclsClusterLinear1	-	Pointer to class Cluster Linear. It deontes the ponter to Class Cluster Linear for Dimension 1 which forms this Cluster Pair.
pclsClusterLinear2	-	Pointer to class Cluster Linear. It deontes the ponter to Class Cluster Linear for Dimension 2 which forms this Cluster Pair.
iRows			-	Integer. It deontes the  no. of rows in vec2dFnValAtChebyNodeInCluster of class Cluster Pair.
iCols				-	Integer. It deontes the  no. of columns in vec2dFnValAtChebyNodeInCluster of class Cluster Pair.
*****************************************************************/
	clsClusterPair::clsClusterPair(int iDimension, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2, int iRows, int iCols){
		int iLoopVari;
        vecpclsClusterLinear.resize(iDimension);
		vecpclsClusterLinear[0] = pclsClusterLinear1;
		vecpclsClusterLinear[1] = pclsClusterLinear2;
		vec2dFnValAtChebyNodeInCluster.resize(iRows);
		for (iLoopVari = 0; iLoopVari < iRows; iLoopVari++){
			vec2dFnValAtChebyNodeInCluster[iLoopVari].reserve(iCols);
		}
	}	

/*****************************************************************
~clsClusterPair - It is the destructor for class ClusterLinear.
*****************************************************************/	
	clsClusterPair::~clsClusterPair() {
	}

/*****************************************************************
get_ptr_LagValAtCollocPtsInCluster - It is the function to return the pointer to FnValAtChebyNodeInCluster variable of the class ClusterLinear.
Return Type - 	Pointer to Vector of Vector Double
Returns 	- 	Pointer to FnValAtChebyNodeInCluster variable
*****************************************************************/		
	vector< vector <double> >* clsClusterPair::get_ptr_FnValAtChebyNodeInCluster(){
		return &vec2dFnValAtChebyNodeInCluster;
	}

/*****************************************************************
get_ptr_ClusterLinear - It is the function to return the pointer to Class ClusterLinear which forms this Cluster Pair.
Agruments:
iDimension 	-	Integer. This variable determines which Cluster Linear is to be returned. 
******************************************************************
Return Type - 	Pointer to Class Cluster Linear
Returns 	-	Pointer to ClusterLinear
*****************************************************************/	
	clsClusterLinear* clsClusterPair::get_ptr_ClusterLinear(int iDimension){
		return vecpclsClusterLinear[iDimension];
	}

/*****************************************************************
print  - It is the function to prnt the class Cluster Pair.
Arguments:
iIndexNumber	-	Integer. This variable is used to indicate the position of this instance of Cluster Pair in the vector of Cluster Pair.
iDegree		-	Integer. This variable is used to indicate the degree of the Clustering Approximation technique. 
iDimension		-	Integer. This variable is used to indicate the Dimension of the Clustering Approximation technique.
*****************************************************************/	
	void clsClusterPair::print(int iIndexNumber, int iDegree, int iDimension){
		int iLoopVari, iLoopVarj;
		
		printf("\nCluster Pair Number = %d\n", iIndexNumber);
		for (iLoopVari = 0; iLoopVari < iDimension; iLoopVari++){
			printf("Linear Cluster in %d Dimension", iLoopVari); 
			vecpclsClusterLinear[iLoopVari]->print(-1);
		}
		printf("\tFunction Values at Chebyshev node in Cluster Pair:\n");
		for (iLoopVari = 0; iLoopVari < vec2dFnValAtChebyNodeInCluster.size(); iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < vec2dFnValAtChebyNodeInCluster[iLoopVari].size(); iLoopVarj++){
				printf("\t%4.3lf", vec2dFnValAtChebyNodeInCluster[iLoopVari][iLoopVarj]);	
			}
			printf("\n");
		}
	}


	
/*****************************************************************
Functions of Class Clustering Approximation
*****************************************************************/

/*****************************************************************
clsClusteringApproximation - It is the constructor for class ClusterApproximation.
Arguments:
iDimension	-	Integer. It deontes the Dimension of the Clustering Approximation Technique.
*****************************************************************/
	clsClusteringApproximation::clsClusteringApproximation(int iDimension){
        this->iDimension = iDimension;
        veciNumberofPts.resize(iDimension);
		vec2dCollocationPts.resize(iDimension);
		vec2pclsClusterLinear.resize(iDimension);
	}

/*****************************************************************
set_DirPath - It is used for setting the Directory Path for storing files to be used by MATLAB to plot data.
Arguments:
pchDirPath	-	Pointer to character. It is a pointer to character which contains Directory Path.
*****************************************************************/	
	void clsClusteringApproximation::set_DirPath(char* pchDirPath){
		char* pchTemp = new char [200];
		strcpy(pchTemp, "mkdir ");
        
        this->pchDirPath = new char[200];
        strcpy(this->pchDirPath, pchDirPath);
        strncat(this->pchDirPath, "/Matlab_Files/", 14);
        strcat(pchTemp, this->pchDirPath);
        system(pchTemp);
    }

/*****************************************************************
write_informationMain - It is used for writing the main Information File which is used by MATLAB. This file contains the all information about Clustering Approximation.
Arguments:
iTotalPoints		-	Integer. This variable denotes Total Number of Different Collocation Paints worked on in the Clustering Approximation Techniques.
iTotalDegree	-	Integer. This variable denotes Total Degrees worked on  in the Clustering Approximation Techniques.
iTotalAdmissible	-	Integer. This variable denotes Total Admissible Values worked on in the Clustering Approximation Techniques.
*****************************************************************/	
	void clsClusteringApproximation::write_informationMain(int iTotalPoints, int iTotalDegree, int iTotalAdmissible){
		char* pchName = new char[80];
		ofstream fileInformationMain;	
		
		printf("Writing File \"InformationMain.txt\" .......\n");
		strcpy(pchName, pchDirPath);
		strncat(pchName, "InformationMain.txt", 19);
		fileInformationMain.open(pchName);
		fileInformationMain << iTotalPoints << "," << iTotalDegree << "," << iTotalAdmissible << "\n";
		fileInformationMain.close();
		printf("done\n");
	}

/*****************************************************************
write_information - It is used for writing the Information File which is used by MATLAB. 
				This file contains information regarding a Clustering Technique for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation.
*****************************************************************/	
	void clsClusteringApproximation::write_information(char* pcharFileName){
		ofstream fileInformation;	
		int iMaxofNumberofPts;
		int iLoopVari;
		
		printf("Writing File \"Information.txt\" .......\n");
		fileInformation.open(pcharFileName);
		fileInformation << veciNumberofPts[0] << "," << veciNumberofPts[1] << "," << veciNumberofPts[1] << ",1,1\n";
		iMaxofNumberofPts = max(veciNumberofPts[0],veciNumberofPts[1]);
		for (iLoopVari = 0; iLoopVari < iMaxofNumberofPts; iLoopVari++){
			if (iLoopVari < vec2dCollocationPts[0].size())
			   fileInformation << vec2dCollocationPts[0][iLoopVari]; 
			fileInformation << ",";
			if (iLoopVari < vec2dCollocationPts[1].size())
			   fileInformation << vec2dCollocationPts[1][iLoopVari] << "," << vecdVal[iLoopVari];
			if (iLoopVari == 0){
				fileInformation << "," << iDegree << "," << dAdmissibilityCoefficient;
				fileInformation << "," << vec2pclsClusterLinear[0].size() << "," << vec2pclsClusterLinear[1].size();
				fileInformation << "," << vecpclsClusterPairNear.size() << "," << vecpclsClusterPairFar.size(); 
			}
			fileInformation << "\n"; 
		}
		fileInformation.close();
		printf("done\n");
	}

/*****************************************************************
write_clusterPairNear - It is used for writing the Cluster Pair Near which is used by MATLAB for plotting Cluster Pair. 
				This file contains information regarding Cluster Pair Near for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation for Cluster Pair Near.
*****************************************************************/
	void clsClusteringApproximation::write_clusterPairNear(char* pcharFileName){
		ofstream fileClusterPair;
		int iLoopVari, iLoopVarj, iLoopVark;
		int iIndexX1, iIndexX2, iIndexY1, iIndexY2;
		double dX1, dX2, dY1, dY2; 
		
		printf("Writing File \"ClusterPairNear.txt\" .......\n");
		fileClusterPair.open(pcharFileName);
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairNear.size(); iLoopVari++){
			iIndexX1 = (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(0))->get_iStart();
			iIndexX2 = (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(0))->get_iEnd();
			iIndexY1 = (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(1))->get_iStart();
			iIndexY2 = (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(1))->get_iEnd();
			for (iLoopVarj = iIndexX1; iLoopVarj <= iIndexX2; iLoopVarj++){
				for (iLoopVark = iIndexY1; iLoopVark <= iIndexY2; iLoopVark++){
					dX1 = vec2dCollocationPts[0][iLoopVarj];
					dY1 = vec2dCollocationPts[1][iLoopVark];
					fileClusterPair << dX1 << "," << dY1 << "\n";				
				}
			}
		}
		fileClusterPair.close();
		printf("done\n");		
	}

/*****************************************************************
write_clusterPairFar - It is used for writing the Cluster Pair Near which is used by MATLAB for plotting Cluster Pair. 
				This file contains information regarding a Cluster Pair Far for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation for Cluster Pair Far.
*****************************************************************/
	void clsClusteringApproximation::write_clusterPairFar(char* pcharFileName){
		ofstream fileClusterPair;
		int iLoopVari;
		double dX1, dX2, dY1, dY2; 
		
		printf("Writing File \"ClusterPairFar.txt\" .......\n");
		fileClusterPair.open(pcharFileName);
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairFar.size(); iLoopVari++){
			dX1 = vec2dCollocationPts[0][(vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(0))->get_iStart()];
			dX2 = vec2dCollocationPts[0][(vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(0))->get_iEnd()];
			dY1 = vec2dCollocationPts[1][(vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(1))->get_iStart()];
			dY2 = vec2dCollocationPts[1][(vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(1))->get_iEnd()];
			fileClusterPair << dX1 << "," << dX2 << "," << dY1 << "," << dY2 << "\n";
		}
		fileClusterPair.close();
		printf("done\n");		
	}

/*****************************************************************
write_clusterLinear - It is used for writing the Cluster Linear which is used by MATLAB for animating Cluster Pair. 
				This file contains information regarding all the Cluster Linear for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation for Cluster Linear.
iDimension		-	Integer. This is used to indicate the Dimension of the Cluster Linear to be written.
*****************************************************************/	
	void clsClusteringApproximation::write_clusterLinear(char* pcharFileName, int iDimension){
		ofstream fileClusterLinear;
		int iLoopVari;
		
		printf("Writing File \"ClusterLinear%d.txt\" .......\n", iDimension);
		fileClusterLinear.open(pcharFileName);
		for (iLoopVari = 0; iLoopVari < vec2pclsClusterLinear[iDimension].size(); iLoopVari++){
			fileClusterLinear << vec2pclsClusterLinear[iDimension][iLoopVari]->get_iStart() << ",";
			fileClusterLinear << vec2pclsClusterLinear[iDimension][iLoopVari]->get_iEnd() << ",";
			fileClusterLinear << vec2pclsClusterLinear[iDimension][iLoopVari]->get_iLeft() << ",";
			fileClusterLinear << vec2pclsClusterLinear[iDimension][iLoopVari]->get_iRight() << "\n";
		}
	}

/*****************************************************************
write_ClusterLinearInClusterPairFar - It is used for writing the Cluster Linear which appear in Cluster Pair Far is used by MATLAB for animating Cluster Pair. 
				This file contains information regarding all the Cluster Linear which appear in Cluster Pair Far for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation for Cluster Linear which appear in Cluster Pair Far.
*****************************************************************/	
	void clsClusteringApproximation::write_clusterLinearInClusterPairFar(char* pcharFileName){
		ofstream fileFinalClusterLinear;
		int iLoopVari;
		
		printf("Writing File \"Final ClusterLinear in ClusterPairFar.txt\" .......\n", iDimension);
		fileFinalClusterLinear.open(pcharFileName);
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairFar.size(); iLoopVari++){
			fileFinalClusterLinear << (vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(0))->get_iStart() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(0))->get_iEnd() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(1))->get_iStart() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(1))->get_iEnd() << "\n";
		}	
	}

/*****************************************************************
write_ClusterLinearInClusterPairNear - It is used for writing the Cluster Linear which appear in Cluster Pair Near is used by MATLAB for animating Cluster Pair. 
				This file contains information regarding all the Cluster Linear which appear in Cluster Pair Near for a particular Number of Collocation points, Degree, and Admissibility Coefficient 
Arguments:
pcharFileName	-	Pointer to Character. This is a Pointer to Character which stores the File Name for writing the information of this Clustering Approximation for Cluster Linear which appear in Cluster Pair Near.
*****************************************************************/	
	void clsClusteringApproximation::write_clusterLinearInClusterPairNear(char* pcharFileName){
		ofstream fileFinalClusterLinear;
		int iLoopVari;
		
		printf("Writing File \"Final ClusterLinear in ClusterPairNear.txt\" .......\n", iDimension);
		fileFinalClusterLinear.open(pcharFileName);
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairNear.size(); iLoopVari++){
			fileFinalClusterLinear << (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(0))->get_iStart() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(0))->get_iEnd() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(1))->get_iStart() << ",";
			fileFinalClusterLinear << (vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(1))->get_iEnd() << "\n";
		}	
	}
	
/*****************************************************************
put_NumberofPts - It is used for inputting the Number of points in each dimension for the Clustering Technique.
				It is only used by the function initialize process.
*****************************************************************/	
	void clsClusteringApproximation::put_NumberofPts() {
        int iLoopVari;
		
		for (iLoopVari = 0; iLoopVari < iDimension; iLoopVari++){
			printf("Enter Number of Collocation Points in %d Dimension : ", iLoopVari + 1);
			scanf("%d", &veciNumberofPts[iLoopVari]);
			vec2dCollocationPts[iLoopVari].resize(veciNumberofPts[iLoopVari]);
		}
		initialize_vector(&vecdFinalAnswer, veciNumberofPts[0], 0);
	}

/*****************************************************************
put_Limits - It is used for inputting the limits of Collocation points for the Clustering Technique.
				It is only used by the function generate and generate_different.
Arguments:
iIndex	-	Integer. It is used for defining the Dimension for the Collocation Points to be inputed.
pdStart	-	Pointer to double. It is a Pointer to the Double variable for denoting the start of the Collocation Points.
pdEnd	-	Pointer to double. It is a Pointer to the Double variable for denoting the end of the Collocation Points.
*****************************************************************/	
	void clsClusteringApproximation::put_Limits(int iIndex, double* pdStart, double* pdEnd){
		printf("Enter the start for Collocation Points for %d dimension ? ", iIndex + 1);
		scanf("%lf", pdStart);
		printf("Enter the end for Collocation Points for %d dimension ? ", iIndex + 1);
		scanf("%lf", pdEnd);
	}

/*****************************************************************
put_Degree - It is used for inputting the Degree for the Clustering Technique.
				It is only used by the function initialize process.
*****************************************************************/	
	void clsClusteringApproximation::put_Degree() {
		printf("Enter degree of Lagrange Polynomial : ");
		scanf("%d", &iDegree);
	}

/*****************************************************************
put_CollocationPts - It is used for inputting the Collocation points for the Clustering Technique.
				It is only used by the function initialize process.
*****************************************************************/	
	void clsClusteringApproximation::put_CollocationPts() {
		int iLoopVari, iLoopVarj;
		
		for (iLoopVari = 0; iLoopVari < iDimension; iLoopVari++){
			printf("Enter Collocation Pts for %d Dimension (%d) :\n", iLoopVari + 1, veciNumberofPts[iLoopVari]);
			for (iLoopVarj = 0; iLoopVarj < veciNumberofPts[iLoopVari]; iLoopVarj++){
				scanf("%lf", &vec2dCollocationPts[iLoopVari][iLoopVarj]);
			}
		}
	}

/*****************************************************************
put_AdmissibilityCoefficient - It is used for inputting the Admissibility Coefficient for the Clustering Technique.
					It is only used by the function initialize process.
*****************************************************************/	
	void clsClusteringApproximation::put_AdmissibilityCoefficient(){
		printf("Enter the Admissibility Coefficient ? ");
		scanf("%lf", &dAdmissibilityCoefficient);
	}	

/*****************************************************************
put_Val - It is used for inputting the Vector C for the Clustering Technique.
		It is only used by the function initialize process.
*****************************************************************/	
	void clsClusteringApproximation::put_Val(){
		int iLoopVari;
		
		vecdVal.resize(veciNumberofPts[1]);
		printf("Enter Value for Multiplication (%d) :\n", veciNumberofPts[1]);
		for (iLoopVari = 0; iLoopVari < veciNumberofPts[1]; iLoopVari++){
			scanf("%lf", &vecdVal[iLoopVari]);
		}
	}

/*****************************************************************
multiply - It is used for multiplication of Matrix and Vector.
		The Matrix and the vector should be conformable. The Pointer to Answer Vector is resized and assigned value 0 before multiplication.
Arguments:
pvec2dTemp1	-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to be multiplied.
pvecdTemp2	-	Pointer to Vector Double. It is a Pointer to vector to be multiplied.
pvecdAns		-	Pointer to Vector Double. It is a Pointer to vector to stroe the answer.
*****************************************************************/	
	void clsClusteringApproximation::multiply (vector< vector<double> >* pvec2dTemp1, vector<double>* pvecdTemp2, vector <double>* pvecdAns){
		int iLoopVari, iLoopVarj;
		
		pvecdAns->resize(pvec2dTemp1->size());
		for (iLoopVari = 0; iLoopVari < pvec2dTemp1->size(); iLoopVari++){
			(*pvecdAns)[iLoopVari] = 0;
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dTemp1)[iLoopVari].size(); iLoopVarj++){
				(*pvecdAns)[iLoopVari] += (*pvec2dTemp1)[iLoopVari][iLoopVarj]*(*pvecdTemp2)[iLoopVarj];
			}
		}
	}

/*****************************************************************
multiply_matrix - It is used for multiplication of Matrix and Matrix.
		The Matrix and the Matrix should be conformable. The Pointer to Answer Matrix is resized and assigned value 0 before multiplication.
Arguments:
pvec2dTemp1	-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to be multiplied.
pvec2dTemp2	-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to be multiplied.
pvec2dAns		-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to stroe the answer.
*****************************************************************/	
	void clsClusteringApproximation::multiply_matrix (vector< vector<double> >* pvec2dTemp1, vector< vector<double> >* pvec2dTemp2, vector < vector<double> >* pvec2dAns){
		int iLoopVari, iLoopVarj, iLoopVark;
		
		initialize_matrix(pvec2dAns, pvec2dTemp1->size(), (*pvec2dTemp2)[0].size(), 0);
		for (iLoopVari = 0; iLoopVari < pvec2dAns->size(); iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dAns)[iLoopVari].size(); iLoopVarj++){
				for (iLoopVark = 0; iLoopVark < pvec2dTemp2->size(); iLoopVark++){
					(*pvec2dAns)[iLoopVari][iLoopVarj] += (*pvec2dTemp1)[iLoopVari][iLoopVark]*(*pvec2dTemp2)[iLoopVark][iLoopVarj];
				}
			}
		}
	}

/*****************************************************************
custom_multiply - It is used for multiplication of Matrix and Vector.
		The Matrix and the Matrix should be conformable. The Pointer to Answer Vector is NOT reaized and NOT assigned value 0 before multiplication.
		The Answer of each Row and Column Multiplication of Matrix and Vector is added in a particular location of Vector Answer. 
Arguments:
pvec2dTemp1		-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to be multiplied.
pvecdTemp2		-	Pointer to Vector Double. It is a Pointer to Vector to be multiplied.
pvecdAns			-	Pointer to Vector Double. It is a Pointer to Vector to stroe the answer.
pclsClusterLinearTemp	-	Pointer to Class Cluster Linear . The start and the end of the Cluster Linear define the start and end pointe for adding the anseer in answer Vector.
*****************************************************************/	
	void clsClusteringApproximation::custom_multiply (vector< vector<double> >* pvec2dTemp1, vector<double>* pvecdTemp2, vector <double>* pvecdAns, clsClusterLinear* pclsClusterLinearTemp){
		int iLoopVari, iLoopVarj;
		int iCounter2;
		
		pvecdAns->resize(pvec2dTemp1->size());
		for (iLoopVari = 0; iLoopVari < pvec2dTemp1->size(); iLoopVari++){
			(*pvecdAns)[iLoopVari] = 0;
			iCounter2 = pclsClusterLinearTemp->get_iStart();
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dTemp1)[iLoopVari].size(); iLoopVarj++){
				(*pvecdAns)[iLoopVari] += (*pvec2dTemp1)[iLoopVari][iLoopVarj]*(*pvecdTemp2)[iCounter2];
				iCounter2++;
			}
		}
	}

/*****************************************************************
copy - It is used for copying of Vector to another Vector.
	The new Vector is resized to match the size of the vector being copied.
Arguments:
pvecdTemp	-	Pointer to Vector Double. It is a Pointer to Vector to be copied.
pvecdAns	-	Pointer to Vector Double. It is a Pointer to Vector storing the copy.
*****************************************************************/	
	void clsClusteringApproximation::copy(vector<double>* pvecdTemp, vector<double>* pvecdAns){
		int iLoopVari;
		
		pvecdAns->resize(pvecdTemp->size());
		for (iLoopVari = 0; iLoopVari < pvecdAns->size(); iLoopVari++){
			(*pvecdAns)[iLoopVari] = (*pvecdTemp)[iLoopVari];
		}
	}

/*****************************************************************
add - It is used for adding of Vector and Vector.
		The Vector and the vector should be of same size. The Pointer to Answer Vector is resized.
Arguments:
pvecdTemp1	-	Pointer to Vector Double. It is a Pointer to Vector to be added.
pvecdTemp2	-	Pointer to Vector Double. It is a Pointer to vector to be added.
pvecdAns		-	Pointer to Vector Double. It is a Pointer to vector to store the answer.
*****************************************************************/	
	void clsClusteringApproximation::add(vector<double>* pvecTemp1, vector<double>* pvecTemp2, vector<double>* pvecAns){
		int iLoopVari;
		
		pvecAns->resize(pvecTemp1->size());
		for (iLoopVari = 0; iLoopVari < pvecAns->size(); iLoopVari++){
			(*pvecAns)[iLoopVari] = (*pvecTemp1)[iLoopVari] + (*pvecTemp2)[iLoopVari] ;
		}
	}

/*****************************************************************
transpose - It is used for transposing of Matrix.
		The Pointer to Answer Matrix is resized.
Arguments:
pvec2dTemp	-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to be transpose.
pvec2dAns		-	Pointer to Vector of Vector Double. It is a Pointer to Matrix to store the transpose.
*****************************************************************/	
	void clsClusteringApproximation::transpose(vector<  vector<double> >* pvec2dTemp, vector< vector<double> >* pvec2dAns){
		int iLoopVari, iLoopVarj;
		
		pvec2dAns->resize((*pvec2dTemp)[0].size());
		for (iLoopVari = 0; iLoopVari < pvec2dAns->size(); iLoopVari++){
			(*pvec2dAns)[iLoopVari].resize(pvec2dTemp->size());
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dAns)[iLoopVari].size(); iLoopVarj++){
				(*pvec2dAns)[iLoopVari][iLoopVarj] = (*pvec2dTemp)[iLoopVarj][iLoopVari];
			}
		}
	}

/*****************************************************************
initialize_vector- It is used for initializing vector.
		The Pointer to Vector is resized according to gven value. All the elements are assigned the same value as iInitialValue. 
Arguments:
pvecdTemp		-	Pointer to Vector Double. It is a Pointer to Vector to be initialized.
iSize			-	Integer. It is variable denoting the size of Vector to be initialized.
iInitialValue	-	Integer. It is varaibale denoting the initial Value for the vector.
*****************************************************************/		
	void clsClusteringApproximation::initialize_vector(vector<double>* pvecdTemp, int iSize, int iInitialValue){
		int iLoopVari;
		
		pvecdTemp->resize(iSize);
		for (iLoopVari = 0; iLoopVari < iSize; iLoopVari++){
			(*pvecdTemp)[iLoopVari] = iInitialValue;
		}
	}	

/*****************************************************************
initialize_matrix- It is used for initializing Matrix.
		The Pointer to Matrix is resized according to gven value. All the elements are assigned the same value as iInitialValue. 
Arguments:
pvecdTemp		-	Pointer to Vector Double. It is a Pointer to Vector to be initialized.
iRows		-	Integer. It is variable denoting the Number of Rows of the Matrix to be initialized.
iCols			-	Integer, It is variable denoting the Number of columns of the Matrix to be initialized.
iInitialValue	-	Integer. It is varaibale denoting the initial Value for the Matrix.
*****************************************************************/	
	void clsClusteringApproximation::initialize_matrix(vector< vector<double> >* pvec2dMatrix, int iRows, int iCols, int iInitialValue){
		int iLoopVari, iLoopVarj;
		
		pvec2dMatrix->resize(iRows);
		for (iLoopVari = 0; iLoopVari < (*pvec2dMatrix).size(); iLoopVari++){
			(*pvec2dMatrix)[iLoopVari].resize(iCols);
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dMatrix)[iLoopVari].size(); iLoopVarj++){
				(*pvec2dMatrix)[iLoopVari][iLoopVarj] = iInitialValue;
			}
		}
	}

/*****************************************************************
add_ClusterLinear - It is used for creating and adding new ClusterLinear to vector containing the pointers to Cluster Linear .
Arguments:
iDimension	-	Integer. It is a variable pointing to the dimension in which the Linear Cluster is to be created and added.
iStartTemp	-	Integer. It is a temporary variable denoting the start index of the collocation points in  the Cluster Linear to be created and added to the vector containing pointers to Cluster Linear.
iEndTemp	-	Integer, It is a temporary variable denoting the end index of the collocation points in  the Cluster Linear to be created and added to the vector containing pointers to Cluster Linear.
iLeftTEmp	-	Integer. It is a temporary variable denoting the value of the Right Child of the Cluster Linear to be created and added to the vector containing pointers to Cluster Linear.
iRightTemp	-	Integer. It is a temporary variable denoting the value of the Left Child of the Cluster Linear to be created and added to the vector containing pointers to Cluster Linear.
*****************************************************************/	
	void clsClusteringApproximation::add_ClusterLinear(int iDimension, int iStartTemp, int iEndTemp, int iLeftTemp, int iRightTemp){
		clsClusterLinear* pclsClusterLinearTemp = new clsClusterLinear(iStartTemp, iEndTemp, iLeftTemp, iRightTemp, iDegree);
		vec2pclsClusterLinear[iDimension].push_back(pclsClusterLinearTemp);
	}
	
/*****************************************************************
build_ClusterTree_Linear - It is used for builiding the cluster tree further from a particular level.
					It is a recursive function which calls itself passing arguments for the lower level to be built. It builts the tree bottom up.
Arguments:
iOffset				-	Integer. It is a variable denoting the offset from the start of the collocation points determing the start of the particular Cluster Linear.
iNumberofPtsInClusterLinear	-	Integer. It is a temporary variable denoting the minimum number of points in  the new Cluster Linear to be created.
						If this number is reached then the cluster is created as leaf. Otherwise as a node.
piDimension			-	Pointer to Integer. It is a pointer to a variable denoting the Dimension of the Cluster Linear  Tree being built.
*****************************************************************/		
	void clsClusteringApproximation::build_ClusterTree_Linear(int iOffset, int iNumberofPtsInClusterLinear, int* piDimension){
		int iTempLeftChild, iTempRightChild, iStart, iEnd;
		int iStartLeft, iPtsInLeft, iStartRight, iPtsInRight;
		
		iStart = iOffset;
		iEnd = iOffset + iNumberofPtsInClusterLinear - 1;
		if (iNumberofPtsInClusterLinear == 1){
			add_ClusterLinear(*piDimension, iStart, iEnd, -1, -1);
		}
		else{
			iStartLeft = iOffset;
			iPtsInLeft =  iNumberofPtsInClusterLinear - (iNumberofPtsInClusterLinear/2);
			build_ClusterTree_Linear(iStartLeft, iPtsInLeft, piDimension);
			iTempLeftChild = vec2pclsClusterLinear[*piDimension].size() - 1;
			
			iStartRight = iOffset + iPtsInLeft;
			iPtsInRight = iNumberofPtsInClusterLinear/2;
			build_ClusterTree_Linear(iStartRight, iPtsInRight, piDimension);
			iTempRightChild = vec2pclsClusterLinear[*piDimension].size() - 1;
            
			add_ClusterLinear(*piDimension, iStart, iEnd, iTempLeftChild, iTempRightChild);
		}
	}

/*****************************************************************
add_ClusterPairNear - It is used for creating and adding new Cluster Pair Near to vector containing the pointers to Cluster Pair Near .
Arguments:
pclsClusterLinear1Temp	-	Pointer to Class Cluster Linear. It is a pointer to class Cluster Linear denoting the Cluster Linear in 1st dimension forming the Cluster Pair Near, which is to be created and added.
pclsClusterLinear2Temp	-	Pointer to Class Cluster Linear. It is a pointer to class Cluster Linear denoting the Cluster Linear in 2nd dimension forming the Cluster Pair Near, which is to be created and added.
*****************************************************************
Return Type	-	Vector of Vector Double.
Returns		-	The pointer to vec2dFnValAtChebyNodeInCluster variable of class ClusterPair
*****************************************************************/		
	void clsClusteringApproximation::add_ClusterPairNear(clsClusterLinear* pclsClusterLinear1Temp, clsClusterLinear* pclsClusterLinear2Temp){
		int iRows, iCols;
		
		iRows = (pclsClusterLinear1Temp->get_iEnd() - pclsClusterLinear1Temp->get_iStart()) + 1;
		iCols = (pclsClusterLinear2Temp->get_iEnd() - pclsClusterLinear2Temp->get_iStart()) + 1;
		clsClusterPair* pclsClusterPairTemp = new clsClusterPair(iDimension, pclsClusterLinear1Temp, pclsClusterLinear2Temp, iRows, iCols);
		vecpclsClusterPairNear.push_back(pclsClusterPairTemp);
	}

/*****************************************************************
add_ClusterPairFar - It is used for creating and adding new Cluster Pair Near to vector containing the pointers to Cluster Pair Far .
Arguments:
pclsClusterLinear1Temp	-	Pointer to Class Cluster Linear. It is a pointer to class Cluster Linear denoting the Cluster Linear in 1st dimension forming the Cluster Pair Far, which is to be created and added.
pclsClusterLinear2Temp	-	Pointer to Class Cluster Linear. It is a pointer to class Cluster Linear denoting the Cluster Linear in 2nd dimension forming the Cluster Pair Far, which is to be created and added.
******************************************************************
Return Type	-	Vector of Vector Double.
Returns		-	The pointer to vec2dFnValAtChebyNodeInCluster variable of class ClusterPair
*****************************************************************/	
	void clsClusteringApproximation::add_ClusterPairFar(clsClusterLinear* pclsClusterLinear1Temp, clsClusterLinear* pclsClusterLinear2Temp){
		int iRows, iCols;
		
		iRows = iDegree + 1;
		iCols = iDegree + 1;
		clsClusterPair* pclsClusterPairTemp = new clsClusterPair(iDimension, pclsClusterLinear1Temp, pclsClusterLinear2Temp, iRows, iCols);
		vecpclsClusterPairFar.push_back(pclsClusterPairTemp);
	}

/*****************************************************************
leaf - It is used for testing whether a ClusterLinear is Leaf or not.
Arguments:
iDimension	-	Integer. It is a temporary variable denoting the dimension of the ClusterLinear.
iIndex	-	Integer. It is a temporary variable denoting the index of ClusterLinear.
******************************************************************
Return Type	-	Boolean.
Returns		-	True if ClusterLinear is Leaf otherwise False
*****************************************************************/	
	bool clsClusteringApproximation::leaf(int iDimension, int iIndex){
		if (vec2pclsClusterLinear[iDimension][iIndex]->get_iLeft() == -1)
			return (true);
		else
			return (false);
	}

/*****************************************************************
lagrange_value_for_Cluster - It is used for finding Lagrange Value for Cluster Pair Far.
Arguments:
pclsClusterLinearTemp	-	Pointer to Class Cluster Linear. It is a pointer to Class ClusterLinear.
dClusterStart		-	Double. It is a temporary variable denoting the start of the ClusterLinear.
dStep			-	Double. It is a temporary variable denoting the step for the Chebyshev nodes for the ClusterLinear
iDimension			-	Integer. It is used to denote the Dimension of the Cluster Linear being worked on.
*****************************************************************/
	void clsClusteringApproximation::lagrange_value_for_Cluster(clsClusterLinear* pclsClusterLinearTemp, double dClusterStart, double dClusterEnd, int iDimension){
		int iLoopVari, iLoopVarj, iLoopVark;
		double dNode, dChebyNode;
		iCounter = 0;
		for (iLoopVari = pclsClusterLinearTemp->get_iStart(); iLoopVari <= pclsClusterLinearTemp->get_iEnd(); iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < iDegree + 1; iLoopVarj++){
				(*(pclsClusterLinearTemp->get_ptr_LagValAtCollocPtsInCluster()))[iCounter][iLoopVarj] = 1;
				dNode = dClusterStart + (0.5*(dClusterEnd - dClusterStart)*(cos(((2*iLoopVarj + 1.0)/(2*(iDegree + 1.0)))*PI) + 1));
				for (iLoopVark = 0; iLoopVark < iDegree + 1; iLoopVark++){
					dChebyNode = dClusterStart + (0.5*(dClusterEnd - dClusterStart)*(cos(((2*iLoopVark + 1.0)/(2*(iDegree + 1.0)))*PI) + 1));
					if (iLoopVarj != iLoopVark){
						(*(pclsClusterLinearTemp->get_ptr_LagValAtCollocPtsInCluster()))[iCounter][iLoopVarj] *= (vec2dCollocationPts[iDimension][iLoopVari] - dChebyNode)/(dNode - dChebyNode);
                    }
				}
			}
			iCounter++;
		}
	}
	
/*****************************************************************
fill_vec2dFnValAtChebyTemp_for_Near - It is used for filling the vec2dFnValAtChebyNodeInCluster for Near ClusterPair.
Arguments:
pvec2dFnValAtChebyTemp	-	Pointer to Vector of Vector of Double. It is the pointer to Matrix which is to be filled.
pclsClusterLinear1		-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 1st Dimension forming the ClusterPair. 
pclsClusterLinear2		-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 2nd Dimension forming the ClusterPair.
*****************************************************************/	
	void clsClusteringApproximation::fill_FnValAtChebyNode_for_Near(vector< vector<double> >* pvec2dFnValAtChebyTemp, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2){
		int iLoopVari, iLoopVarj, iCounter;
		
		iCounter = 0;
		for (iLoopVari = pclsClusterLinear1->get_iStart(); iLoopVari <= pclsClusterLinear1->get_iEnd(); iLoopVari++){
			for (iLoopVarj = pclsClusterLinear2->get_iStart(); iLoopVarj <= pclsClusterLinear2->get_iEnd(); iLoopVarj++){
				(*pvec2dFnValAtChebyTemp)[iCounter].push_back(kernel_function(vec2dCollocationPts[0][iLoopVari], vec2dCollocationPts[1][iLoopVarj]));
			}
			iCounter++;
		}
	}

/*****************************************************************
fill_vec2dFnValAtChebyTemp_for_Far - It is used for filling the vec2dFnValAtChebyNodeInCluster for Far ClusterPair.
Arguments:
pvec2dFnValAtChebyTemp	-	Pointer to Vector of Vector of Double. It is the pointer to Matrix which is to be filled.
pclsClusterLinear1		-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 1st Dimension forming the ClusterPair. 
pclsClusterLinear2		-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 2nd Dimension forming the ClusterPair.
*****************************************************************/	
	void clsClusteringApproximation::fill_FnValAtChebyNode_for_Far(vector< vector<double> >* pvec2dFnValAtChebyTemp, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2){
		int iLoopVari, iLoopVarj;
		double dChebyshevX, dChebyshevY;
		double dClusterStart1,dClusterStart2, dClusterEnd1, dClusterEnd2;
		
		dClusterStart1 = vec2dCollocationPts[0][pclsClusterLinear1->get_iStart()];
		dClusterStart2 = vec2dCollocationPts[1][pclsClusterLinear2->get_iStart()];
		dClusterEnd1 = vec2dCollocationPts[0][pclsClusterLinear1->get_iEnd()];
		dClusterEnd2 = vec2dCollocationPts[1][pclsClusterLinear2->get_iEnd()];
		
		for (iLoopVari = 0; iLoopVari < iDegree + 1; iLoopVari++){
			dChebyshevX = dClusterStart1 + (0.5*(dClusterEnd1 - dClusterStart1)*(cos(((2*iLoopVari + 1.0)/(2*(iDegree + 1.0)))*PI) + 1));
			for (iLoopVarj = 0; iLoopVarj < iDegree + 1; iLoopVarj++){
				dChebyshevY = dClusterStart2 + (0.5*(dClusterEnd2 - dClusterStart2)*(cos(((2*iLoopVarj + 1.0)/(2*(iDegree + 1.0)))*PI) + 1));
				(*pvec2dFnValAtChebyTemp)[iLoopVari].push_back(kernel_function(dChebyshevX, dChebyshevY));
			}
		}
		
	}
	
/*****************************************************************
check_and_fill_LagVal - It is used for checking if Lagrange Values for the Cluster Linear have been calculated or not.
				In case they have not been calculated, they are calculated.
Arguments:
pclsClusterLinear	-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear which has to checked for Lagrange Values. 
iDimension		-	Integer. The dimension of Cluster Linear to be worked on.
*****************************************************************/	
	void clsClusteringApproximation::check_and_fill_LagVal(clsClusterLinear* pclsClusterLinear, int iDimension){
		double dClusterStart, dClusterEnd;
		
		dClusterStart = vec2dCollocationPts[iDimension][pclsClusterLinear->get_iStart()];
		dClusterEnd = vec2dCollocationPts[iDimension][pclsClusterLinear->get_iEnd()];
		if ((pclsClusterLinear->get_ptr_LagValAtCollocPtsInCluster())->size() == 0){
			pclsClusterLinear->increase_size_LagValAtCollocPtsInCluster(iDegree);
			lagrange_value_for_Cluster(pclsClusterLinear, dClusterStart, dClusterEnd, iDimension);
			if (iDimension == 1){
				transpose(pclsClusterLinear->get_ptr_LagValAtCollocPtsInCluster() ,&vec2dTempForTranspose);
				custom_multiply(&vec2dTempForTranspose, &vecdVal,pclsClusterLinear->get_ptr_vecdW(), pclsClusterLinear);	
			}
		}
	}

/*****************************************************************
tree_traverse - It is used for travelling deeper into the Cluster Pair Tree to build Cluster Pairs.
Arguments:
pclsClusterLinear1	-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 1st dimension whose children are to be checked for making Cluster Pair. 
pclsClusterLinear2	-	Pointer to Class Cluster Linear. It is used to define the ClusterLinear in 2nd dimension whose children are to be checked for making Cluster Pair. 
*****************************************************************/	
	void clsClusteringApproximation::tree_traverse(clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2){
		int iVari, iVarj;
		iVari = pclsClusterLinear1->get_iLeft();
		 while (iVari != -1){
			iVarj = pclsClusterLinear2->get_iLeft();
			while (iVarj != -1){
				build_ClusterPair(iVari, iVarj);
				if (iVarj == pclsClusterLinear2->get_iLeft())
					iVarj = pclsClusterLinear2->get_iRight();
				else
					iVarj = -1;
			}
			if (iVari == pclsClusterLinear1->get_iLeft())
				iVari = pclsClusterLinear1->get_iRight();
			else
				iVari = -1;
		}	
	}
	
/*****************************************************************
build_ClusterPair - It is used for builiding the Cluster Pair Tree further from a particular level.
					It is a recursive function which calls itself passing arguments for the lower level to be built. It builts the tree top down.
Arguments:
iIndex1	-	Integer. It is used to define the index of the ClusterLinear in the 1st Dimension. 
iIndex2	-	Integer. It is used to define the index of the ClusterLinear in the 2st Dimension.
*****************************************************************/	
	void clsClusteringApproximation::build_ClusterPair (int iIndex1, int iIndex2){
		
		clsClusterLinear *pclsClusterLinear1, *pclsClusterLinear2;	
		
		pclsClusterLinear1 = vec2pclsClusterLinear[0][iIndex1];
		pclsClusterLinear2 = vec2pclsClusterLinear[1][iIndex2];
		if (leaf(0,iIndex1) || leaf(1,iIndex2)){
			add_ClusterPairNear(pclsClusterLinear1, pclsClusterLinear2);
			(pclsClusterLinear1->get_ptr_AppearsIn())->push_back(-vecpclsClusterPairNear.size());
			(pclsClusterLinear2->get_ptr_AppearsIn())->push_back(-vecpclsClusterPairNear.size());
		}
		else if (admissible_ClusterPair(iIndex1, iIndex2)){
			add_ClusterPairFar(pclsClusterLinear1, pclsClusterLinear2);
			(pclsClusterLinear1->get_ptr_AppearsIn())->push_back(vecpclsClusterPairFar.size());
			(pclsClusterLinear2->get_ptr_AppearsIn())->push_back(vecpclsClusterPairFar.size());
		}
		else{
			tree_traverse(pclsClusterLinear1, pclsClusterLinear2);
		}
	}

/*****************************************************************
compute_S - It is used for computing S for each ClusterLinear in 2nd Dimension.
Arguments:
pbAtleastOneFar	-	Pointer to Boolean. It is used to indicate if the ClusterLinear of 1st Dimension appears in ClusterPair otleast once or not. 
iPairNumber	-	Integer. It is used to indicate the index of ClusterPair Near in a vector containing all the cluster pairs near .
pvecdS		-	Pointer to Vector Double. It is used to indicate the value of Vector S for ClusterLinear of 1st Dimension.
*****************************************************************/
	void clsClusteringApproximation::compute_S(bool* pbAtLeastOneFar, int iPairNumber, vector<double>* pvecdS){
		vector<double> vecdTemp;
		vector<double> vecdSTemp;
		clsClusterLinear* pclsClusterLinearTemp;
				
		*pbAtLeastOneFar = true;
		pvec2dFnValAtChebyTemp = vecpclsClusterPairFar[iPairNumber]->get_ptr_FnValAtChebyNodeInCluster();
		pclsClusterLinearTemp = vecpclsClusterPairFar[iPairNumber]->get_ptr_ClusterLinear(1);
		multiply(pvec2dFnValAtChebyTemp, pclsClusterLinearTemp->get_ptr_vecdW(), &vecdTemp);
		copy(pvecdS, &vecdSTemp);
		add(&vecdTemp, &vecdSTemp, pvecdS);		
	}
	
/*****************************************************************
compute_final_answer_for_near - It is used for computing the Final Answer for ClusterPair Near.
Arguments:
iIndex		-	Integer. It is used to indicate the index of ClusterLinear in a vector containing all the ClusterLinear.
iPairNumber	-	Integer. It is used to indicate the index of ClusterPair Near in a vector containing all the cluster pairs near .
*****************************************************************/	
	void clsClusteringApproximation::compute_final_answer_for_near(int iIndex, int iPairNumber){
		int iCounter, iLoopVark;
		vector<double> vecdTemp;
		clsClusterLinear* pclsClusterLinearTemp;		
	
		pclsClusterLinearTemp = vecpclsClusterPairNear[iPairNumber]->get_ptr_ClusterLinear(1);
		pvec2dFnValAtChebyTemp = vecpclsClusterPairNear[iPairNumber]->get_ptr_FnValAtChebyNodeInCluster();
		custom_multiply(pvec2dFnValAtChebyTemp, &vecdVal, &vecdTemp, pclsClusterLinearTemp);
		iCounter = 0;
		for (iLoopVark = vec2pclsClusterLinear[0][iIndex]->get_iStart(); iLoopVark <= vec2pclsClusterLinear[0][iIndex]->get_iEnd(); iLoopVark++){
			vecdFinalAnswer[iLoopVark] += vecdTemp[iCounter];
			iCounter++;
		}	
	}

/*****************************************************************
compute_final_answer_for_far - It is used for computing the Final Answer for ClusterLinear in 1st Dimension.
Arguments:
iIndex	-	Integer. It is used to indicate the index of Cluster Linear Near in a vector containing all the cluster pairs near .
pvecdS	-	Pointer to Vector Double. It is used to indicate the value of Vector S for ClusterLinear of 1st Dimension.
*****************************************************************/	
	void clsClusteringApproximation::compute_final_answer_for_far(int iIndex, vector<double>* pvecdS){
		int iLoopVarj, iCounter;
		vector<double> vecdTemp;
		
		multiply(vec2pclsClusterLinear[0][iIndex]->get_ptr_LagValAtCollocPtsInCluster(), pvecdS, &vecdTemp);
		iCounter = 0;
		for (iLoopVarj = vec2pclsClusterLinear[0][iIndex]->get_iStart(); iLoopVarj <= vec2pclsClusterLinear[0][iIndex]->get_iEnd(); iLoopVarj++){
			vecdFinalAnswer[iLoopVarj] += vecdTemp[iCounter];
			iCounter++;
		}	
	}
	
/*****************************************************************
finalAnswer - It is used for find the Final Answer after building the Cluster Pairs.
*****************************************************************/	
	void clsClusteringApproximation::finalAnswer(){
		int iLoopVari, iLoopVarj, iPairNumber;
		bool bAtLeastOneFar;
		vector<int>* pveciAppearsInTemp;
		vector<double> vecdS;
		printf("\nCalculating Final Answer..............");
		for (iLoopVari = 0; iLoopVari < vec2pclsClusterLinear[0].size(); iLoopVari++){
			pveciAppearsInTemp = vec2pclsClusterLinear[0][iLoopVari]->get_ptr_AppearsIn();
			if (pveciAppearsInTemp->size() <= 0){
				continue;
			}
			clk3Start = clock();
			initialize_vector(&vecdS, iDegree + 1, 0);
			bAtLeastOneFar = false;
			for (iLoopVarj = 0; iLoopVarj < pveciAppearsInTemp->size(); iLoopVarj++){
				iPairNumber = abs((*pveciAppearsInTemp)[iLoopVarj]) - 1;
				if ((*pveciAppearsInTemp)[iLoopVarj] > 0){
					compute_S(&bAtLeastOneFar, iPairNumber, &vecdS);
					clk3Stop = clock();
					clk3Diff += clk3Stop - clk3Start;
				}
				else{
					clk4Start = clock();
					compute_final_answer_for_near(iLoopVari, iPairNumber);
					clk4Stop = clock();
					clk4Diff += clk4Stop - clk4Start;
				}
			}
			if (bAtLeastOneFar == true){
				clk3Start = clock();
				compute_final_answer_for_far(iLoopVari, &vecdS);
				clk3Stop = clock();
				clk3Diff += clk3Stop - clk3Start;				
			}
		}
	}

/*****************************************************************
make_matrix - It is used for making matrix for calculation of Frobenius Error.
*****************************************************************/		
	void clsClusteringApproximation::make_matrix(){
		int iLoopVari, iLoopVarj, iPairNumber;
		vector<int>* pveciAppearsInTemp;
		clsClusterLinear* pclsClusterLinearTemp;
		
		printf("\nMaking Matrix..............");
		for (iLoopVari = 0; iLoopVari < vec2pclsClusterLinear[0].size(); iLoopVari++){
			pveciAppearsInTemp = vec2pclsClusterLinear[0][iLoopVari]->get_ptr_AppearsIn();
			if (pveciAppearsInTemp->size() <= 0){
				continue;
			}
			for (iLoopVarj = 0; iLoopVarj < pveciAppearsInTemp->size(); iLoopVarj++){
				iPairNumber = abs((*pveciAppearsInTemp)[iLoopVarj]) - 1;
				if ((*pveciAppearsInTemp)[iLoopVarj] > 0){
					pvec2dFnValAtChebyTemp = vecpclsClusterPairFar[iPairNumber]->get_ptr_FnValAtChebyNodeInCluster();
					pclsClusterLinearTemp = vecpclsClusterPairFar[iPairNumber]->get_ptr_ClusterLinear(1);
					make_matrix_for_far(pvec2dFnValAtChebyTemp, vec2pclsClusterLinear[0][iLoopVari], pclsClusterLinearTemp);
				}
				else{
					pvec2dFnValAtChebyTemp = vecpclsClusterPairNear[iPairNumber]->get_ptr_FnValAtChebyNodeInCluster();
					pclsClusterLinearTemp = vecpclsClusterPairNear[iPairNumber]->get_ptr_ClusterLinear(1);
					make_matrix_for_near(pvec2dFnValAtChebyTemp, vec2pclsClusterLinear[0][iLoopVari], pclsClusterLinearTemp);
				}
			}
		}
	}

/*****************************************************************
print_ClusterLinear - It is used for printing a particular ClusterLinear on screen. Used for debugging Purposes.
Arguments:
iDimension		-	Integer. It is used to denote the dimension of the ClusterLinear to be printed.
iIndexNumber 	-	Integer, It is used to denote the index number of ClusterLinear to be printed in the vector conating all the ClusterLinear.
*****************************************************************/		
	void clsClusteringApproximation::print_ClusterLinear(int iDimension, int iIndexNumber){
		vec2pclsClusterLinear[iDimension][iIndexNumber]->print(iIndexNumber);
	}

/*****************************************************************
print_ClusterLinear - It is used for printing a particular ClusterLinear on screen. Used for debugging Purposes.
*****************************************************************/	
	void clsClusteringApproximation::print_ClusterLinear(){
		int iLoopVari, iLoopVarj;
		
		for (iLoopVari = 0; iLoopVari < iDimension; iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < vec2pclsClusterLinear[iLoopVari].size(); iLoopVarj++){
				vec2pclsClusterLinear[iLoopVari][iLoopVarj]->print(iLoopVarj);
			}
		}
	}

/*****************************************************************
print_ClusterPairFar - It is used for printing a particular ClusterPair Far on screen. Used for debugging Purposes.
Arguments:
iIndexNumber 	-	Integer, It is used to denote the index number of ClusterPair Far to be printed in the vector conating all the ClusterPair Far.
*****************************************************************/	
	void clsClusteringApproximation::print_ClusterPairFar(int iIndexNumber){
		vecpclsClusterPairFar[iIndexNumber]->print(iIndexNumber, iDegree, iDimension);
	}

/*****************************************************************
print_ClusterPairNear - It is used for printing a particular ClusterPair Near on screen. Used for debugging Purposes.
Arguments:
iIndexNumber 	-	Integer, It is used to denote the index number of ClusterPair Near to be printed in the vector conating all the ClusterPair Near.
*****************************************************************/		
	void clsClusteringApproximation::print_ClusterPairNear(int iIndexNumber){
		vecpclsClusterPairNear[iIndexNumber]->print(iIndexNumber, iDegree, iDimension);
	}	

/*****************************************************************
print_ClusterPairFar - It is used for printing all ClusterPair Far on screen. Used for debugging Purposes.
*****************************************************************/	
	void clsClusteringApproximation::print_ClusterPairFar(){
		for (int iLoopVari = 0; iLoopVari < vecpclsClusterPairFar.size(); iLoopVari++){
			vecpclsClusterPairFar[iLoopVari]->print(iLoopVari, iDegree, iDimension);
		}	
	}

/*****************************************************************
print_ClusterPairNear - It is used for printing all ClusterPair Near on screen. Used for debugging Purposes.
*****************************************************************/		
	void clsClusteringApproximation::print_ClusterPairNear(){
		for (int iLoopVari = 0; iLoopVari < vecpclsClusterPairNear.size(); iLoopVari++){
			vecpclsClusterPairNear[iLoopVari]->print(iLoopVari, iDegree, iDimension);
		}	
	}

/*****************************************************************
print - It is used for printing a vector on screen. Used for debugging Purposes.
Arguments:
pchName 	-	Pointer to Character, It is used to denote the name of the vector to be printed.
pveciPrint	-	Pointer to Vector Integer. It is a pointer to vector to be printed.
*****************************************************************/		
	void clsClusteringApproximation::print(char* pchName, vector<int>* pveciPrint){
		int iLoopVari;
		
		printf("%s", pchName);
		for (iLoopVari = 0; iLoopVari < pveciPrint->size(); iLoopVari++){
			printf("\t%d", (*pveciPrint)[iLoopVari]);
		}
		printf("\n");
	}

/*****************************************************************
print - It is used for printing a vector on screen. Used for debugging Purposes.
Arguments:
pchName 	-	Pointer to Character, It is used to denote the name of the vector to be printed.
pvecdPrint	-	Pointer to Vector Double. It is a pointer to vector to be printed.
*****************************************************************/	
	void clsClusteringApproximation::print(char* pchName, vector<double>* pvecdPrint){
		int iLoopVari;
		
		printf("%s", pchName);
		for (iLoopVari = 0; iLoopVari < pvecdPrint->size(); iLoopVari++){
			printf("\t%6.1f", (*pvecdPrint)[iLoopVari]);
		}
		printf("\n");
	}

/*****************************************************************
print - It is used for printing matrix on screen. Used for debugging Purposes.
Arguments:
pchName 		-	Pointer to Character, It is used to denote the name of the vector to be printed.
pvec2dPrint	-	Pointer to Vector of Vector Double. It is a pointer to vector to be printed.
*****************************************************************/	
	void clsClusteringApproximation::print(char* pchName, vector< vector<double> >* pvec2dPrint){
		int iLoopVari, iLoopVarj;
		
		printf("%s", pchName);
		for (iLoopVari = 0; iLoopVari < pvec2dPrint->size(); iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < (*pvec2dPrint)[iLoopVari].size(); iLoopVarj++){
				printf("\t%4.3lf", (*pvec2dPrint)[iLoopVari][iLoopVarj]);
			}
			printf("\n");
		}	
	}

/*****************************************************************
clear_pointers - It is used to free memory from pointers.
*****************************************************************/	
	void clsClusteringApproximation::clear_pointers(){
		int iLoopVari, iLoopVarj;
		
		printf("Cleaning Pointers.......");
		for (iLoopVari = 0; iLoopVari < vec2pclsClusterLinear.size(); iLoopVari++){
			for (iLoopVarj = 0; iLoopVarj < vec2pclsClusterLinear[iLoopVari].size(); iLoopVarj++){
				delete vec2pclsClusterLinear[iLoopVari][iLoopVarj];
			}
			vec2pclsClusterLinear[iLoopVari].resize(0);
		}
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairFar.size(); iLoopVari++){
			delete vecpclsClusterPairFar[iLoopVari];
		}
		vecpclsClusterPairFar.resize(0);
		
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairNear.size(); iLoopVari++){
			delete vecpclsClusterPairNear[iLoopVari];
		}
		vecpclsClusterPairNear.resize(0);
		printf("done\n");
	}

/*****************************************************************
initialize_process - It is used for initialization process. It asks for differen values frm the user and initializes matrix and vectotrs. It is in case when all data is to be inputed by user.
*****************************************************************/	
	void clsClusteringApproximation::initialize_process(){
		put_NumberofPts();
		put_Degree();
		put_AdmissibilityCoefficient();
		put_CollocationPts();
		put_Val();
		initialize_vector(&vecdFinalAnswerDirect, veciNumberofPts[0], 0);
		initialize_matrix(&vec2dMatrix, veciNumberofPts[0], veciNumberofPts[1], 0);		
	}

/*****************************************************************
preprocess - This method performs the preprocessing for clustering approximation technique.
*****************************************************************/		
	void clsClusteringApproximation::preprocess(){
		int iLoopVari, iLoopVarj;
		for (iLoopVari = 0; iLoopVari < iDimension; iLoopVari++){
			printf("\nStarting Building Cluster Trees for %d dimension..............", iLoopVari + 1);
			build_ClusterTree_Linear(0, veciNumberofPts[iLoopVari], &iLoopVari);
		}
		printf("\nStarting Building Cluster Pairs..............");
		build_ClusterPair(vec2pclsClusterLinear[0].size() - 1, vec2pclsClusterLinear[1].size() - 1);
    }

/*****************************************************************
preprocess - This method performs the preprocessing for clustering approximation technique.
*****************************************************************/		
	void clsClusteringApproximation::preprocess_compute(){
		int iLoopVari;
		clsClusterLinear *pclsClusterLinear1, *pclsClusterLinear2;	
		
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairNear.size(); iLoopVari++){
			pclsClusterLinear1 = vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(0);
			pclsClusterLinear2 = vecpclsClusterPairNear[iLoopVari]->get_ptr_ClusterLinear(1);
			pvec2dFnValAtChebyTemp = vecpclsClusterPairNear[iLoopVari]->get_ptr_FnValAtChebyNodeInCluster();
			clk2Start = clock();
			fill_FnValAtChebyNode_for_Near(pvec2dFnValAtChebyTemp, pclsClusterLinear1, pclsClusterLinear2);
			clk2Stop = clock();
			clk2Diff += clk2Stop - clk2Start;			
		}
		for (iLoopVari = 0; iLoopVari < vecpclsClusterPairFar.size(); iLoopVari++){
			pclsClusterLinear1 = vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(0);
			pclsClusterLinear2 = vecpclsClusterPairFar[iLoopVari]->get_ptr_ClusterLinear(1);
			pvec2dFnValAtChebyTemp = vecpclsClusterPairFar[iLoopVari]->get_ptr_FnValAtChebyNodeInCluster();
			clk2Start = clock();
			fill_FnValAtChebyNode_for_Far(pvec2dFnValAtChebyTemp, pclsClusterLinear1, pclsClusterLinear2);
			clk2Stop = clock();
			clk2Diff += clk2Stop - clk2Start;			
			check_and_fill_LagVal(pclsClusterLinear1, 0);
			clk1Start = clock();
			check_and_fill_LagVal(pclsClusterLinear2, 1);
			clk1Stop = clock();
			clk1Diff += (clk1Stop - clk1Start);
		}

	}
	
/*****************************************************************
process - This method performs the clustering approximation technique.
*****************************************************************/		
	void clsClusteringApproximation::process(){
		int iLoopVari, iLoopVarj;
		
		printf("Start..............");
		printf("\n");
		print("Number of Pts: ", &veciNumberofPts);
		printf("Degree: %d\n", iDegree);
		printf("Admissibility Coefficient: %lf\n", dAdmissibilityCoefficient);
		preprocess();
		preprocess_compute();
		finalAnswer();
		printf("\n");
	}

/*****************************************************************
count - It is used for counting the max occurences of ClusterLinear in ClusterPairs Far and Near.
Arguments:
piMaxFar 			-	Pointer to Integer. It is used to denote the name of the vector to be printed.
piMaxNear			-	Pointer to Integer. It is a pointer to vector to be printed.
pvecpclsClusterLinear		-	Pointer to vector of Pointers of ClusterLinear. The vector containg all the pointer to cluster Linear of particular dimension.
piNumberofClusterLinearInFar	-	Pointer to Integer. It is a pointer to a variable which stores the Number of Cluster Linear which were used in forming the ClusterPairFar.
piNumberofClusterLinearInNear	-	Pointer to Integer. It is a pointer to a variable which stores the Number of Cluster Linear which were used in forming the ClusterPairNear.
*****************************************************************/
	void clsClusteringApproximation::count(int* piMaxFar, int* piMaxNear, vector<clsClusterLinear*>* pvecpclsClusterLinear, int* piNumberofClusterLinearInFar, int* piNumberofClusterLinearInNear){
		int iLoopVari, iLoopVarj;
		int iFar, iNear;
		clsClusterLinear* pclsClusterLinear;
		vector<int>* pveciAppearsIn;
		bool bAtleastOneFar;
		bool bAtleastOneNear;
		*piMaxFar = 0;
		*piMaxNear = 0;
		*piNumberofClusterLinearInNear = 0;
		*piNumberofClusterLinearInFar = 0;

		for (iLoopVari = 0; iLoopVari < pvecpclsClusterLinear->size(); iLoopVari++){
			pclsClusterLinear = (*pvecpclsClusterLinear)[iLoopVari];
			pveciAppearsIn = pclsClusterLinear->get_ptr_AppearsIn();
			iFar = 0;
			iNear = 0;
			bAtleastOneFar = false;
			bAtleastOneNear = false;
			for (iLoopVarj = 0; iLoopVarj < pveciAppearsIn->size(); iLoopVarj++){
				if ((*pveciAppearsIn)[iLoopVarj] < 0){
					iNear++;
					bAtleastOneNear = true;
				}
				else{
					iFar++;
					bAtleastOneFar = true;
				}
			}
			if (bAtleastOneNear == true)
				(*piNumberofClusterLinearInNear)++;
			 if (bAtleastOneFar == true)
				(*piNumberofClusterLinearInFar)++;
			if (iNear > *piMaxNear)
				*piMaxNear = iNear;
			if (iFar > *piMaxFar)
				*piMaxFar = iFar;
		}		
}

/*****************************************************************
count - It is used for counting the max occurences of ClusterLinear in ClusterPairs Far and Near.
	It calculates Maximum from both dimension, and stores the Maximum from both dimension.
*****************************************************************/
	void clsClusteringApproximation::count(){
		int iMaxFar1, iMaxNear1;
		int iMaxFar2, iMaxNear2;

		count(&iMaxFar1, &iMaxNear1, &vec2pclsClusterLinear[0], &iNumberofClusterLinearInFar1, &iNumberofClusterLinearInNear1);
		count(&iMaxFar2, &iMaxNear2, &vec2pclsClusterLinear[1], &iNumberofClusterLinearInFar2, &iNumberofClusterLinearInNear2);
		if (iMaxFar1 > iMaxFar2)
			iMaxFar = iMaxFar1;
		else
			iMaxFar = iMaxFar2;
		if (iMaxNear1 > iMaxNear2)
			iMaxNear = iMaxNear1;
		else
			iMaxNear = iMaxNear2;
	}

/*****************************************************************
calculate_answer_direct - It is used for calculating the answer in the direct way. It is used for error calculation.
*****************************************************************/	
	double clsClusteringApproximation::calculate_answer_direct(){
		   return (calculate_answer_direct(&vecdFinalAnswerDirect));
    }       

/*****************************************************************
calculate_answer_direct - It is used for calculating the answer in the direct way. It is used for error calculation.
Arguments:
pvecdFinalAnswerDirect	-	Pointer to Vector Double. This is used for storing the Final Answer by Direct Calculation.
*****************************************************************/	
	double clsClusteringApproximation::calculate_answer_direct(vector<double>* pvecdFinalAnswerDirect){
		int iLoopVari, iLoopVarj;
		double dTemp1, dTempAnswer;
		double dErrorFrobenius;
		double dError;
		double dNormTrue;
		dErrorFrobenius = 0;
		dNormTrue = 0;
		
		for (iLoopVari = 0; iLoopVari <veciNumberofPts[0]; iLoopVari++){
			dTempAnswer = 0;
			for (iLoopVarj = 0; iLoopVarj < veciNumberofPts[1]; iLoopVarj++){
				dTemp1 = kernel_function(vec2dCollocationPts[0][iLoopVari], vec2dCollocationPts[1][iLoopVarj]);
				dTempAnswer += dTemp1*vecdVal[iLoopVarj];
				dError = fabs(vec2dMatrix[iLoopVari][iLoopVarj] - dTemp1);
				vec2dMatrixTemp[iLoopVari][iLoopVarj] = dTemp1;
				dErrorFrobenius += pow(dError,2);
				dNormTrue += pow(dTemp1,2); 
			}
			(*pvecdFinalAnswerDirect)[iLoopVari] = dTempAnswer;
		}
		dErrorFrobenius = fabs(dErrorFrobenius)/(veciNumberofPts[0]*veciNumberofPts[1]);
        return(dErrorFrobenius);		
	}

/*****************************************************************
make_matrix_for_far - It is used for making Matrix for Cluster Pair Far. It is used for error calculation.
Arguments:
pvec2dX			-	Pointer to Vector of Vector Double. This is a temporary variable denoting the Function value at chebyshev nodes in the Cluster Pair Far.
pclsClusterLinear1	-	Pointer to Class ClusterLinear. This is a temporary variable denoting the ClusterLinear in 1st dimension forming the Cluser Pair Far.
pclsClusterLinear2	-	Pointer to Class ClusterLinear. This is a temporary variable denoting the ClusterLinear in 2nd dimension forming the Cluser Pair Far.
*****************************************************************/	
	void clsClusteringApproximation::make_matrix_for_far(vector< vector <double> >* pvec2dX, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2){
		vector< vector<double> > vec2dAns;
		vector< vector<double> > vec2dTemp;
		vector< vector<double> >* pvec2dV1;
		vector< vector<double> >* pvec2dV2;
		
		pvec2dV1 = pclsClusterLinear1->get_ptr_LagValAtCollocPtsInCluster();
		pvec2dV2 = pclsClusterLinear2->get_ptr_LagValAtCollocPtsInCluster();
		transpose(pvec2dV2, &vec2dAns);
		multiply_matrix(pvec2dX, &vec2dAns, &vec2dTemp);
		multiply_matrix(pvec2dV1, &vec2dTemp, &vec2dAns);
		make_matrix_for_near(&vec2dAns, pclsClusterLinear1, pclsClusterLinear2);
	}

/*****************************************************************
make_matrix_for_near - It is used for making Matrix for Cluster Pair Near. It is used for error calculation.
Arguments:
pvec2dX			-	Pointer to Vector of Vector Double. This is a temporary variable denoting the Function value at chebyshev nodes in the Cluster Pair Near.
pclsClusterLinear1	-	Pointer to Class ClusterLinear. This is a temporary variable denoting the ClusterLinear in 1st dimension forming the Cluser Pair Near.
pclsClusterLinear2	-	Pointer to Class ClusterLinear. This is a temporary variable denoting the ClusterLinear in 2nd dimension forming the Cluser Pair Near.
*****************************************************************/	
	void clsClusteringApproximation::make_matrix_for_near(vector< vector <double> >* pvec2dX, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2){
		int iLoopVari, iLoopVarj;
        int iStart1, iStart2;
        iStart1 = pclsClusterLinear1->get_iStart();
		iStart2 = pclsClusterLinear2->get_iStart();
		
		for (iLoopVari = pclsClusterLinear1->get_iStart(); iLoopVari <= pclsClusterLinear1->get_iEnd(); iLoopVari++){
			for (iLoopVarj = pclsClusterLinear2->get_iStart(); iLoopVarj <= pclsClusterLinear2->get_iEnd(); iLoopVarj++){
				vec2dMatrix[iLoopVari][iLoopVarj] += (*pvec2dX)[iLoopVari - iStart1][iLoopVarj - iStart2];
			}
		}		
	}

/*****************************************************************
calculate_error_inf - It is used for calculating the error in the infinity norm. It is used for error calculation.
*****************************************************************/	
	double clsClusteringApproximation::calculate_error_inf(){
            return(calculate_error_inf(&vecdFinalAnswerDirect, &vecdFinalAnswer));
    }

/*****************************************************************
calculate_error_inf - It is used for calculating the error in the infinity norm. It is used for error calculation.
Arguments:
pvecdTemp1	-	Pointer to Vector Double. This is one of the vectors for finidng the infinity norm error.
pvecdTemp2	-	Pointer to Vector Double. This is one of the vectors for finidng the infinity norm error. 
*****************************************************************/	
	double clsClusteringApproximation::calculate_error_inf(vector<double>* pvecdTemp1, vector<double>* pvecdTemp2){
		int iLoopVari;
		double dMaxError;
		double dErrorRelative;
		dMaxError = 0;
		
		for (iLoopVari = 0; iLoopVari < pvecdTemp1->size(); iLoopVari++){
			if (fabs((*pvecdTemp1)[iLoopVari]) < EPS2)
				dErrorRelative = fabs((*pvecdTemp1)[iLoopVari] - (*pvecdTemp2)[iLoopVari])/EPS2;
			else
				dErrorRelative = fabs((*pvecdTemp1)[iLoopVari] - (*pvecdTemp2)[iLoopVari])/fabs((*pvecdTemp1)[iLoopVari]);
			if(dErrorRelative > dMaxError)
				dMaxError = dErrorRelative; 
		}
		return (dMaxError*100);
	}	
	
/*****************************************************************
generate - For generating data in case of single entry.
*****************************************************************/
	void clsClusteringApproximation::generate(){
		double dStart1, dStart2, dEnd1, dEnd2;
		
		put_NumberofPts();
		put_Degree();
		put_Limits(0, &dStart1, &dEnd1);
		put_Limits(0, &dStart2, &dEnd2);
        put_AdmissibilityCoefficient();
		printf("\n\n");
        generate_CollocationPts(dStart1, dEnd1, veciNumberofPts[0], &vec2dCollocationPts[0], get_iOperation(0));
		generate_CollocationPts(dStart2, dEnd2, veciNumberofPts[1], &vec2dCollocationPts[1], get_iOperation(1));
		generate_Val();
		initialize_vector(&vecdFinalAnswer, veciNumberofPts[0], 0);
		initialize_vector(&vecdFinalAnswerDirect, veciNumberofPts[0], 0);
		initialize_matrix(&vec2dMatrix, veciNumberofPts[0], veciNumberofPts[1], 0);
		printf("Generating done\n");
	}

/*****************************************************************
generate_Val - For generating Vector C.
*****************************************************************/	
	void clsClusteringApproximation::generate_Val(){
		int iLoopVari;
		
		printf("Generating Values.......");
        vecdVal.resize(veciNumberofPts[1]);
		srand(time(0));
		for (iLoopVari = 0; iLoopVari < veciNumberofPts[1]; iLoopVari++){
			vecdVal[iLoopVari] = rand();
		}
		printf("done\n");
	}

/*****************************************************************
generate_different - For generating data in case of multiple Number of Points, Degree, Admissibility Coefficient.
				Upto iNumberofTimes computation for Time Complexity Computation  and processing and writing all the results
*****************************************************************/		
	void clsClusteringApproximation::generate_different(){
		double dStart1, dStart2, dEnd1, dEnd2;
		char* pchName;
		char* pchPts;
		char* pchAdmissiblity;
		char* pchExt = ".txt";
		int iNumberofTimes;
		int iNumberofPtsStart, iNumberofPtsEnd;
		int iNumberofPtsIncrements;
		vector<int> veciNumberofPtsStep;
		int iDegreeStart, iDegreeEnd, iDegreeStep;
		double dAdmissibleStart, dAdmissibleEnd, dAdmissibleStep;
		int iNumberofPts;
		int iCounterPts, iCounterDegree, iCounterAdmissibility;
		int* piCounterPtsIncrement;
		int iLoopVari;
		clock_t clkStart, clkStop, clkDiff;
		time_t timeStart, timeStop, timeDiff;
		double dErrorInf, dErrorFrobenius;
		ofstream fileCount;
		ofstream fileError;
		ofstream fileTime;
		
		iNumberofTimes = get_iNumberofTimes();
		pchName = new char[200];		
		pchPts = new char[2];
		pchAdmissiblity = new char[2];
		
		veciNumberofPtsStep.resize(1);
		printf("Enter the range for Number of Points and Increment ? ");
		scanf("%d %d %d", &iNumberofPtsStart, &iNumberofPtsEnd, &veciNumberofPtsStep[0]);
		piCounterPtsIncrement = new int;
		*piCounterPtsIncrement = 1;
		if (veciNumberofPtsStep[0] == 0){
                                   piCounterPtsIncrement = &iCounterPts;                                   
		   printf("Enter the Number of Increments ? ");
		   scanf("%d", &iNumberofPtsIncrements);
		   veciNumberofPtsStep.resize(iNumberofPtsIncrements);
		   for (iLoopVari = 0; iLoopVari < iNumberofPtsIncrements; iLoopVari++){
               printf("Enter Increments (%d) ? ", iLoopVari + 1);
               scanf("%d", &veciNumberofPtsStep[iLoopVari]);    
           }
        }    
        put_Limits(0, &dStart1, &dEnd1);
		put_Limits(1, &dStart2, &dEnd2);
		printf("Enter the range for Degree and Increment ? ");
		scanf("%d %d %d", &iDegreeStart, &iDegreeEnd, &iDegreeStep);
		printf("Enter the range for Admissibility Coefficient and Increment ? ");
		scanf("%lf %lf %lf", &dAdmissibleStart, &dAdmissibleEnd, &dAdmissibleStep);
		strcpy(pchName, pchDirPath);
		strncat(pchName, "Count.txt", 9);
		fileCount.open(pchName);
		strcpy(pchName, pchDirPath);
		strncat(pchName, "Error.txt", 9);
		fileError.open(pchName);
		strcpy(pchName, pchDirPath);
		strncat(pchName, "Time.txt", 8);
		fileTime.open(pchName);		
		iCounterPts = 0;
		for(iNumberofPts = iNumberofPtsStart; iNumberofPts <=iNumberofPtsEnd; iNumberofPts += veciNumberofPtsStep[*piCounterPtsIncrement - 1]){
			veciNumberofPts[0] = iNumberofPts;
			veciNumberofPts[1] = iNumberofPts;
			generate_CollocationPts(dStart1, dEnd1, veciNumberofPts[0], &vec2dCollocationPts[0], get_iOperation(0));
			generate_CollocationPts(dStart2, dEnd2, veciNumberofPts[1], &vec2dCollocationPts[1], get_iOperation(1));	
			generate_Val();
			iCounterDegree = 0;
			for (iDegree = iDegreeStart; iDegree <= iDegreeEnd; iDegree += iDegreeStep){
				iCounterAdmissibility = 0;
				dAdmissibilityCoefficient = dAdmissibleStart;
				while (dAdmissibilityCoefficient <= dAdmissibleEnd){
					timeDiff = 0;
					clkDiff = 0;
					clk1Diff = 0;
					clk2Diff = 0;
					clk3Diff = 0;
					clk4Diff = 0;
					for (iLoopVari = 0; iLoopVari < iNumberofTimes; iLoopVari++){
						initialize_vector(&vecdFinalAnswer, veciNumberofPts[0], 0);
						initialize_vector(&vecdFinalAnswerDirect, veciNumberofPts[0], 0);
						initialize_matrix(&vec2dMatrix, veciNumberofPts[0], veciNumberofPts[1], 0);
						initialize_matrix(&vec2dMatrixTemp, veciNumberofPts[0], veciNumberofPts[1], 0);
						initialize_vector(&vecdVectorTemp, veciNumberofPts[0], 0);						
						timeStart = time(NULL);
						clkStart = clock();
						process();
						timeStop = time(NULL);
						clkStop = clock();					
						timeDiff += timeStop - timeStart;
						clkDiff += clkStop - clkStart;
						if (iLoopVari != iNumberofTimes - 1)
						clear_pointers();
					}
					fileTime << iNumberofPts << "," << iDegree << "," << dAdmissibilityCoefficient << "," << timeDiff/iNumberofTimes << "," << clkDiff/iNumberofTimes << ",";
					fileTime << clk1Diff/iNumberofTimes  << "," << clk2Diff/iNumberofTimes << "," << clk3Diff/iNumberofTimes << "," << clk4Diff/iNumberofTimes << "\n";	
					make_matrix();
					count();
					fileCount << iNumberofPts << "," << iDegree << "," << dAdmissibilityCoefficient << "," << iMaxNear << "," << iMaxFar << ",";
					fileCount << iNumberofClusterLinearInNear1 << "," << iNumberofClusterLinearInFar1 << "," << iNumberofClusterLinearInNear2 << "," << iNumberofClusterLinearInFar2 << "\n";
					if (iDegree == iDegreeStart){
						*pchPts = char(65 + iCounterPts);
						*pchAdmissiblity = char(65 + iCounterAdmissibility);
						strcpy(pchName, pchDirPath);
						strncat(pchName,"Information", 11);
						strncat(pchName, pchPts, 1);
						strncat(pchName, pchAdmissiblity, 1);
						strcat(pchName, pchExt);				
						write_information(pchName);
						strcpy(pchName, pchDirPath);
						strcpy(pchName, pchDirPath);
						strncat(pchName,"ClusterLinear0_", 15);
						strncat(pchName, pchPts, 1);
						strncat(pchName, pchAdmissiblity, 1);
						strcat(pchName, pchExt);				
						write_clusterLinear(pchName, 0);
						strcpy(pchName, pchDirPath);
						strncat(pchName,"ClusterLinear1_", 15);
						strncat(pchName, pchPts, 1);
						strncat(pchName, pchAdmissiblity, 1);
						strcat(pchName, pchExt);				
						write_clusterLinear(pchName, 1);
						strcpy(pchName, pchDirPath);
						strncat(pchName,"ClusterPairNear_", 32);
						strncat(pchName, pchPts, 1);
						strncat(pchName, pchAdmissiblity, 1);
						strcat(pchName, pchExt);				
						write_clusterLinearInClusterPairNear(pchName);
						strcpy(pchName, pchDirPath);
						strncat(pchName,"ClusterPairFar_", 31);
						strncat(pchName, pchPts, 1);
						strncat(pchName, pchAdmissiblity, 1);
						strcat(pchName, pchExt);				
						write_clusterLinearInClusterPairFar(pchName);										
					}
					dErrorFrobenius = calculate_answer_direct(&vecdFinalAnswerDirect);
					dErrorInf = calculate_error_inf(&vecdFinalAnswerDirect, &vecdFinalAnswer);
					fileError << iNumberofPts << "," << iDegree << "," << dAdmissibilityCoefficient << "," << dErrorFrobenius << "," << dErrorInf << "\n";
					multiply (&vec2dMatrixTemp, &vecdVal, &vecdVectorTemp);
					printf("The error by Frobenius Norm is %lf & and by D-inf Norm is %lf\n", dErrorFrobenius, dErrorInf);
					dAdmissibilityCoefficient += dAdmissibleStep;
					iCounterAdmissibility++;
					printf("\n");
					clear_pointers();
				}
				iCounterDegree++;				
			}
			iCounterPts++;
		}
		write_informationMain(iCounterPts, iCounterDegree, iCounterAdmissibility);
		fileCount.close();
		//getch();
		delete pchName;
		delete pchPts;
		delete piCounterPtsIncrement;
	}
