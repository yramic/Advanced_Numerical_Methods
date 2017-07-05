#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

using namespace std;

/*****************************************************************
Class Cluster Linear - For creating instances of Cluster Linear.
Variables:
iStart					-	Integer. Indicates the Start of Linear Cluster as index of Collocation Point
iEnd						-	Integer. Indicates the End of Linear Cluster as index of Collocation Point
iLeft						-	Integer. Indicates the Left child of Linear Cluster as index number in Vector of Class Cluster Linear
iRight					-	Integer. Indicates the Right child of Linear Cluster as index number in Vector of Class Cluster Linear
vec2dLagValAtCollocPtsInCluster	-	Vector of Vector Double. Double dimension Vector for storing Lagrange Value at Collocation Points in Cluster. 
							Rows correspond to Collocation Points.
							Columns correspond to Degree.
veciAppearsIn				-	Vector of Integer. Single dimension Vector for storing the index of ClusterPairs in which this Cluster Linear appears.
vecdW					-	Vector of Double. Single dimension Vector for storing the value of W for the Cluster Linear of Dimension 1.
							W is defined by multiplication of ClusterPair`s Function Value Matrix and Lagrange Value at Collocation Points in Cluster.
							This is defined only in Cluster Linear of Dimesnion 1.
*****************************************************************/
class clsClusterLinear {
    int iStart;
    int iEnd;
    int iLeft;
    int iRight;
	vector < vector <double> > vec2dLagValAtCollocPtsInCluster;
    vector<int> veciAppearsIn;
	vector <double> vecdW;
    
	public:
	clsClusterLinear(int iTempStart, int iTempEnd, int iTempLeft, int iTempRight, int iDegree);
	~clsClusterLinear();
	void increase_size_LagValAtCollocPtsInCluster(int iDegree);
	int get_iStart();
	int get_iEnd();
	int get_iLeft();
	int	get_iRight();
	vector< vector <double> >* get_ptr_LagValAtCollocPtsInCluster();
    vector<int>* get_ptr_AppearsIn();
	vector<double>* get_ptr_vecdW();
	void print(int iIndexNumber);
};

/*****************************************************************	
Class Cluster Pair -  For creating instances of Cluster Pair.
Variables:
vecpclsClusterLinear			-	Vector of Pointer of Class Cluster Linear. Single Dimension Vector for storing pointers of type Class Cluster Linear. 
							It is used to indicate the Cluster Linear used in forming the Cluster Pair. 
vec2dFnValAtChebyNodeInCluster	-	Vector of Vector Double. Double Dimension Vector for storing the kernel function values at Chebyshnev Nodes in Cluster for Cluster Pair Far.
							Double Dimension Vector for storing the kernel function values at Collocation Points in Cluster Pair for Cluster Pair Near.
*****************************************************************/
class clsClusterPair {
    vector <clsClusterLinear*> vecpclsClusterLinear;
    vector< vector <double> > vec2dFnValAtChebyNodeInCluster;

	public:	
	clsClusterPair(int iDimension, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2, int iRows, int iCols);
	~clsClusterPair();
	vector< vector <double> >* get_ptr_FnValAtChebyNodeInCluster();
	clsClusterLinear* get_ptr_ClusterLinear(int iDimension);
	void print(int iIndexNumber, int iDegree, int iDimension);
};

/*****************************************************************	
Class Clustering Approximation - Class used for doing Clustering Approximation.
Variables:
vec2pclsClusterLinear		-	Vector of Vector of Pointer of Class Cluster Linear. It is used for storing the pointers to the Cluster Linear created.
						Rows indicate the Dimension Number.
						First Row contains the Linear Cluster formed in first Dimension.
						Second Row contains the Linear Cluster formed in Second Dimension and so on......
vecpclsClusterPairNear		-	Vector of Pointer of Class Cluster Pair. It is used for storing the pointers to the Cluster Pair Near created.
						It stores pointers to all the Cluster Pairs of Near type created.
vecpclsClusterPairFar		-	Vector of Pointer of Class Cluster Pair. It is used for storing the pointers to the Cluster Pair Far created.
						It stores pointers to all the Cluster Pairs of Far type created.
vec2dCollocationPts		-	Vector of Vector Double. Double Dimension Vector for storing the Collocation Points of different dimensions.
						Rows indicate the Dimension Number.
						First Row contains the Collocation Points in first Dimension.
						Second Row contains the Collocation Points in Second Dimension and so on......
vecdVal				-	Vector Double. It is used for storing the value of Vector C.
pvec2dFnValAtChebyTemp	-	Pointer to Vector of Vector Double. It is used in various functions to store the pointer of vec2dFnValAtChebyTemp ( a variable of Class Cluster Pair).
vec2dTempForTranspose	-	Vector to Vector of Double. It is a temporary double dimension vector used for storing Transpose of Matrices.
veciNumberofPts			-	Vector Integer. It is used for storing Number of Points in each Dimension.
vecdFinalAnswer			-	Vector Double. It is used for storing the Final Answer of the Computation.
vecdFinalAnswerDirect		-	Vector Double. It is used for storing the Final Answer computed by Direct Calculation.
vec2dMatrix			-	Vector of Vector Double. It is used for storing the matrix containing the values of Kernel Function at Collocation Points.
						It is required for Frobenius error calculation.
vec2dMatrixTemp		-	Vector of Vector Double. It is used for storing the matrix containing the approximated values of Kernel Function at Collocation Points computed by the clustering approxomation technique.
						It is required for Frobenius error calculation.
vecdVectorTemp			-	Vecto Double. It is a temporary Vector used for storing Vector.
iDegree				-	Integer. It is used for storing the degree of current Clustering Approximation Procedure.
iDimension				-	Integer. It is used for storing the highest dimension of current Clustering Approximation Procedure. In the procedure implemented it is always 2
iCounter				-	Integer. It is a temporary variable used for various times in the functions.
dAdmissibilityCoefficient	-	Double. It is the admissibility coefficient required in the admissibility criteria.
dClusterStart1			-	Double. It is the collocation point at the start of the cluster pair in dimension 1.
dClusterStart2			-	Double. It is the collocation point at the start of the cluster pair in dimension 2.
dStep1				-	Double. It is the step size for calculating the chebyshev nodes in dimension 1.
dStep2				-	Double. It is the step size for calculating the chebyshev nodes in dimension 2.
pchDirPath				-	Pointer to character. It is a pointer for stoing the directory path for storing files required by MATLAB to plot data. 
iMaxNear				-	Integer. It is the maximum number of occurences of cluster linear in Cluster Pair Far.
iMaxFar				-	Integer. It is the maximum number of occurences of cluster linear in Cluster Pair Near.
*****************************************************************/
class clsClusteringApproximation {
	vector< vector <clsClusterLinear*> > vec2pclsClusterLinear;
	vector<clsClusterPair*> vecpclsClusterPairNear;
	vector<clsClusterPair*> vecpclsClusterPairFar;
	vector< vector <double> > vec2dCollocationPts;
	vector<double> vecdVal;
	vector< vector <double> >* pvec2dFnValAtChebyTemp;
	vector< vector <double> > vec2dTempForTranspose;
	vector<int> veciNumberofPts;
	vector<double> vecdFinalAnswer;
	vector<double> vecdFinalAnswerDirect;
	vector< vector<double> > vec2dMatrix;
	vector < vector<double> > vec2dMatrixTemp;
    vector <double> vecdVectorTemp; 
	int iDegree;
	int iDimension;
	int iCounter;
	int iNumberofClusterLinearInFar1;
	int iNumberofClusterLinearInFar2;
	int iNumberofClusterLinearInNear1;
	int iNumberofClusterLinearInNear2;

	double dAdmissibilityCoefficient;
	double dClusterStart1, dClusterStart2;
	double dStep1, dStep2;
	char* pchDirPath;
	int iMaxNear, iMaxFar;
	clock_t clk1Start, clk1Stop, clk1Diff;
	clock_t clk2Start, clk2Stop, clk2Diff;
	clock_t clk3Start, clk3Stop, clk3Diff;
	clock_t clk4Start, clk4Stop, clk4Diff;

	public:
	clsClusteringApproximation(int iDimension);
	void set_DirPath(char* pchDirPath);
	void generate_Val();
	void generate_CollocationPts(double dStart, double dEnd, int iNumberofPts, vector<double>* pvecdCollocationPts, int iOperation);
	void write_informationMain(int iTotalDegree, int iTotalPoints, int iTotalAdmissible);
	void write_information(char* pcharFileName);
	void write_clusterPairNear(char* pcharFileName);
	void write_clusterPairFar(char* pcharFileName);
	void write_clusterLinear(char* pcharFileName, int iDimension);
	void write_clusterLinearInClusterPairFar(char* pcharFileName);
	void write_clusterLinearInClusterPairNear(char* pcharFileName);
	void put_NumberofPts();
	void put_Limits(int iIndex, double* pdStart, double* pdEnd);
	void put_Degree();
	void put_CollocationPts();
	void put_Val();
	void put_AdmissibilityCoefficient();
	void multiply (vector< vector<double> >* pvec2dTemp1, vector<double>* pvecdTemp2, vector <double>* pvecdAns);
	void multiply_matrix (vector< vector<double> >* pvec2dTemp1, vector< vector<double> >* pvec2dTemp2, vector < vector<double> >* pvec2dAns);
    void custom_multiply (vector< vector<double> >* pvec2dTemp1, vector<double>* pvecdTemp2, vector <double>* pvecdAns, clsClusterLinear* pclsClusterLinearTemp);
	void copy(vector<double>* pvecTemp, vector<double>* pvecAns);
	void add(vector<double>* pvecTemp1, vector<double>* pvecTemp2, vector<double>* pvecAns);
	void transpose(vector<  vector<double> >* pvecTemp, vector< vector<double> >* pvecAns);
	void initialize_vector(vector<double>* pvecdTemp, int iSize, int iInitialValue);
	void initialize_matrix(vector< vector<double> >* pvec2dMatrix, int iRows, int iCols, int iInitialValue); 
	void add_ClusterLinear(int iDimension, int iStartTemp, int iEndTemp, int iLeftTemp, int iRightTemp);
	void build_ClusterTree_Linear(int iOffset, int iNumberofPtsInClusterLinear, int* piDimension);
	void add_ClusterPairNear(clsClusterLinear* pclsClusterLinear1Temp, clsClusterLinear* pclsClusterLinear2Temp);
	void add_ClusterPairFar(clsClusterLinear* pclsClusterLinear1Temp, clsClusterLinear* pclsClusterLinear2Temp);
	bool leaf(int iDimension, int iIndex);
	bool admissible_ClusterPair(int iIndex1, int iIndex2);
	double kernel_function(double dX, double dY);
	void lagrange_value_for_Cluster(clsClusterLinear* pclsClusterLinearTemp, double dClusterStart, double dStep, int iDimension);
	void fill_FnValAtChebyNode_for_Near(vector< vector<double> >* pvec2dFnValAtChebyTemp, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2);
	void fill_FnValAtChebyNode_for_Far(vector< vector<double> >* pvec2dFnValAtChebyTemp, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2);
	void check_and_fill_LagVal(clsClusterLinear* pclsClusterLinear, int iDimension);	
	void tree_traverse(clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2);
	void build_ClusterPair (int iIndex1, int iIndex2);
	void compute_S(bool* pbAtLeastOneFar, int iPairNumber, vector<double>* pvecdS);
	void compute_final_answer_for_near(int iIndex, int iPairNumber);
	void compute_final_answer_for_far(int iIndex, vector<double>* pvecdS);
	void finalAnswer();
	void make_matrix();	
	void print_ClusterLinear(int iDimension, int iIndexNumber);
	void print_ClusterLinear();
	void print_ClusterPairFar(int iIndexNumber);
	void print_ClusterPairNear(int iIndexNumber);
	void print_ClusterPairFar();
	void print_ClusterPairNear();
	void print(char* pchName, vector<int>* pveciPrint);
	void print(char* pchName, vector<double>* pvecdPrint);
	void print(char* pchName, vector< vector<double> >* pvec2dPrint);
	void clear_pointers();
	void initialize_process();
	void generate();
	void generate_different();
	void preprocess();
	void preprocess_compute();
    void process();
	void count(int* piMaxFar, int* piMaxNear, vector<clsClusterLinear*>* pvecpclsClusterLinear, int* piNumberofClusterLinearInFar, int* piNumberofClusterLinearInNear);
	void count();
	double calculate_answer_direct();
	double calculate_answer_direct(vector<double>* vecdFinalAnswerDirect);
	void make_matrix_for_far(vector< vector <double> >* pvec2dX, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2);
	void make_matrix_for_near(vector< vector <double> >* pvec2dX, clsClusterLinear* pclsClusterLinear1, clsClusterLinear* pclsClusterLinear2);
	double calculate_error_inf(vector<double>* pvecdTemp1, vector<double>* pvecdTemp2);
	double calculate_error_inf();
	int get_iOperation(int iDimension);
	int get_iNumberofTimes();
};

#endif // CLUSTER_HPP
