#include "cluster.hpp"

int main(){
    int iAnswer;
	clsClusteringApproximation clsClusteringApproximation1(2);
	printf("Enter your Choice...\n");
	//printf("\t1)Generate Collocation Points in specified Interval uniformly\n");
	printf("\t1)Enter Data from keyboard\n");
	printf("\t2)Generate Collocation Points in specified Interval uniformly, with varying degree and Number of Points\n");
	printf("Your Choice ? ");
	scanf("%d", &iAnswer);
    clsClusteringApproximation1.set_DirPath("/scratch");
	switch(iAnswer){
		//case 1:
		//	clsClusteringApproximation1.generate();
		//	break;
		case 1:
			clsClusteringApproximation1.initialize_process();
			break;
		case 2:
			clsClusteringApproximation1.generate_different();
			return(0);
			break;
	}
	
	clsClusteringApproximation1.process();
	clsClusteringApproximation1.make_matrix();
	clsClusteringApproximation1.write_informationMain(1, 1, 1);
	clsClusteringApproximation1.write_information("InformationAA.txt");
	clsClusteringApproximation1.write_clusterLinear("ClusterLinear0_AA.txt", 0);
	clsClusteringApproximation1.write_clusterLinear("ClusterLinear1_AA.txt", 1);
    clsClusteringApproximation1.write_clusterPairFar("ClusterPair_FarAA.txt");
	clsClusteringApproximation1.write_clusterPairNear("ClusterPair_NearAA.txt");
	
	
	//dErrorFrobenius = clsClusteringApproximation1.calculate_answer_direct();
	//dErrorInf = clsClusteringApproximation1.calculate_error_inf();
	//printf("The error by Frobenius Norm is %lf & and by L-inf Norm is %lf\n", dErrorFrobenius, dErrorInf);
	
	//////////////////////clsClusteringApproximation1.generate_different();
    //getch();    		
	clsClusteringApproximation1.clear_pointers();
	//system("matlab -r test");
	//getch();
	
	return(0);

}
