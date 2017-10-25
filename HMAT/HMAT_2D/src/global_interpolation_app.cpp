#include "../include/global_interpolation_app.hpp"
#include "../include/kernel.hpp"
#include "../include/cheby.hpp"
#include <iostream>
#include <Eigen/Dense>

/*GlobalInterpolationApp::GlobalInterpolationApp(GlobalSmoothKernel gskernel, std::vector<Point> pp, int n):
    GSK_(gskernel), PPoints_(pp)
{}

GlobalInterpolationApp::GlobalInterpolationApp(ConstantKernel ckernel, std::vector<Point> pp, int n):
    CK_(ckernel), PPoints_(pp)
{}

GlobalInterpolationApp::GlobalInterpolationApp(PolynomialKernel ckernel, std::vector<Point> pp, int n):
    PK_(ckernel), PPoints_(pp)
{}
GlobalInterpolationApp::GlobalInterpolationApp(GaussKernel ckernel, std::vector<Point> pp, int n):
    GK_(ckernel), PPoints_(pp)
{}*/
GlobalInterpolationApp::GlobalInterpolationApp(Kernel* ckernel, std::vector<Point> pp, int n):
    K_(ckernel), PPoints_(pp)
{}
// approximate matrix-vector multiplication
Eigen::VectorXd GlobalInterpolationApp::mvProd(Eigen::VectorXd& c, unsigned deg)
{
    double maxX,minX,maxY,minY;
    maxX = PPoints_.begin()->getX();
    minX = PPoints_.begin()->getX();
    maxY = PPoints_.begin()->getY();
    minY = PPoints_.begin()->getY();
    for (std::vector<Point>::iterator it=PPoints_.begin()+1; it!=PPoints_.end(); it++){
        if (it->getX()>maxX) {
            maxX = it->getX();
        }
        else if (it->getX()<minX) {
            minX = it->getX();
        }
        if (it->getY()>maxY) {
            maxY = it->getY();
        }
        else if (it->getY()<minY) {
            minY = it->getY();
        }
    }
    Cheby cbx(minX, maxX, deg);
    Cheby cby(minY, maxY, deg);
    Eigen::VectorXd tkx = cbx.getNodes(); // Chebyshew nodes for x axis
    //std:: cout << "tkx" << tkx << std::endl;
    Eigen::VectorXd wkx = cbx.getWghts(); // weights of Lagrange polynomial for x axis
    //std:: cout << "wkx" << wkx << std::endl;
    Eigen::VectorXd tky = cby.getNodes(); // Chebyshew nodes for y axis
    //std:: cout << "tky" << tky << std::endl;
    Eigen::VectorXd wky = cby.getWghts(); // weights of Lagrange polynomial for y axis
    //std:: cout << "wky" << wky << std::endl;

    int ppts = PPoints_.size();   // farfield like computation
    Eigen::MatrixXd V = Eigen::MatrixXd::Constant(ppts, (deg+1)*(deg+1), 1);
    for(unsigned i=0; i<=ppts-1; ++i) { // calculation of Vx combined with Vy
        for(unsigned j1=0; j1<=deg; ++j1) {
            for(unsigned k1=0; k1<j1; ++k1) {
                for(unsigned j2=0; j2<=deg; ++j2) {
                    V(i,j1*(deg+1) + j2) *= (PPoints_[i].getX() - tkx[k1]);


                }
            }
            // Skip "k1 == j1"
            for(unsigned k1=j1+1; k1<=deg; ++k1) {
                for(unsigned j2=0; j2<=deg; ++j2) {
                    V(i,j1*(deg+1) + j2) *= (PPoints_[i].getX() - tkx[k1]);

                }
            }
            for(unsigned j2=0; j2<=deg; ++j2) {
                for(unsigned k2=0; k2<j2; ++k2) {
                    V(i,j1*(deg+1) + j2) *= (PPoints_[i].getY() - tky[k2]);
                }
                // Skip "k2 == j2"
                for(unsigned k2=j2+1; k2<=deg; ++k2) {
                    V(i,j1*(deg+1) + j2) *= (PPoints_[i].getY() - tky[k2]);
                }
                V(i,j1*(deg+1) + j2) *= wkx(j1)*wky(j2);
            }

        }
    }
    Eigen::VectorXd Vc((deg+1)*(deg+1));
    Eigen::MatrixXd X((deg+1)*(deg+1),(deg+1)*(deg+1));
    Vc = V.transpose() * c;
    for(int i=0; i<=deg; i++){
        for(int j=0; j<=deg; j++){
            for(int k=0; k<=deg; k++){
                for(int l=0; l<=deg; l++){
                    X(i*(deg+1)+j,k*(deg+1)+l) = (*K_)(cbx.getNodes()[i],cby.getNodes()[j],cbx.getNodes()[k],cby.getNodes()[l]);
                }
            }
        }
    }

    Eigen::VectorXd f;
    f = V * X * Vc;
    return f;
}

