#include "abelintegralequation.h"

int main() {
    /* SAM_LISTING_BEGIN_1 */
    {
        auto u = [](double t) { return 2./M_PI*sqrt(t); };
        auto y = [](double t) { return t; };
        
        double tau = 0.01;
        size_t N = round(1./tau);
        VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
        VectorXd u_ex(N+1);
        for(int i=0; i<N+1; ++i) {
            u_ex(i) = u(grid(i));
        }
    
        cout << "\nSpectral Galerkin\n" << endl;
        for(int p=2; p<=10; ++p) {
            VectorXd u_app = poly_spec_abel(y, p, tau);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "p = " << p << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
    }
    /* SAM_LISTING_END_1 */

    /* SAM_LISTING_BEGIN_4 */
    {
        auto u = [](double t) { return 2./M_PI*sqrt(t); };
        auto y = [](double t) { return t; };

        cout << "\n\nConvolution Quadrature, Implicit Euler\n"  << endl;
        for(int N=16; N<=2048; N*=2) {

            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = u(grid(i));
            }
            
            VectorXd u_app = cq_ieul_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
        
        cout << "\n\nConvolution Quadrature, BDF-2\n"  << endl;
        for(int N=16; N<=2048; N*=2) {
            
            VectorXd grid = VectorXd::LinSpaced(N+1,0.,1.);
            VectorXd u_ex(N+1);
            for(int i=0; i<N+1; ++i) {
                u_ex(i) = 2./M_PI*sqrt(grid(i));
            }
            
            VectorXd u_app = cq_bdf2_abel(y, N);
            VectorXd diff  = u_ex - u_app;
            double err_max = diff.cwiseAbs().maxCoeff();
            cout <<   "N = " << N << setw(15)
                 << "Max = "
                 << scientific << setprecision(3)
                 << err_max << endl;
        }
    }
    /* SAM_LISTING_END_4 */
return 0;
}