# Examples for 2D-Parametric BEM

* These examples illustrate the use of 2D-Parametric BEM library as well as validating it by using known problems. 

## Single Layer Galerkin Matrix example

* In this example, the 2D-Parametric BEM code is used to find the Galerkin Matrix for the bi-linear form induced by Single Layer BIO, in the case of a unit circle in 2-D centered at the origin.
* For a circle with uniform mesh, the obtained Galerkin Matrix is expected to be a Circulant Matrix therefore "Fourier vectors" are its eigenfunctions.
* The Fourier harmonics are eigenfunctions of the Single Layer potential operator on the circle for eigenvalues k^{-1}. Thus, the sorted eigenvalues for the obtained Galerkin Matrix are expected to behave approximately like k^{-1}
* The executable can be build using the command "make single_layer_test" from the build directory. Using the python script "plot_sl.py", the inverse of sorted eigenvalues can be plotted and the expected behaviour can be verified.

## Double integration of Logarithmic Kernel

* In BEM, double integrals of a singular kernel are encountered everywhere. This example deals with the double integral of the Logarithmic kernel "log|x-y|" and demonstrates the need for transformations tackling both the variables instead of just one for numerical integration. In particular, it demonstrates the problem when the double integration is performed by integrating the inner integral first.
* The example has two parts: 
	* Plotting the integral of the Logarithmic kernel "log |x-y| " with respect to x in the range [-1,1], as a function of y. This function of y is what's obtained from the inner integral in the double integral. 
	* Finding the type of convergence by integrating the analytic expression obtained for the inner integral, using Gauss Quadrature for two sets of y values: [-1,1] & [1,2].
The first case has a line of singularity which runs along the diagonal of the parameter domain, whereas the second case contains a point singularity at one of the corners of the parameter domain.
* To build the executable, use the command "make log" from the build directory. The python script "plot.py" can be used to obtain the relevant plots. 


