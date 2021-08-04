This is the code in MATLAB associated with the manuscript: “Feasibility conditions of ecological models: Unfolding the links between model parameters”  by Mohammad AlAdwani, Serguei Saavedra.

We provided four examples: one for two-species Lotka-Volterra (LV) model with simple higher-order terms, a second for two-species  LV model with type III functional responses, a third for LV model with simple higher-order terms, and a fourth for LV model with higher-order interactions. To run the code you additionally need to (1) download PHCpack libraries from http://www.math.uic.edu/~jan/PHClab1.0.4.tar.gz (2) add add into the folder the exec file: http://homepages.math.uic.edu/~jan/download.html (3) run the code inside the folder where you saved PHCpack libraries

The code provides the number of free-equilibrium points and the probability of feasibility. The code is fully commented for its use and potential extensions.

The codes run MATLAB symbolic toolbox and each supplements the associated detailed example in the manuscript’s appendix. For extension and changing parameter values (fixed parameters) and ranges (varied parameters):

In Example 1, no parameter restriction was imposed.
In Example 2, change the parameter values and ranges in lines 732-742 and change the plotted quantities in lines 780-799.
In Example 3, change the parameter values and  in lines 442-461 and change the plotted quantities in lines 500-520.
In Example 4, change the parameter values in lines 6-18 &  parameter ranges (variables) in 748-753 then change the plotted quantities in lines 828-847.
