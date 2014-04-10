##Research Code for Fitting Multiple Gaussian Graphical Model


###Compilation

You need to have [R][R-website] installed before compiling the code.
You have to install R from source on a linux machine.

To compile the C/C++ code, just type `make` in the root directory. To
recompile, first type `make clean`, then type `make`. 



###Example
To run the 
provided example, just excute `example.R` file in R  

~~~ R
source("src/R/example.R")
~~~ 

The solution path for convex and non-convex method are stored in 
a list object `sol_path`, consisting of two arrays: `sol_nonconvex`
 and `sol_convex`, with the first one refering to the solution path
 using truncated $$L_1$$ penalty and the second one refering to that of 
 $$L_1$$ penalty. The array indices are arranged in a way such that 
 `sol_nonconvex[,,j,i]` is refering to the solution corresponding 
 to grid point `(Lambda1.vec[i], Lambda2.vec[j])`.   

[R-website]: http://www.r-project.org/
