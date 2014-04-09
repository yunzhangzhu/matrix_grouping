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
a list object `sol_path`.

[R-website]: http://www.r-project.org/
