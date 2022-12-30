# Overview

The goal of this project was to simulate the dynamics of a 
cube boucing around inside a larger cube. This can be though of
as similar to a single die sitting on a table, and a cup being 
placed over it and shaken. A more in depth description can be found
in the PDF writup, `writup.pdf`

This problem, especially with all the possible impacts would be 
very difficult and tedious to solve using Newton's laws directly, 
however using the Euler-Lagrange equations, the problem is much
more doable. This is however somewhat computationally intensive.

All code for the simulation can be found in the `src` directory

- `main.ipynb` is the primary code that does everything. 

- `helpers.py` defines some useful helper functions.

- The `latex-report` directory is the LaTeX source code used to 
create the writup.

# Demo 

![final_demo](https://user-images.githubusercontent.com/45540813/210030830-d86044cf-20b7-4337-b2ec-2c696aa04525.gif)
