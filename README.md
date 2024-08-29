# antibiotic-combination-selection-dynamics
Data and code for the paper "Private benefit of β-lactamase dictates selection dynamics of combination antibiotic treatment"

Data and code are organized per figure. Matlab code (.m) files in the outermost folder are generally used to generate final figures. Matlab codes in inner folders are sometimes used to generate intermediate data (.mat, .csv, or .xlsx). Jupyter Notebook (.ipynb) files may also be used to generate intermediate data. Data is usually stored in .xlsx or .csv format. Microscopy images are not included for size. 

For Figure 5CD and S7-9, data in "isolates" and ODE systems in "simulation" are used to run fitting in "scipy_optimization". Fitting results are stored in "estimates". Fit is checked and simulations using estimated parameters are saved in "para_evaluation". This data is converted to .csv or .mat and used to run GRprocessing.m and Figure_5CD in the outermost folder. 
