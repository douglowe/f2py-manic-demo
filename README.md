# f2py-manic-demo
Testing the construction of Python modules from complex Fortran model code.

Directory Structure:

library_files:  Contains compilation and library directories for MANIC
|-> MANIC_SourceCode: Compilation directory for MANIC fortran model
|-> MANIC_Python_Interface: Compilation directory for MANIC-Python interface

fortran_model_run_dir: Contains input file needed for running MANIC testcase
|-> example_output: Contains example output (and log file) from successful model run
