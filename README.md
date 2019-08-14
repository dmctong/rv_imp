# Code from Tong and Hernandez, 2019

The code contained in this github is intended for transparency in the scientific process. The code itself is written directly for use on the UCSF computational cluster running SGE, and is therefore not suited to general usage.

Here is a brief rundown of the code:

* pipeline1 simulates the genetic data using the Tennessen demographic model (2012) and msprime (Kelleher et al).
* pipeline1a does some file type manipulations and pre-calculations
* pipeline2 simulates the Uricchio phenotype model (Uricchio et al, 2016)
* pipeline3 does the imputation using IMPUTE4 (Marchini) and runs RVATs using rvtests
* pipeline4 compiles the results
