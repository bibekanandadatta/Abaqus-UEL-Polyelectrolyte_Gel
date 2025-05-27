# Abaqus-UEL-Polyelectrolyte_Gel

Abaqus/Standard user element subroutines for coupled electro-chemo-mechanics of polyelectrolyte hydrogel. Publication for the work is under-preparation, and the repository is not currently available to public.


## Requirements


## Brief description of the repository

### dev
This repository contains the same source code as the src directory except the additional source code 'test_uel_pegel.for' for numerically testing the implementation.
As of now, numerical testing is available for axisymmetric and plane strain elements. 

To compile the code using Intel Fortran compiler, use:

```
ifort \Qmkl -o pegel test_uel_pegel.for
```

To run the code, use: 
```
.\pegel
```


### docs
This directory contains the pre-print of the article and the supplementary material. Supplementary material provides brief description of the different modules and a simple guideline on how to use the code in Abaqus/Standard. The directory contains a PDF version of Abaqus/Standard documentation for user element. For full documentation, visit: https://help.3ds.com. You will need to open an account to access the online documentation.

### src
This directory contains the source code of the implementation.



### test
This directory contains the Abaqus input files (.inp) of a few different examples as presented in the published article. Although we performed parametric studies for the paper, we did not include those input files as they are created from the baseline cases provided here.


## Source codes



## Abaqus test files