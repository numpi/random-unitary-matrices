# Sampling the eigenvalues of random unitary and orthogonal matrices

This repository contains the software for the manuscript [1], which describes a method for sampling the eigenvalues of Haar-distributed unitary and orthogonal matrices from the groups U(n), O(n), SU(n), and SO(n). The software is for MATLAB, and contains scripts replicating the experiments included in the paper.

## Using the software
The software relies on the unitary QR solver found in EISCOR [2], a Fortran 90 library for core-chasing eigenvalue algorithms. A script called ```compile_eiscor.m``` is included that automatically downloads and compile EISCOR and the MEX interface.

To get started, you need to download the code in this repository and run the command
```compile_eiscor```; the output should look similar to the following.
```
>> compile_eiscor
If necessary, I will download and install eiscor. 
You need to have the following packages installed: git, gfortran, make. 

Should I proceed? [yn]: y
Building with 'gfortran'.
MEX completed successfully.
```

Then, you may use the function ```e = sample_eigs(n, true)``` to sample the eigenvalues of matrices in U(n), or ```e = sample_eigs(n, false)``` to sample eigenvalues in O(n). In addition, it is possible to sample the conditional probability obtained by prescribing the
determinant equal to any complex number of modulus one. The most general syntax is as follows:
```
>> e = sample_eigs(n, complexoutput, det, useeiscor);
```
where:
 * ```n``` is the dimension of the unitary or orthogonal matrices to sample, and therefore the number of eigenvalues computed.
 * ```complexoutput``` specified if the complex or real groups are desired.
 * ```det``` is the assigned determinant. ```det = 1``` samples from the special groups SU(n) and SO(n); ```det = []``` samples from U(n) and O(n), and is the default value.
 * ```useeiscor``` forces the function to use EISCOR (if ```true```) or the standard QR solver (if ```false```); this is mostly useful to make complexity comparisons, or to skip the check for the EISCOR mex file that takes a non-trivial amount of time when dealing with small dimensions.
 
## Bugs

If you find bugs, please either open an issue on Github, or drop an email to [leonardo.robol@unipi.it](mailto:leonardo.robol@unipi.it) or [massimiliano.fasi@gmail.com](mailto:massimiliano.fasi@gmail.com).

## References
The algorithms implemented in the code are described in detail in [1]. 

[1] M. Fasi & L. Robol. Sampling the eigenvalues of random orthogonal and unitary matrices. In preparation.

[2] EISCOR - A Fortran 90 package for core-chasing eigenvalue algorithms. Github repository at [https://github.com/eiscor/eiscor](https://github.com/eiscor/eiscor).
