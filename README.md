# Abaqus UEL Elasticity
 
 
This repository contains the Fortran source code for small strain isotropic linear elastic user element (UEL) subroutine and example input files for Abaqus/Standard. Standard displacement-based finite element formulation has been adopted in this subroutine. The purpose of the project is to make users familiar with developing the UEL subroutine in Abaqus/Standard using a standard textbook formulation. This source code contains necessary subroutines related to element operations and matrix algebra.


> [!WARNING]
> This repository is not meant to be a complete guideline or tutorial for Abaqus user element (UEL) subroutines. This Abaqus feature is intended for advanced users and requires an understanding of finite element formulation and substantial programming. Users are recommended to consult textbooks on finite element analysis, Fortran programming, and Abaqus user documentation before developing their UEL subroutines.



## Obtaining the file

If you have `git` installed, you can clone the repository to your local machine using
```bash
git clone https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity.git
```
Alternatively, you can download the files in a `zip` folder in this repository using the `code` drop-down menu on the top right corner

To receive updates in this repository, you can also `fork` the repository and sync as updates are deployed, develop your code by creating a separate branch, and propose updates using the `pull` and `merge` features of GitHub.



## Description of the repository

| File name     | Description   |
| :---------    | :-----------  |
| `uel_mech.for` | is the Fortran source code that implements the isotropic linear elastic user element. The main `UEL` subroutine performs all the initial checks but the main calculations are performed in a subsequent subroutine. The source code includes additional subroutines with Lagrangian interpolation functions for 4 types of 2D continuum elements (Tri3, Tri6, Quad4, and Quad8) and 4 types of 3D continuum elements (Tet4, Tet10, Hex8, Hex20) and Gaussian quadratures with reduced and full integration schemes. Body force and traction boundary conditions were not been included in this implementation, however, these can be applied by overlaying standard Abaqus elements on the user elements (to be discussed in the **Modeling in Abaqus** section). Since Abaqus/Viewer does not provide native support for visualizing user elements, an additional layer of elements with the same element connectivity has been created and results at the integration points of the elements are stored using the `UVARM` subroutine. |
| `addElemMech.py` | is a Python code that modifies a simple input file and adds the overlaying dummy elements on the user elements. For complicated input files, this will not work properly and modification of this code will be required (optional). |
| `<>.inp` | are the example input files prepared to be executed with the user element subroutine. Since the user-defined elements share the same topology as one of the Abaqus built-in elements, those models were built in Abaqus/CAE and then exported as input files. Later those input files were modified to include keywords and data to include user element definitions, properties, and overlaying dummy elements. |
| `runAbq.ps1` | is a PowerShell batch file that can execute the user subroutine and specified input file from the PowerShell terminal (optional).
| Theory.pdf | is a brief summary of the theory and algorithm used to implement the UEL source code. |
| Abaqus Docs.pdf | A collection of publicly available Abaqus documentation in PDF format which are related to Abaqus UEL. The web versions of these documents are available at https://help.3ds.com. |



## Modeling in Abaqus

Since the implemented user elements have the same topology as the built-in Abaqus elements, users can build a primary model in Abaqus/CAE and then export the input (`.inp`) file. Once the input file is available, as described in the Abaqus documentation, the following information needs to be modified.

### Element definition

The first modification is to introduce the definition of the user element being used in the analysis with the element name, coordinates, nodes, and the number of real and integer properties required for that element.


### Properties

For isotropic linear elastic elements developed in this UEL, users need to specify the following properties:
- Young's modulus, $E$
- Poisson's ratio, $\nu$
- Number of integration points, `nInt`
- Number of post-processed variables, `nPostVars`


### Dummy elements for visualization

To visualize the results, an additional set of built-in Abaqus elements with the same element connectivity as the user element has been created in the input file. These additional elements (so-called dummy elements) have negligible elastic properties and thus will not affect the results. If you are using a reduced integration element from the user subroutine, then use the same type of element from Abaqus as dummy elements.



## Configuring Abaqus and executing the subroutines

To run user subroutines in Abaqus, you will need to install Intel Visual Studio and Intel oneAPI package and link them with Abaqus. Follow [this tutorial](https://www.bibekanandadatta.com/blog/2021/link-intel-and-vs-abaqus-2020/) if you have not done it before.

Make sure that the source code and input file are in the same directory. Using the `Abaqus command line terminal` or `cmd terminal` or `PowerShell terminal`, you can execute the following command from the directory to execute the subroutine.

```bash
abaqus interactive double analysis job=<your_job_name> input=<input_file_name.inp> user=<source_code.for>
```
Specify the variable names (inside < >) in the above command as needed. For additional information on executing user subroutines, check the Abaqus user manual.

If you use the PowerShell-based terminal, you can also execute the subroutine by running the `runAbq.ps1` file. Make sure to check the input file name in the file.
```bash
./runAbq
```


## Documentation

Users should consult Abaqus documentation, a standard finite element textbook for solid mechanics, and textbooks and tutorials related to Fortran programming. If further documentation is required, it can be made available upon request.


## Citation

In case you use this subroutine for educational or research purposes, please cite this source.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11074933.svg)](https://doi.org/10.5281/zenodo.11074933)
