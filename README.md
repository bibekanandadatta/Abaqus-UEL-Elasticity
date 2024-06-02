# Abaqus UEL Elasticity
 
 
This repository contains the Fortran source code for small strain isotropic linear elastic user element (UEL) subroutine and example input files for Abaqus/Standard. Standard displacement-based finite element formulation has been adopted in this subroutine. The purpose of the project is to make users familiar with developing the UEL subroutine in Abaqus/Standard using a standard textbook formulation. This source code contains necessary subroutines related to element operations and matrix algebra.


> [!WARNING]
> Abaqus UEL subroutine feature is intended for advanced users and requires the understanding of finite element formulations and substantial programming. Users are recommended to consult textbooks on finite element analyses, Fortran programming, and Abaqus user documentation before developing their UEL subroutines. 



## Obtaining the files

If you have `git` installed, you can clone the repository to your local machine using
```bash
git clone https://github.com/bibekananda-datta/Abaqus-UEL-Elasticity.git
```

You can also `fork` the repository and sync as updates are deployed, develop your code by creating a separate branch, and propose updates using the `pull` and `merge` features of GitHub.

Alternatively, you can download the files in a `zip` folder in this repository using the `code` drop-down menu on the top right corner. In this approach, you will not receive any bug fixes and updates.



## Description of the repository

All the source codes are located in the `src` subdirectory and the Abaqus test cases are located in the `tests` subdirectory. The documentations are available in the `docs` subdirectory. Compiling the source code requires Intel oneMKL library. Compiling the source code requires the LAPACK library from the Intel oneMKL package. See below for the details.

|   File name  |  Description  |
| ----------   | ------------- |
| `uel_mech.for` | is the Fortran source code that implements the isotropic linear elastic user element. The main `UEL` subroutine performs all the initial checks but the main calculations are performed in a subsequent subroutine. The source code includes additional subroutines with Lagrangian interpolation functions for 4 types of 2D continuum elements (Tri3, Tri6, Quad4, and Quad8) and 4 types of 3D continuum elements (Tet4, Tet10, Hex8, Hex20) and Gaussian quadratures with reduced and full integration schemes. Body force and traction boundary conditions were not been included in this implementation, however, these can be applied by overlaying standard Abaqus elements on the user elements (to be discussed in the **Modeling in Abaqus** section). Since Abaqus/Viewer does not provide native support for visualizing user elements, an additional layer of elements with the same element connectivity has been created and results at the integration points of the elements are stored using the `UVARM` subroutine. |
| `<some_module>.for` | These are the utility files with different Fortran module that are included in the main source file using `include <filename.ext>` statement at the beginning of the main source code. |
| `addElemMech.py` | is a Python code in the `tests` directory that modifies a simple Abaqus input file and adds the overlaying dummy elements on the user elements. For complicated input files, this will not work properly and modification of this code will be required (optional). |
| `<...>.inp` | are the example input files prepared to be executed with the user element subroutine. Since the user-defined elements share the same topology as one of the Abaqus built-in elements, those models were built in Abaqus/CAE and then exported as input files. Later those input files were modified to include keywords and data to include user element definitions, properties, and overlaying dummy elements. |
| `abaqus_v6.env` | is the Abaqus environment file which adds the additional compiling option for the Intel oneMKL library. This needs to be in the same directory as the Abaqus jobs. |
| `runAbq.ps1` | is a PowerShell batch file in the `tests` directory that can execute the user subroutine and specified input file from the PowerShell terminal (optional). |
| elastic_elem.pdf | is a summary of the theory and algorithm used to implement the provided source code. |
| Abaqus_docs.pdf | is a collection of publicly available Abaqus documentation in PDF format related to the Abaqus UEL subroutine. The web versions of these documents are available at https://help.3ds.com. |



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

To visualize the results, an additional set of built-in Abaqus elements with the same element connectivity as the user element has been created in the input file. These additional elements (so-called dummy elements) have negligible elastic properties and thus will not affect the results. If you are using a reduced integration element from the user subroutine, then use the analogous elements from Abaqus as the dummy elements.



## Configuring Abaqus and executing the subroutine

To run user subroutines in Abaqus, you will need to install Intel Visual Studio and Intel oneAPI package and link them with Abaqus. Follow [this blog tutorial](https://www.bibekanandadatta.com/blog/2021/link-intel-and-vs-abaqus-2020/) if you have not done it before. Additionally, [see this blog post](https://www.bibekanandadatta.com/blog/2024/lapack-Intel-Fortran-Abaqus/) to learn how to link and use LAPACK library from tge Intel oneMKL library to Abaqus user subroutines.

Navigate to the `tests` directory. Open the `Abaqus command line terminal` or `cmd terminal` or `PowerShell terminal`, you can execute the following command from the directory to execute the subroutine. Make sure to use the right directory for the main source code.

```bash
abaqus interactive double analysis job=<your_job_name> input=<input_file_name.inp> user=../src/uel_mech.for
```
Specify the variable names (inside < >) in the above command as needed. For additional information on executing user subroutines, check the Abaqus user manual.

If you use the PowerShell-based terminal, you can also execute the subroutine by running the `runAbq.ps1` file. Make sure to check the input file name in the file.
```bash
./runAbq
```



## Citation

If you use this repository (documentation or source code), please consider citing this from the following:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11075088.svg)](https://doi.org/10.5281/zenodo.11075088)

APA format:
```
Datta, B. (2024, April 28). An Abaqus user element (UEL) implementation of linear elastostatics. Zenodo. https://doi.org/10.5281/zenodo.11075088.
```

BibTeX:
``` bibtex
@misc{dattaAbaqusUserElement2024,
  author       = {Datta, Bibekananda},
  title        = {{An Abaqus user element (UEL) implementation of linear elastostatics}},
  month        = apr,
  year         = 2024,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.11075088},
  url          = {https://doi.org/10.5281/zenodo.11075088}
}
```
