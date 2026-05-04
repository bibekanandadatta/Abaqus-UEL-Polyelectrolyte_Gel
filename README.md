# Abaqus-UEL-Polyelectrolyte_Gel

Abaqus/Standard user element (UEL) subroutine for coupled electro-chemo-mechanics of polyelectrolyte hydrogels.

## Obtaining the repository

If you have `git` installed, you can clone the repository to your local machine using

```bash
git clone https://github.com/bibekananda-datta/Abaqus-UEL-Polyelectrolyte_Gel.git
```

You can also fork the repository and sync your fork as updates are deployed. For testing and development, create a separate branch.

Alternatively, you can download the repository as a `.zip` file using the **Code** drop-down menu on the top-right corner of the GitHub repository page. Note that if you download the repository as a `.zip` file, you will not automatically receive future bug fixes or updates.

## Requirements to execute the code

**Abaqus:** Executing the user element (UEL) subroutine for polyelectrolyte hydrogels requires Abaqus/Standard. Compiling and linking Abaqus user subroutines requires a research or commercial Abaqus license; student or teaching licenses do not provide this capability.

**Intel oneAPI:** To compile and link the code, users need to install Intel oneAPI, including both the Base Toolkit and the HPC Toolkit. The Intel Fortran compiler (`ifort` or `ifx`) is included in the Intel oneAPI HPC Toolkit. The code also relies on LAPACK subroutines from the Intel oneAPI Math Kernel Library (oneMKL), which is included in the Base Toolkit.

**Microsoft Visual Studio:** Microsoft Visual Studio is required to link the compiled object files with Abaqus on Windows.

## Configuring Abaqus and executing the subroutine

To run user subroutines in Abaqus, you need to install Microsoft Visual Studio and Intel oneAPI and link them with Abaqus. Follow [this blog tutorial](https://www.bibekanandadatta.com/blog/2021/link-intel-and-vs-abaqus-2020/) if you have not done this before. Additionally, see [this blog post](https://www.bibekanandadatta.com/blog/2024/lapack-Intel-Fortran-Abaqus/) to learn how to link and use LAPACK from Intel oneMKL in Abaqus user subroutines.

Make sure that the user subroutine and input file are in the appropriate working directory. Using the Abaqus command-line terminal, `cmd`, or PowerShell, you can execute the subroutine using

```bash
abaqus interactive double analysis ask_delete=off job=<your_job_name> input=<input_file_name.inp> user=../src/uel_pegel.for
```

Specify the variable names inside `< >` as needed. For additional information on executing user subroutines, consult the Abaqus user manual.

> [!WARNING]
> Starting in 2025, Intel discontinued the classic Fortran compiler (`ifort`), and newer Intel oneAPI releases contain the newer `ifx` compiler. This code has not yet been tested with newer versions of the Intel compiler. As of Abaqus 2025, official support for the `ifx` compiler is not available. However, the installation, configuration, and execution procedure is expected to remain similar.

## Available features

### Elements

The implemented code currently supports the following elements. Following Abaqus convention, user elements have the prefix `U` followed by a numerical tag.

- `U1`: Three-dimensional 4-node tetrahedral element.
- `U2`: Three-dimensional 8-node hexahedral element.
- `U3`: Two-dimensional 3-node triangular axisymmetric element.
- `U4`: Two-dimensional 4-node quadrilateral axisymmetric element.
- `U5`: Two-dimensional 3-node triangular plane strain element.
- `U6`: Two-dimensional 4-node quadrilateral plane strain element.

The element formulations are implemented in the `pegel_element` module. General 3D elements (`U1`–`U2`) and 2D plane strain elements (`U5`–`U6`) are implemented in the `pegel_general` subroutine. Axisymmetric elements (`U3`–`U4`) are implemented in the `pegel_axisymmetric` subroutine.

> [!NOTE]
> Plane stress elements are not available.

### Material model

The `pegel_material` module contains the subroutine implementing the constitutive laws and their tangents for polyelectrolyte hydrogels. Currently, the polyelectrolyte gel constitutive model includes a pre-swollen Neo-Hookean elastomer, Flory-Huggins mixing potential, and dilute ionic species mixture.

### Error logging

Errors and warnings generated during the computation are written to an `<aaErr_job_name>.dat` file. This file is opened using the `UEXTERNALDB` subroutine provided by Abaqus.

### Visualization

To visualize the results, an additional set of built-in Abaqus elements with the same connectivity as the user elements is created in the input file. These additional elements, also called dummy elements, have negligible elastic properties and therefore do not affect the results. Post-processing variables are transferred to these dummy elements using the Abaqus `UVARM` subroutine. Users need to define the variables for post-processing in the input file.

## Brief description of the repository

### `src`

This directory contains the Fortran source code for the user element implementation.

- **Driver:** `uel_pegel.for` is the main driver file containing the Abaqus-provided interfaces for the `UEL`, `UVARM`, and `UEXTERNALDB` subroutines. Abaqus requires a single user subroutine file to be provided during compilation. The remaining source files are included in the main file using `include <filename.for>`.

- **Modules:** The main driver uses the following modules:

  - `global_parameters.for`: Sets global variables and constants accessed by other modules and subroutines.
  - `error_logging.for`: Contains subroutines for writing error and warning messages to the screen and to an error-log `.dat` file.
  - `post_processing.for`: Sets global variables required for post-processing of user element results through the `UVARM` subroutine.
  - `lagrange_element.for`: Contains subroutines for evaluating Lagrangian shape functions used in element-level operations.
  - `gauss_quadrature.for`: Performs numerical integration of the element tangent stiffness matrix and residual vector at each integration point.
  - `surface_integration.for`: Contains subroutines for surface integration used to implement traction- and pressure-type loading. These routines are not used in the current implementation.
  - `nonlinear_solver.for`: Contains Newton-Raphson-based nonlinear solver subroutines for single-equation and multiple-equation systems.
  - `linear_algebra.for`: Contains subroutines for linear algebra operations.
  - `solid_mechanics.for`: Contains subroutines for managing vectors and tensors following the conventions used in this code.
  - `pegel_element.for`: Contains subroutines for evaluating element formulations for general 2D, 3D, and axisymmetric elements.
  - `pegel_material.for`: Contains subroutines for evaluating the material constitutive behavior and the tangents required to compute the element tangent matrix and residual vector.

### `test`

This directory contains Abaqus input files (`.inp`) for several examples presented in the published article. The following three subdirectories contain the specific examples presented in the manuscript:

- `sun_calibration_validation`
- `bilayer_bending`
- `confined_compression`

Although parametric studies were performed for the paper, those input files are not included because they were generated from the provided baseline cases.

### `utils`

This directory contains a PowerShell script to run Abaqus jobs from a PowerShell window and print the status file in a separate PowerShell window or tab. It also contains a simple Python script for adding dummy elements for visualization in a simple mesh.

### `dev`

This directory contains the same source code as the `src` directory, along with the additional source file `test_uel_pegel.for` for numerical testing of the implementation.

Currently, numerical testing is available for axisymmetric and plane strain elements.

To compile the code using the Intel Fortran compiler, use

```bash
ifort /Qmkl -o pegel test_uel_pegel.for
```

To run the code, use

```bash
.\pegel
```

### `docs`

This directory contains the preprint of the article and supplementary material. The supplementary material provides a brief description of the different modules and detailed user guidelines on how to prepare the input file and use the code in Abaqus/Standard. For full Abaqus documentation, visit <https://help.3ds.com>. You may need to create an account to access the online documentation.
