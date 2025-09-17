thesis
==============

This repository contains all source files and supporting code for my PhD thesis. It includes:
- **LaTeX** files---complete thesis source, including chapters, appendices, and references.
- **C++ programs**, **Mathematica notebooks**, and **FORM scripts**---I used to validate analytic expressions.
- **Plotting utilities**---scripts and programs for generating the figures and plots appearing in the thesis.

Feel free to use, modify, and redistribute for your own projects.

## Prerequisites

- A complete `LaTeX` distribution (e.g., texlive).

- The `make` build automation utility.

## Workflow

The suggested workflow is:

1. Clone this repository:

    `git clone git@github.com:tomtom1-4/phd-thesis.git`.

2. Create a new branch:

    `git checkout -b mythesis`.

3. If you want to build the thesis, go into the `LaTeX` directory

    `cd LaTeX`

3. First-time complete build of the thesis:

    `make -s thesis`.

4. Utilities for the generation of the figures can be found in the `figures` subdirectory.

5. Programs I used to validate analytic expressions can be found in the `programs` subdirectory.

## Commands

* Remove all the generated files except the PDFs:

    `make -s clean`.

* Remove all the generated files:

    `make -s cleanall`.

* Compile the thesis:

    `make -s thesis`.

    This command embeds the correct sequence of commands relevant for the proper compilation of the thesis, including analytical index, bibliography, and title page.

* Update the analytical index:

    `make -s index && make -s tex`.

## Credits

* [Lorenzo Pantieri](http://www.lorenzopantieri.net),

    author of *thesis* package which is aimed to improve some ty­po­graph­i­cal points of the *Clas­sicTh­e­sis* style.

    Documentation and features explaination ca be found [here](http://ftp.uniroma2.it/TeX/macros/latex/contrib/thesis/thesis.pdf).

* [André Miede](http://www.ctan.org/author/miede),

    author of *ClassicThesis* package ([link](http://ctan.mirror.garr.it/mirrors/CTAN/macros/latex/contrib/classicthesis/ClassicThesis.pdf)) which pro­vides an el­e­gant lay­out de­signed in homage to Bringhurst’s "The Ele­ments of Ty­po­graphic Style".

* [Carlo Tasillo](https://www.linkedin.com/in/dr-carlo-tasillo-7bb98a266/),

    made me aware of the above packages, and I used some of his modifications in [thesis-template-DESY](https://github.com/tasicarl/thesis-template-DESY) as inspiration for my own project.
