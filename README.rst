=====================================================
Verification and validation for air/water flow models
=====================================================

https://github.com/erdc/air-water-vv

Organization of the test set
----------------------------

The tests are divided into two-dimensional and three-dimensional
tests, with additional subdivisions of those directores as needed. A
test consists of a single directory including

- A README file giving a brief overview of the problem and references
  to published work
- All data files describing the geometry and physical parameters of
  the problem
- A set of input files for a code
  (e.g. https://github.com/erdc/proteus)
- A script for results postprocessing (Optional)
- A script for running regression tests

Checklist for adding new test problems
------------------------

- Clone or fork this repository
- Create a new branch
- Add a new subdirectory with the contents of your test case
- Fit the new test case in the existing folder structure, by creating a new category if applicable
- Document the case by creating a README file, and preferably include appropriate sketches
- Make your case modular by including in the context any parameters subject to change
- Carefully document the context parameters and other parts of the code
- Create a test to check that the case runs without errors 
- Create test to check that results compare well with some reference results (e.g. experiments)
- Include the path of the test script in .travis.yaml file
- Push the branch to the remote repository
- Make a pull request, preferebly by recommending a reviewer
- The reviewer will ensure that all steps above have been followed
- Once accepted, include your name in contributors

Contributors
------------
- Chris Kees, Aggelos Dimakopoulos, Eleni Zve, Matthew Farthing, Aron Ahmadia, Ido Akkerman, Matt Malej, Roham Bakhtyar, Giovanni Cozzuto, Tristan de Lataillade, Giorgos Makrygiannis, Louis Maurel


