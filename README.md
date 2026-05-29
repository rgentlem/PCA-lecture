# PCA-lecture
A lecture on using PCA in Computational Biology

## Toronto 2026 CBHC version

`PCA2_quarto.qmd` is the Quarto source for the version of this lecture presented in Toronto in 2026 at the CBHC conference.

Rendered outputs and intermediate build artifacts, such as HTML, PowerPoint, and `_files/` directories, are generated from the source document and are not intended to be tracked in the repository.

## Build notes

In order to build the document you need to first run both `SetUpGO.R` and `OSCAmultisample.R`.

The latter takes quite a while to run, and you might instead just edit out those slides from the lecture.

There is a file, `package_versions.csv` that contains a list of all the packages you need. 
They are obtained by running `sessionInfo()` after having sourced both of the other scripts.
