# CyclingSignatures.jl

CyclingSignatures.jl is a Julia package to identify and classify oscillations in time series.

**Note** Tutorial and Docs are still work in progress.

## Quickstart Guide

Download and install Julia. To view the examples, download the repository, open the main folder in VSCode, execute in the Pkg-Repl

    pkg> activate .
    pkg> instantiate

and run either of examples/PaperExamples/lorenz-plotting.jl, examples/PaperExamples/dw-plotting.jl, examples/PaperExamples/dadras-plotting.jl line by line. Important: The .jld2 files are stored using git lfs. To download them using git, git lfs needs to be installed.

To run the experiments, run examples/PaperExamples/lorenz-experiment.jl, examples/PaperExamples/double-well-experiment.jl or examples/PaperExamples/dadras-experiment.jl either as scripts or line by line in VSCode. They will overwrite the existing xyz-results.jld2 files and can then be viewed using the plot files.
