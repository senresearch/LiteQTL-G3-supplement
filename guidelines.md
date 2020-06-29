# Structure of the supplementary material for reproducible research project.

### Directory structure:
- code
  - Makefile (This file include commands that runs everything to generate analysis result, figures, tables etc. Please see below how to write your Makefile)
  - compare_time.jl
  - preprocess-data.R
  - generate_figure.jl
  - ...
- data
  - raw data
    - .gitignore (don't need to track data)
  - processed data
    - .gitignore (don't need to track data)
- figure
  - .gitignore (don't want to track images)

README.md ( what is in the repo and how to get to results shown in paper)

### How to write your Makefile.
