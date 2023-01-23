# TODO: fix this readme, explain what's done etc
# Genesis-DB: a database for autonomous laboratory systems
This repo contains code for the demonstration in the paper Genesis-DB: a database for autonomous laboratory systems

[fig](figs/paper-fig.pdf)

To run hLICORN algorithm install R (tested with version 4.2.1), e.g. using conda
```
$ conda create --name R && \\
    conda activate R && \\
    conda install -c conda-forge r-base=4.2.1
```
After installing dependencies from `renv.lock`, run algorithm (from the hlicorn directory) with the following command
```
$ Rscript differentialExpression.R
```

For the Python scripts to interface with the database and visualise the GRN, install the dependencies in `requirements.txt` (tested with Python 3.8.13), e.g. with the following commands

```
$ conda create --name db-demo python=3.8.13 && \\
    conda activate db-demo && \\
    pip install -r requirements.txt
```
