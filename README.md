# TODO: fix this readme, explain what's done etc
# Genesis-DB: a database for autonomous laboratory systems
This repo contains code for the demonstration in the paper Genesis-DB: a database for autonomous laboratory systems

## Install python dependencies
Using Python 3.8.13 the dependencies can be installed from the `requirements.txt` file, e.g. using conda and the following commands:
```
$ conda create --name db-demo python=3.8.13 && \\
    conda activate db-demo && \\
    pip install -r requirements.txt
```

## Generate RDF representation of data
From the `scripts` directory, run
```
$ python ccp_ttl.py
```
to generate the file `ccp.ttl`, modelling the experimental conditions, in the `data` directory. Also run the following command
```
$ python trans_ttl.py
```
to generate the transcirptomics experimental results in the files `transcriptomics<0-4>.py` in the `data` directory. Note, the generated data will include the high temperature experiments.

## Add data to database
Clone [this](https://github.com/TW-Genesis/genesis-database-system) repository and copy the `.ttl` files generated above to the `data` directory in `genesis-database-system`. Then follow the instructions to build the docker image, load the data, and start the database server. 

## Get gene counts
From the `scripts` directory, run
```
$ python csv_from_query.py
```

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
