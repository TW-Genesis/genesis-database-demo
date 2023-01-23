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
which will generate a copy of `sce_RNA_RAW_counts.csv` based on the data in the database, which is the input to the hLICORN algorithm.

## Finding the GRN using hLICORN
To run hLICORN algorithm install R (tested with version 4.2.1), e.g. using conda
```
$ conda create --name R && \\
    conda activate R && \\
    conda install -c conda-forge r-base=4.2.1
```
From the `hlicorn` directory, running the following command
```
$ Rscript differentialExpression.R
```
will install dependencies and save regulatory networks to `reg_frames<ignored condition>.Rdata` in the `data` directory. Chose if conditions should be ignored on lines `15-19` in `differentialExpression.Rdata` - ignoring conditions corresponds to missing experiments.

## Visualise GRN
Visualise the GRN by running
```
$ python plot_grn.py
```
from the `scripts` directory. Specify which GRN to visualise using the `cond` flag, e.g. as `--cond SCT' to ignore the high temperature experiment. Make sure the corresponding `reg_frames.Rdata` file has been generated. For more help on options run `python plot_grn.py --help`.

## Visualise experimental conditions
Run
```
$ python plot_cond.py
```
to visualise the explored experimental conditions. Ignore specified conditions using the `ignore_process` flag, e.g. as `--ignore_process SCT`. For more help on options run `python plot_cond.py --help`.
