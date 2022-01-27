# Diagnosticator-ASILO

#### This is the ASILO dependency of Diagnosticator local app

ASILO runs the algorithm that actually performs ACMG annotation and variants prioritization
it runs in a docker container which waits for input files to be detected within the shared
docker volume and, once detected, launch the analysis on them
input files are the outputs of [VEP-filter] (https://github.com/cccnrc/diagnosticator-VEP-filter) Diagnosticator dependency

this application basically consists of a script `waiter-launcher-v0.sh` which constantly runs
as soon as a file ending in: `*.asilo_input` (created by [VEP-filter] (https://github.com/cccnrc/diagnosticator-VEP-filter) docker)  is detected
it operates the analysis with the filters specified in that input file and after the analysis
run it moves the `*.asilo_input` file to `*.asilo_input.output` and creates `ASILO.NEW`
which are used by Diagnosticator rq-worker to understand the most recent analysis and results



### DEVELOPMENT instructions
```
### clone github repo
git clone https://github.com/cccnrc/diagnosticator-asilo.git
APP_DIR=$( realpath ./diagnosticator-asilo )
cd $APP_DIR

### change whatever in the application files
atom ./waiter-launcher-v0.sh

### rebuild the docker
docker build -t asilo:0.3 .

### try to run the docker (mount a volume in which you have a VCF file)
#     - you can use VCF-EXAMPLE directory in this repo
#     - you also have kidney_acmg59.gl genelist file here
docker run --rm -it --name DX-ASILO \
  -v ${APP_DIR}/VCF-EXAMPLE:/INPUT \
  asilo:0.3 /bin/bash

### from another terminal window create the *.asilo_input file to start the analysis
docker exec -it DX-ASILO /bin/bash
echo -e "DIAGNOSTICATOR-TUTORIAL.vcf\tkidney_acmg59.gl\t10E-5" > /INPUT/try0.asilo_input

### you will see analysis logs in the first terminal and
#     if everything works fine the script will generate:
#      - DIAGNOSTICATOR-TUTORIAL.csv    ### converted CSV of the input VCF
#      - analisi_result                 ### DIR with the results to be loaded in Diagnosticator
#      - ASILO.NEW                      ### flag to store datetime of the last analysis
#      - try0.asilo_input.output        ### from try0.asilo_input
```

### UPDATE github with your new branch
```
cd $APP_DIR
git branch <your-name>-development
git checkout <your-name>-development
git add .
git commit -m "<your-name>-development ..."
git push https://github.com/cccnrc/diagnosticator-asilo.git <your-name>-development
```
