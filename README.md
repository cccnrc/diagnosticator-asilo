# this is the ASILO dependency of Diagnosticator local app

##### ASILO runs the algorithm that actually performs ACMG annotation and variants prioritization
###### it runs in a docker container which waits for input files to be detected within the shared
###### docker volume and, once detected, launch the analysis on them
###### input files are the outputs of VEP-FILTER Diagnosticator dependency

###### this application basically consists of a script `waiter-launcher-v0.sh` which constantly runs
###### as soon as a file ending in: `*.asilo_input` (created by VEP-FILTERING docker)  is detected
###### it operates the analysis with the filters specified in that input file and after the analysis
###### run it moves the `*.asilo_input` file to `*.asilo_input.output` and creates `ASILO.NEW`
###### which are used by Diagnosticator rq-worker to understand the most recent analysis and results

#### DEVELOPMENT instructions
```
# git init
git add .
git commit -m "diagnosticator asilo v0.2"
git branch -M main
git remote add origin https://github.com/cccnrc/diagnosticator-asilo.git
git push -u origin main
# git reset --soft HEAD~1      #### avoid last commit if not pushed yet

APP_DIR=<your-ASILO-directory>
```
