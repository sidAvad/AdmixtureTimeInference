#$ -S /bin/bash      # Use bash
#$ -N jupyter-notebook   # name for job in qstat/output filename prefix
#$ -j y              # output from stdout/err in one file
#$ -m ae             # send email on abort or exit of job
#$ -M sna53@cornell.edu   # email address to send to
#$ -l h_vmem=100G     # use 12G of memory instead of the default (4G)
#$ -q long_term.q@cbsubscb09 # use the long_term queue (for jobs > 24 hours)

##$ -pe bscb 4        # number of cores to use
##$ -binding linear:1 # use hyperthreading (never hurts) -- set # same as -pe #
##$ -q gpu.q
#### For task arrays:
##$ -t 1-40           # Run task numbers 1-8, task variable is $SGE_TASK_ID
####$ -tc 5          # number of tasks to run concurrently (if heavy I/O)


/programs/bin/labutils/mount_server cbsubscb09 /storage
#cd siddharth/AdmixtureTimeInference/

# get tunneling info
XDG_RUNTIME_DIR=""
port=$(shuf -i8000-9999 -n1)
node=$(hostname -s)
user=$(whoami)
cluster=$(hostname -f)

date

# print tunneling instructions jupyter-log
echo -e "
MacOS or linux terminal command to create your ssh tunnel:
ssh -N -L ${port}:${node}:${port} ${user}@${cluster}.tc.cornell.edu
   
For more info and how to connect from windows, 
   see research.computing.yale.edu/jupyter-nb
Here is the MobaXterm info:

Forwarded port:same as remote port
Remote server: ${node}
Remote port: ${port}
SSH server: ${cluster}.tc.cornell.edu
SSH login: $user
SSH port: 22

Use a Browser on your local machine to go to:
localhost:${port}  (prefix w/ https:// if using password)
"

# load modules or conda environments here
# e.g. farnam:
# module load Python/2.7.13-foss-2016b 
source activate root

# DON'T USE ADDRESS BELOW. 
# DO USE TOKEN BELOW
jupyter lab --no-browser --port=${port} --ip=${node}

