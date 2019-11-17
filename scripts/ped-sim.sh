#$ -S /bin/bash      # Use bash
#$ -N ped-sim_admixture-times  # name for job in qstat/output filename prefix
#$ -o /fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/$JOB_ID.ped-sim
#$ -j y              # output from stdout/err in one file
#$ -m ae             # send email on abort or exit of job
#$ -l h_vmem=4G     # use 256G of memory instead of the default (4G)
#$ -q regular.q # use the regular queue (for jobs > 24 hours)

/programs/bin/labutils/mount_server cbsubscb09 /storage
cd siddharth/AdmixtureTimeInference/

OUT_DIR=~/siddharth/AdmixtureTimeInference/Output/ped-sim
MAP_SPF=/fs/cbsubscb09/storage/resources/genetic_maps/refined_mf.simmap #mapfile for ped-sim run 
MAP_AVG=~/siddharth/AdmixtureTimeInference/Resources/sex-averaged_mf.simmap 

t=3
ntypes=16
DEF_FILE=~/siddharth/AdmixtureTimeInference/Resources/admixture-times_${t}gens.def

#Run ped-sim with poisson model and sex averaged map
echo 'running ped-sim in with sex averaged maps and poisson crossover model'
echo "${t}gens is being simulatied" && cat ${DEF_FILE}

ped-sim --bp -d ${DEF_FILE} -m ${MAP_AVG} -i /dev/null --pois -o ${OUT_DIR}/admixture-time_${t}gens_poisson_sexavg --founder_ids 

#Run ped-sim with interference model and sex-specific map 
echo 'running ped-sim with sex specific maps and interference crossover model'
echo "${t}gens is being simulated" && cat ${DEF_FILE}

ped-sim --bp -d ${DEF_FILE} -m ${MAP_SPF} -i /dev/null --intf ~/siddharth/programs/ped-sim/interfere/nu_p_campbell.tsv -o ${OUT_DIR}/admixture-time_${t}gens_intf_sexspf --founder_ids 





