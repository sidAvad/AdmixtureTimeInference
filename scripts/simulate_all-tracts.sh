#$ -S /bin/bash      # Use bash
#$ -N simulate_all-tracts  # name for job in qstat/output filename prefix
#$ -o /fs/cbsubscb09/storage/siddharth/AdmixtureTimeInference/simulate_all-tracts_2gens.$JOB_ID
#$ -j y              # output from stdout/err in one file
#$ -m ae             # send email on abort or exit of job
#$ -l h_vmem=4G     # use 256G of memory instead of the default (4G)
#$ -q regular.q # use the regular queue (for jobs < 24 hours)


/programs/bin/labutils/mount_server cbsubscb09 /storage
cd siddharth/AdmixtureTimeInference/

t=2
ntypes=4
mode=all

RES_DIR=~/siddharth/AdmixtureTimeInference/Results/all-tract-lengths
IN_DIR=~/siddharth/AdmixtureTimeInference/Output/ped-sim
pedsim_postfix=poisson_sexavg

#Subset input bp file to get focal samples and write out file 
BP_FILE=${IN_DIR}/admixture-time_"$t"gens_${pedsim_postfix}
((gens=$t+1)) 

grep "g$gens" "$BP_FILE".bp > "$BP_FILE"-focal.bp

#Convert input bp file to genetic positions and write outfile 
python AdmixtureTimeInference/phys2gen.py -i "$BP_FILE"-focal.bp 


for ((typ=1;typ<=$ntypes;typ++))
do
    echo "number of simulated generations = ${t}"
    echo "admixture type : ${typ}"

    #Get homozygous tract lenghts for current admixture type and write output to Results 
    python AdmixtureTimeInference/get_homozyg_tracts.py -i "$BP_FILE"-focal.gen -a "$typ" -l Resources/${t}_admixturelabels.txt -o ${RES_DIR}/${t}gens/admixture-time_"$t"gens_type${typ}_${pedsim_postfix}.res -m ${mode}
done    
