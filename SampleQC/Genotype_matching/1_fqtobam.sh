#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-2
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12

set -e 

genome=canfam4 #shorthand name for genome (for bam name)
workingdir=/scratch/bmarine2/Triad_mapping
genome_path=/data/CEM/smacklab/genomes/bwa/UU_Cfam_GSD_1.0_ROSY #path to bwa indexed genome
genome_prefix=UU_Cfam_GSD_1.0_ROSY.fa #indexed genome prefix
output_path=${workingdir}/nvidia_fq2bam_mapping/bams #path to directory for alignments
fastq_path=/scratch/nsnyderm/dap_fastqs/triad_fastqs ## path to fastqs, but could be arg1 during job submission
fastq_list=/scratch/bmarine2/Triad_mapping/triad_150_samples_20240517_fromTerra-240528-first10.tsv ## path to fastq list where each entry/row contains metadata for DAP fastqs (SRA or terra table)
## path to dog_id IDs/prefixes (in this case, the dog_id variable will contain dog ID of DAP dog)
dog_ids=(`awk '{print $3}' $fastq_list`)
platforms=(`awk '{print $1}' $fastq_list`) # fastq files downloaded from Terra are named according to their platform IDs
thread=12

for i in "${!platforms[@]}";
do
        dog_id="${dog_ids[i]}"
        platform="${platforms[i]}"

        #echo $dog_id
        #echo $platform
        #echo $fastq_path

        r1=`ls "$fastq_path"/"$platform"_fastq-r1.fastq.gz`
        r2=`ls "$fastq_path"/"$platform"_fastq-r2.fastq.gz`

        echo $r1
        echo $r2
        echo ${genome_path}/${genome_prefix}

        ## run the alignment
        apptainer exec --nv -B ${workingdir},${genome_path},${fastq_path},${PWD} /packages/apps/simg/parabricks-4.0.sif \
        pbrun fq2bam \
        --ref ${genome_path}/${genome_prefix} \
        --in-fq ${r1} ${r2} \
        --read-group-sm ${dog_id} \
        --bwa-options='-K 100000000 -Y' \
        --out-bam ${output_path}/${dog_id}.${genome}.bam

done

         
## calculate samtools flagstat and coverage
module load samtools-1.16-gcc-11.2.0

for i in "${!platforms[@]}";
do
        dog_id="${dog_ids[i]}"
        platform="${platforms[i]}"

        mkdir -p ${output_path}/flagstat
        mkdir -p ${output_path}/samtools_coverage
        mkdir -p ${output_path}/mosdepth_coverage

        samtools flagstat -@ ${thread} ${output_path}/${dog_id}.${genome}.bam > ${output_path}/flagstat/${dog_id}.${genome}.bam.flagstat

        samtools coverage  \
                -o ${output_path}/samtools_coverage/${dog_id}.${genome}.cov \
                ${output_path}/${dog_id}.${genome}.bam

done
