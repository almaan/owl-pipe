#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1

#SBATCH -t 12:00:00
#SBATCH -J aa-pipe
#SBATCH --mail-user alma.andersson@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH -e /home/satac/slurm/job-%J.err
#SBATCH -o /home/satac/slurm/job-%J.out


module load python/3.6.4
module load samtools/1.6
module load star/2.5.4a

export PYTHONPATH=""
export PYTHONHOME=""

#set -e

source /home/alma.andersson/.bashrc

INPUTDIR="/home/alma.andersson/satac/data/20201130-satac/20201130-bcl2fastq"
OUTDIR="/home/alma.andersson/satac/data/20201130-satac/PROC"

# Make sure to specify correct genome

#mouse
MAP=/fastdisk/mouse/GRCm38_86v2/StarIndex
ANN=/fastdisk/mouse/GRCm38_86v2/annotation/gencode.vM11.annotation.gtf

#human
#MAP=/fastdisk/human/GRCh38_86v2/StarIndex
#ANN=/fastdisk/human/GRCh38_86v2/annotation/gencode.v25.annotation.gtf

# Barcodes settings
ID="/home/satac/rsc/omni-v1-barcode-file.txt"

# INDEX LIST
IDX="/home/satac/data/test-pipe/index-names.txt"

cd /home/satac/scripts/my_pipe

python3 pipeline.py  -r1 $INPUTDIR/R1.fastq.gz\
    -r2 $INPUTDIR/R2.fastq.gz\
    -o $OUTDIR\
    -bc $ID\
    -in $IDX\
    -ow\
    -gd $MAP\
    -an $ANN
