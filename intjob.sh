#! /bin/sh
qsub -I -q skystd -l ncpus=4:mem=4gb,walltime=2:00:00
module load anaconda3/5.1.0-gcc/8.3.1
source activate mda
python inputwriter.py