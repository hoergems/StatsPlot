#!/bin/sh
#
#SBATCH --job-name=plot_job
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4096
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hoergems@gmail.com

source /home/hoe01h/.bash_profile
cd /data/hoe01h/Downloads/StatsPlot/
python plot_stats.py -d /datastore/hoe01h/dubin/randomScene/0 -s -cf
python plot_stats.py -d /datastore/hoe01h/dubin/randomScene/10 -s -cf
python plot_stats.py -d /datastore/hoe01h/dubin/randomScene/20 -s -cf
python plot_stats.py -d /datastore/hoe01h/dubin/randomScene/30 -s -cf
