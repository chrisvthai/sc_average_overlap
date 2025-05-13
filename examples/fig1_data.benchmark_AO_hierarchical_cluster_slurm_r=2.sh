#!/bin/bash
#SBATCH --job-name=BenchmarkAOHierarchical_r=2
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=p_foran
#SBATCH --mem=128000
#SBATCH --export=ALL
#SBATCH --output=slurm.%N.%j.out

echo "Start: $(date)"

echo "Resolution: 2.0"
python benchmark_AO_hierarchical_cluster.py --resolution=2.0 --n_markers=100 --n_iters=100
python benchmark_AO_hierarchical_cluster.py --resolution=2.0 --n_markers=50 --n_iters=100
python benchmark_AO_hierarchical_cluster.py --resolution=2.0 --n_markers=25 --n_iters=100
python benchmark_AO_hierarchical_cluster.py --resolution=2.0 --n_markers=10 --n_iters=100
python benchmark_AO_hierarchical_cluster.py --resolution=2.0 --n_markers=5 --n_iters=100

echo "End: $(date)"
