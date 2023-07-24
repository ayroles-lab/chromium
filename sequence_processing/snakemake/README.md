snakemake --use-conda --cluster-config cluster_config.json --cluster "sbatch --job-name={cluster.job-name} --output={cluster.output} --error={cluster.error} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --time={cluster.time} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user}" --jobs 400 --conda-frontend conda --rerun-incomplete --latency-wait 240


squeue -u $USER | awk '{print $1}' | xargs -n 1 scancel

snakemake --delete-all-output -c24  --rerun-incomplete

rm logs/slurm/*