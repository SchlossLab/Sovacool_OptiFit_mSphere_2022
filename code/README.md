Before submitting a job to the cluster, make sure the email, account, etc. are correct.

slurm:
```
#SBATCH --account=your_account
#SBATCH --mail-user=your_email
```

pbs:
```
#PBS -A your_allocation
#PBS -M uniqname@umich.edu
```
