# HelloWorld HPC Cluster Tutorial

This simple tutorial will walk you through running your first HPC job using SSH and SLURM job scheduler.

---

## 1. Create `helloworld.py` on your **local computer** using your favorite text editor. TextEdit and Notepad are two good options.

Create a Python script named `helloworld.py` with the following contents:

```python
# helloworld.py
import time

for i in range(10):
    print(f"Hello from HPC! Iteration {i+1}")
    time.sleep(120)  # Sleep for 120 seconds so we can see the job in the queue
```

---

## 2. Upload `helloworld.py` to the cluster

Use `scp`:

```bash
scp helloworld.py your_utln@login.pax.tufts.edu:~/
```

Or upload using **Open OnDemand**, if available.

---

## 3. SSH into the cluster

```bash
ssh your_utln@login.pax.tufts.edu
```

---

## 4. Create a SLURM batch file using `nano` or `vi`

Run:

```bash
nano helloworld_job.sh
```

Paste this content:

```bash
#!/bin/bash
#SBATCH --job-name=hello_job
#SBATCH --output=hello_output.log
#SBATCH --error=hello_error.log
#SBATCH --time=00:02:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load python  # Load Python if required by your system
python helloworld.py
```

Save and exit with `Ctrl+O`, `Enter`, `Ctrl+X`.

---

## 5. Submit your batch job

```bash
sbatch helloworld_job.sh
```

---

## 6. Check job status

```bash
squeue --me
```

---

## 7. Check output and error logs

After the job completes, view logs:

```bash
cat hello_output.log
cat hello_error.log
```

---

## 8. Check job efficiency with `seff`

Find your job ID and run:

```bash
seff JOB_ID
```

Example:

```bash
seff 123456
```

---

✅ You’ve now run your first HPC job using SLURM and SSH!
