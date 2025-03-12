#######
 Bbmap
#######

**************
 Introduction
**************

Bbmap is a short read aligner, as well as various other bioinformatic
tools.

**********
 Versions
**********

-  38.93
-  38.96

**********
 Commands
**********

-  addadapters.sh
-  a_sample_mt.sh
-  bbcountunique.sh
-  bbduk.sh
-  bbest.sh
-  bbfakereads.sh
-  bbmap.sh
-  bbmapskimmer.sh
-  bbmask.sh
-  bbmerge-auto.sh
-  bbmergegapped.sh
-  bbmerge.sh
-  bbnorm.sh
-  bbqc.sh
-  bbrealign.sh
-  bbrename.sh
-  bbsketch.sh
-  bbsplitpairs.sh
-  bbsplit.sh
-  bbstats.sh
-  bbversion.sh
-  bbwrap.sh
-  calcmem.sh
-  calctruequality.sh
-  callpeaks.sh
-  callvariants2.sh
-  callvariants.sh
-  clumpify.sh
-  commonkmers.sh
-  comparesketch.sh
-  comparevcf.sh
-  consect.sh
-  countbarcodes.sh
-  countgc.sh
-  countsharedlines.sh
-  crossblock.sh
-  crosscontaminate.sh
-  cutprimers.sh
-  decontaminate.sh
-  dedupe2.sh
-  dedupebymapping.sh
-  dedupe.sh
-  demuxbyname.sh
-  diskbench.sh
-  estherfilter.sh
-  explodetree.sh
-  filterassemblysummary.sh
-  filterbarcodes.sh
-  filterbycoverage.sh
-  filterbyname.sh
-  filterbysequence.sh
-  filterbytaxa.sh
-  filterbytile.sh
-  filterlines.sh
-  filtersam.sh
-  filtersubs.sh
-  filtervcf.sh
-  fungalrelease.sh
-  fuse.sh
-  getreads.sh
-  gi2ancestors.sh
-  gi2taxid.sh
-  gitable.sh
-  grademerge.sh
-  gradesam.sh
-  idmatrix.sh
-  idtree.sh
-  invertkey.sh
-  kcompress.sh
-  khist.sh
-  kmercountexact.sh
-  kmercountmulti.sh
-  kmercoverage.sh
-  loadreads.sh
-  loglog.sh
-  makechimeras.sh
-  makecontaminatedgenomes.sh
-  makepolymers.sh
-  mapPacBio.sh
-  matrixtocolumns.sh
-  mergebarcodes.sh
-  mergeOTUs.sh
-  mergesam.sh
-  msa.sh
-  mutate.sh
-  muxbyname.sh
-  normandcorrectwrapper.sh
-  partition.sh
-  phylip2fasta.sh
-  pileup.sh
-  plotgc.sh
-  postfilter.sh
-  printtime.sh
-  processfrag.sh
-  processspeed.sh
-  randomreads.sh
-  readlength.sh
-  reducesilva.sh
-  reformat.sh
-  removebadbarcodes.sh
-  removecatdogmousehuman.sh
-  removehuman2.sh
-  removehuman.sh
-  removemicrobes.sh
-  removesmartbell.sh
-  renameimg.sh
-  rename.sh
-  repair.sh
-  replaceheaders.sh
-  representative.sh
-  rqcfilter.sh
-  samtoroc.sh
-  seal.sh
-  sendsketch.sh
-  shred.sh
-  shrinkaccession.sh
-  shuffle.sh
-  sketchblacklist.sh
-  sketch.sh
-  sortbyname.sh
-  splitbytaxa.sh
-  splitnextera.sh
-  splitsam4way.sh
-  splitsam6way.sh
-  splitsam.sh
-  stats.sh
-  statswrapper.sh
-  streamsam.sh
-  summarizecrossblock.sh
-  summarizemerge.sh
-  summarizequast.sh
-  summarizescafstats.sh
-  summarizeseal.sh
-  summarizesketch.sh
-  synthmda.sh
-  tadpipe.sh
-  tadpole.sh
-  tadwrapper.sh
-  taxonomy.sh
-  taxserver.sh
-  taxsize.sh
-  taxtree.sh
-  testfilesystem.sh
-  testformat2.sh
-  testformat.sh
-  tetramerfreq.sh
-  textfile.sh
-  translate6frames.sh
-  unicode2ascii.sh
-  webcheck.sh

*************
 Example job
*************

Adjust slurm options based on job requirements (`slurm cheat sheet
<https://slurm.schedmd.com/pdfs/summary.pdf>`_):

.. code::

   #!/bin/bash
   #SBATCH -p partitionName  # batch, gpu, preempt, mpi or your group's own partition
   #SBATCH -t 1:00:00  # Runtime limit (D-HH:MM:SS)
   #SBATCH -N 1   # Number of nodes
   #SBATCH -n 1   # Number of tasks per node
   #SBATCH -c 4   # Number of CPU cores per task
   #SBATCH --mem=8G       # Memory required per node
   #SBATCH --job-name=bbmap       # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load bbmap/XXXX ### Latest version is recommended.
