#########
 Bbtools
#########

**************
 Introduction
**************

BBTools is a suite of fast, multithreaded bioinformatics tools designed
for analysis of DNA and RNA sequence data.

**********
 Versions
**********

-  39.00

**********
 Commands
**********

-  addadapters.sh
-  addssu.sh
-  adjusthomopolymers.sh
-  alltoall.sh
-  analyzeaccession.sh
-  analyzegenes.sh
-  analyzesketchresults.sh
-  applyvariants.sh
-  a_sample_mt.sh
-  bbcms.sh
-  bbcountunique.sh
-  bbduk.sh
-  bbest.sh
-  bbfakereads.sh
-  bbmap.sh
-  bbmapskimmer.sh
-  bbmask.sh
-  bbmerge-auto.sh
-  bbmerge.sh
-  bbnorm.sh
-  bbrealign.sh
-  bbrename.sh
-  bbsketch.sh
-  bbsplitpairs.sh
-  bbsplit.sh
-  bbstats.sh
-  bbversion.sh
-  bbwrap.sh
-  bloomfilter.sh
-  calcmem.sh
-  calctruequality.sh
-  callgenes.sh
-  callpeaks.sh
-  callvariants2.sh
-  callvariants.sh
-  clumpify.sh
-  commonkmers.sh
-  comparegff.sh
-  comparesketch.sh
-  comparessu.sh
-  comparevcf.sh
-  consect.sh
-  consensus.sh
-  countbarcodes.sh
-  countgc.sh
-  countsharedlines.sh
-  crossblock.sh
-  crosscontaminate.sh
-  cutgff.sh
-  cutprimers.sh
-  decontaminate.sh
-  dedupe2.sh
-  dedupebymapping.sh
-  dedupe.sh
-  demuxbyname.sh
-  diskbench.sh
-  estherfilter.sh
-  explodetree.sh
-  fetchproks.sh
-  filterassemblysummary.sh
-  filterbarcodes.sh
-  filterbycoverage.sh
-  filterbyname.sh
-  filterbysequence.sh
-  filterbytaxa.sh
-  filterbytile.sh
-  filterlines.sh
-  filterqc.sh
-  filtersam.sh
-  filtersilva.sh
-  filtersubs.sh
-  filtervcf.sh
-  fixgaps.sh
-  fungalrelease.sh
-  fuse.sh
-  gbff2gff.sh
-  getreads.sh
-  gi2ancestors.sh
-  gi2taxid.sh
-  gitable.sh
-  grademerge.sh
-  gradesam.sh
-  icecreamfinder.sh
-  icecreamgrader.sh
-  icecreammaker.sh
-  idmatrix.sh
-  idtree.sh
-  invertkey.sh
-  kapastats.sh
-  kcompress.sh
-  keepbestcopy.sh
-  khist.sh
-  kmercountexact.sh
-  kmercountmulti.sh
-  kmercoverage.sh
-  kmerfilterset.sh
-  kmerlimit2.sh
-  kmerlimit.sh
-  kmerposition.sh
-  kmutate.sh
-  lilypad.sh
-  loadreads.sh
-  loglog.sh
-  makechimeras.sh
-  makecontaminatedgenomes.sh
-  makepolymers.sh
-  mapPacBio.sh
-  matrixtocolumns.sh
-  mergebarcodes.sh
-  mergeOTUs.sh
-  mergepgm.sh
-  mergeribo.sh
-  mergesam.sh
-  mergesketch.sh
-  mergesorted.sh
-  msa.sh
-  mutate.sh
-  muxbyname.sh
-  partition.sh
-  phylip2fasta.sh
-  pileup.sh
-  plotflowcell.sh
-  plotgc.sh
-  postfilter.sh
-  printtime.sh
-  processfrag.sh
-  processhi-c.sh
-  processspeed.sh
-  randomgenome.sh
-  randomreads.sh
-  readlength.sh
-  readqc.sh
-  reducesilva.sh
-  reformatpb.sh
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
-  rqcfilter2.sh
-  rqcfilter.sh
-  runhmm.sh
-  samtoroc.sh
-  seal.sh
-  sendsketch.sh
-  shred.sh
-  shrinkaccession.sh
-  shuffle2.sh
-  shuffle.sh
-  sketchblacklist2.sh
-  sketchblacklist.sh
-  sketch.sh
-  sortbyname.sh
-  splitbytaxa.sh
-  splitnextera.sh
-  splitribo.sh
-  splitsam4way.sh
-  splitsam6way.sh
-  splitsam.sh
-  stats.sh
-  statswrapper.sh
-  streamsam.sh
-  subsketch.sh
-  summarizecontam.sh
-  summarizecoverage.sh
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
-  unzip.sh
-  vcf2gff.sh
-  webcheck.sh
-  Xcalcmem.sh

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
   #SBATCH --job-name=bbtools     # Job name
   #SBATCH --mail-type=FAIL,BEGIN,END     # Send an email when job fails, begins, and finishes
   #SBATCH --mail-user=your.email@tufts.edu       # Email address for notifications
   #SBATCH --error=%x-%J-%u.err   # Standard error file: <job_name>-<job_id>-<username>.err
   #SBATCH --output=%x-%J-%u.out  # Standard output file: <job_name>-<job_id>-<username>.out

   module purge   ### Optional, but highly recommended.
   module load bbtools/XXXX       ### Latest version is recommended.
