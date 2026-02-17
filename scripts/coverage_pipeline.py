#!/usr/bin/env python

import os
import pandas as pd
import sys
import argparse
import pickle
from genome_coverage_checker.util import *
import time

start_time = time.time()

# 0. Get command line arguments
parser = argparse.ArgumentParser(description='This script is to check which taxa reads have been assigned to by Kraken, pull out these reads, download reference genomes for the taxa, and map the reads to the reference genomes.')
parser.add_argument('--processors', dest='n_proc', default=1,
                    help="Number of processors to use for all multiprocessing steps where this is possible")
parser.add_argument('--fastq_dir', dest='fastq_dir', default=None,
                    help="The directory containing the fastq files")
parser.add_argument('--kraken_kreport_dir', dest='kraken_kreport_dir', default=None,
                    help="The directory containing the kraken kreport files")
parser.add_argument('--kraken_outraw_dir', dest='kraken_outraw_dir', default=None,
                    help="The directory containing the kraken outraw files")
parser.add_argument('--output_dir', dest='output_dir', default=None,
                    help="The directory to save the output files to")
parser.add_argument('--genome_dir', dest='genome_dir', default=None,
                    help="The directory to save the downloaded genomes to. The default is for them to be added to the output_dir")
parser.add_argument('--bowtie2_db_dir', dest='bowtie2_db_dir', default=None,
                    help="The directory to save the bowtie2 database to. The default is for them to be added to the output_dir")
parser.add_argument('--read_lim', dest='read_lim', default=None,
                    help="Look at all taxa with > this number of reads mapped by Kraken")
parser.add_argument('--read_mean', dest='read_mean', default=None,
                    help="Look at all taxa with > this number of mean reads mapped by Kraken within sample groupings")
parser.add_argument('--assembly_folder', dest='assembly_folder', default=None,
                    help="The folder containing the assembly summaries for bacteria and archaea")
parser.add_argument('--sample_metadata', dest='sample_metadata', default=None,
                    help="Location of the sample metadata file. It is expected that this will be a CSV (comma separated) file with a header and two columns. The first column contains sample names and the second contains the sample groupings")
parser.add_argument('--species', dest='species', default=None,
                    help="Location of the species file. An optional file containing a list of taxonomy ID's or species names to include")
parser.add_argument('--project_name', dest='project_name', default=None,
                    help="Name of this project")
parser.add_argument('--rerun', dest='rerun', default=False, action='store_true',
                    help="If this flag is added, it will re-run everything regardless of whether the folder and a checkpoint already exist")
parser.add_argument('--all_domains', dest='all_domains', default=False, action='store_true',
                    help="The default for coverage checker is for only the genomes of prokaryotes to be downloaded and checked. If you'd like to include all genomes then add this flag.")
parser.add_argument('--representative_only', dest='representative_only', default=False, action='store_true',
                    help="The default for coverage checker is to use all genomes. If you'd like to limit the checking to only representative and reference genomes then add this flag.")
parser.add_argument('--bowtie2_setting', dest='bowtie2_setting', default='sensitive',
                    help="The default bowtie2 setting to use. Options are very-fast, fast, sensitive or very-sensitive. Default is sensitive (same as Bowtie2 default)")
parser.add_argument('--skip_coverage', dest='skip_coverage', default=False, action='store_true',
                    help="Whether to skip getting coverage for all reads across the genomes. This can take a little while, so you can skip it if you don't think that this output is useful to you.")
parser.add_argument('--skip_cleanup', dest='skip_cleanup', default=False, action='store_true',
                    help="If you want to skip the cleanup. This will keep the intermediate files containing some of the commands run, e.g. for making bowtie2 databases. They may be helpful if you're trying to troubleshoot issues.")
parser.add_argument('--duplicate_check', dest='duplicate_check', default=False, action='store_true',
                    help="If you want to do the check for duplicates within the fastq files. Note that this step can take a while if you have a lot of samples - it was mainly added because you'll get some weird results if you have duplicate reads in your files. This can happen if you rerun coverage checker using the same output folder.")
parser.add_argument('--grouped_samples_only', dest='grouped_samples_only', default=False, action='store_true',
                    help="If you only want to run coverage checker with the grouped samples (i.e. by metadata variable or overall). The default is to run coverage checker individually on each sample, but if you only want the overall results, it will save on computation time to run coverage checker with this option.")
parser.add_argument('--no_grouped_samples', dest='no_grouped_samples', default=False, action='store_true',
                    help="If you only want to run coverage checker on the individual samples. The default is to run coverage checker individually on each sample as well as on the groups.")
parser.add_argument('--coverage_program', dest='coverage_program', default='Bowtie2', choices=['Minimap2', 'Bowtie2', 'Both'],
                    help="Which of the programs to use for getting coverage across the genome. Default is Bowtie2.")
parser.add_argument('--mapq_threshold', dest='mapq_threshold', default=None, choices=range(0,256),
                    help="The threshold to use to determine coverage within the genomes. Note that this is not used if you have set --skip_coverage.")
parser.add_argument('--identity_threshold', dest='identity_threshold', default=None, choices=range(0,101),
                    help="The threshold to use to determine coverage within the genomes. Note that this is not used if you have set --skip_coverage.")


#Read in the command line arguments
args = parser.parse_args()
n_proc, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, genome_dir, bowtie2_db_dir = int(args.n_proc), args.fastq_dir, args.kraken_kreport_dir, args.kraken_outraw_dir, args.output_dir, args.genome_dir, args.bowtie2_db_dir
if '.' in output_dir:
  sys.exit("No . can be in the output_dir name. Please use only _ and -")
if genome_dir == None:
  genome_dir = output_dir+'/genomes/'
if bowtie2_db_dir == None:
  bowtie2_db_dir = output_dir+'/bowtie2_db/'
read_lim, read_mean, assembly_folder = args.read_lim, args.read_mean, args.assembly_folder
if assembly_folder == None:
  assembly_folder = output_dir
sample_metadata, species, project_name, rerun, all_domains, representative_only, skip_coverage, skip_cleanup, duplicate_check, bowtie2_setting, grouped_samples_only, no_grouped_samples, coverage_program, mapq_threshold, identity_threshold = args.sample_metadata, args.species, args.project_name, args.rerun, args.all_domains, args.representative_only, args.skip_coverage, args.skip_cleanup, args.duplicate_check, args.bowtie2_setting, args.grouped_samples_only, args.no_grouped_samples, args.coverage_program, args.mapq_threshold, args.identity_threshold
wd = os.getcwd()
if read_lim == None: read_lim = 0
else: read_lim = int(read_lim)
if read_mean == None: read_mean = 0
else: read_mean = int(read_mean)
if mapq_threshold != None:
  try:
    float(mapq_threshold)
  except:
    sys.exit("mapq_threshold must be a number")
if identity_threshold != None:
  try:
    float(identity_threshold)
  except:
    sys.exit("identity_threshold must be a number")
if project_name == None:
  sys.exit("You must set --project_name")
    
if grouped_samples_only and no_grouped_samples:
  sys.exit("You cannot set both grouped_samples_only and no_grouped_samples. These contradict eachother. Please choose one only and try running again.")
    
if coverage_program in ['Minimap2', 'Both']:
  sys.stdout.write("Please note that we have had issues with running Minimap2 where too many threads are used despite the default to be to use 1 thread for each file (so this should be multiplied by the number of threads that you have set.\n")
  sys.stdout.write("You can stop the run if you are concerned that this may happen or that this could cause you problems.\n")
  sys.stdout.flush()

o = sys.stdout
time_string = str(time.ctime(start_time)).replace(':', '-')
print('Logging all output to '+'Genome Coverage Checker log '+time_string+'.txt\n')
print('Check this file at any point to see where Genome Coverage Checker is in the pipeline.\n')

f = open('Genome Coverage Checker log '+time_string+'.txt', 'w')
sys.stdout = f

#check whether we've already run this and which checkpoint we're at
if rerun:
  cp = '0'
elif os.path.exists(output_dir+'/checkpoint.txt'):
  cp = get_checkpoint(output_dir)
  sys.stdout.write("Got information from previous run: already at "+cp+"\n\n")
  sys.stdout.flush()
else:
  cp = '0'

#if not rerun and cp not 0, check whether the arguments given are the same?
if cp != '0':
  if os.path.exists(output_dir+'/pickle_intermediates/args.pickle'):
    with open(output_dir+'/pickle_intermediates/args.pickle', 'rb') as f:
      all_args = pickle.load(f)
    wd, n_proc, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name, genome_dir, bowtie2_db_dir, all_domains, representative_only, skip_coverage, skip_cleanup, duplicate_check, bowtie2_setting, grouped_samples_only, no_grouped_samples, coverage_program = all_args
    
# 1. Run the initial checks
if cp == '0':
  sys.stdout.write("Running initial checks\n")
  sys.stdout.flush()
  wd, n_proc, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name, genome_dir, bowtie2_db_dir = run_initial_checks(wd, n_proc, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, genome_dir, bowtie2_db_dir, coverage_program, skip_coverage) #run all of the initial checks to ensure that all of the folders and files exist before starting to try and run anything
  all_args = [wd, n_proc, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, output_dir, assembly_folder, read_lim, read_mean, sample_metadata, species, project_name, rerun, md, samples, taxid_name, genome_dir, bowtie2_db_dir, all_domains, representative_only, skip_coverage, skip_cleanup, duplicate_check, bowtie2_setting, grouped_samples_only, no_grouped_samples, coverage_program]
  save_pickle(all_args, output_dir+'/pickle_intermediates/args.pickle')
  cp = update_checkpoint(output_dir, "1_initial_checks_run")
  sys.stdout.write("Completed check-point 1 initial checks\n\n")
  sys.stdout.flush()
else:
  sys.stdout.write("Skipping initial checks because check-point wasn't 0\n\n")
  sys.stdout.flush()
  

# 2. Read in all kreports and combine them to get the genomes that we'll need
if cp == "1_initial_checks_run":
  sys.stdout.write("Reading in kreports and combining them\n")
  sys.stdout.flush()
  group_samples, taxid, kreports = get_kreports(samples, kraken_kreport_dir, output_dir, md, read_lim, read_mean, project_name, taxid_name)
  names, objects = ['group_samples', 'taxid', 'kreports'], [group_samples, taxid, kreports]
  for n in range(len(names)):
    save_pickle(objects[n], output_dir+'/pickle_intermediates/'+names[n]+'.pickle')
  cp = update_checkpoint(output_dir, "2_combined_kreports")
  sys.stdout.write("Completed check-point 2 combined kreports\n\n")
  sys.stdout.flush()
else:
  print("Skipping 2_combined_kreports because check-point wasn't 1_initial_checks_run\n\n")
  names, objects = ['group_samples', 'taxid', 'kreports'], []
  for n in range(len(names)):
    with open(output_dir+'/pickle_intermediates/'+names[n]+'.pickle', 'rb') as f:
      objects.append(pickle.load(f))
  group_samples, taxid, kreports = objects

if grouped_samples_only:
  sys.stdout.write("You are running Genome Coverage Checker using the --grouped_samples_only option, giving %s taxa and %s samples. This gives %s taxon-sample combinations to be run.\n" % (len(taxid), len(group_samples), len(taxid)*len(group_samples)))
  sys.stdout.write("If you think this will take a long time, consider changing the --read_lim or --read_mean options.\n\n")
elif no_grouped_samples:
  sys.stdout.write("You are running Genome Coverage Checker with %s taxa and %s samples. This gives %s taxon-sample combinations to be run.\n" % (len(taxid), len(samples), len(taxid)*len(samples)))
  sys.stdout.write("If you think this will take a long time, consider stopping and re-running with the --grouped_samples_only option, or with a higher --read_lim or --read_mean option set.\n")
  sys.stdout.write("Using the --grouped_samples_only option would give %s samples and %s taxon-sample combinations.\n\n" % (len(group_samples), len(group_samples)*len(taxid)))
else:
  sys.stdout.write("You are running Genome Coverage Checker with %s taxa and %s samples. This gives %s taxon-sample combinations to be run.\n" % (len(taxid), len(group_samples)+len(samples), len(taxid)*(len(group_samples)+len(samples))))
  sys.stdout.write("If you think this will take a long time, consider stopping and re-running with the --grouped_samples_only option, or with a higher --read_lim or --read_mean option set.\n")
  sys.stdout.write("Using the --grouped_samples_only option would give %s samples and %s taxon-sample combinations.\n\n" % (len(group_samples), len(group_samples)*len(taxid)))
sys.stdout.flush()  

# 3. Download all genomes
if cp == "2_combined_kreports":
  sys.stdout.write("Downloading all genomes\n")
  sys.stdout.flush()  
  taxid = download_genomes(taxid, assembly_folder, output_dir, all_domains, representative_only, n_proc, genome_dir)
  save_pickle(taxid, output_dir+'/pickle_intermediates/taxid.pickle')
  cp = update_checkpoint(output_dir, "3_downloaded_genomes")
  sys.stdout.write("Completed check-point 3 downloaded genomes\n\n")
  sys.stdout.flush()  
else:
  sys.stdout.write("Skipping 3_downloaded_genomes because check-point wasn't 2_combined_kreports\n\n")
  sys.stdout.flush()  
  with open(output_dir+'/pickle_intermediates/taxid.pickle', 'rb') as f:
    taxid = pickle.load(f)

# 4. Extract the reads for each taxonomy ID
if cp == "3_downloaded_genomes":
  sys.stdout.write("Extracting reads for all taxonomy ID's\n")
  sys.stdout.flush()
  extract_reads(taxid, output_dir, samples, fastq_dir, kraken_kreport_dir, kraken_outraw_dir, n_proc)
  cp = update_checkpoint(output_dir, "4_extracted_reads")
  sys.stdout.write("Completed check-point 4 extracted all reads\n")
  sys.stdout.write("%s read files exist in the directory\n\n" % (len(os.listdir(output_dir+'/reads_mapped'))))
  sys.stdout.flush()
else:
  sys.stdout.write("Skipping 4_extracted_reads because check-point wasn't 3_downloaded_genomes\n\n")
  sys.stdout.flush()

# 5. Combine the files for each taxonomy ID (across multiple samples for the sample groupings)
if cp == "4_extracted_reads":
  sys.stdout.write("Combining files for each taxonomy ID\n")
  sys.stdout.flush()
  all_files = combine_convert_files_paf(taxid, output_dir, samples, group_samples, n_proc, genome_dir, duplicate_check, grouped_samples_only, no_grouped_samples)
  save_pickle(all_files, output_dir+'/pickle_intermediates/all_files.pickle')
  cp = update_checkpoint(output_dir, "5_combined_files")
  sys.stdout.write("Completed check-point 5 combined files\n")
  sys.stdout.write("There are now %s read files in the directory - this is the number of comparisons for Bowtie2/Minimap2\n\n" % (len(os.listdir(output_dir+'/reads_mapped'))))
  sys.stdout.flush()
else:
  sys.stdout.write("Skipping 5_combined_files because check-point wasn't 4_extracted_reads\n\n")
  sys.stdout.flush()
  with open(output_dir+'/pickle_intermediates/all_files.pickle', 'rb') as f:
    all_files = pickle.load(f)

# 6. Make bowtie2 databases & run bowtie2 for all taxonomy ID's
if cp == "5_combined_files":
  if coverage_program in ['Bowtie2', 'Both']:
    sys.stdout.write("Running Bowtie2\n")
    sys.stdout.flush()
    make_bowtie2_databases(taxid, output_dir, n_proc, bowtie2_db_dir, genome_dir)
    run_bowtie2_paf(all_files, taxid, output_dir, n_proc, bowtie2_setting, bowtie2_db_dir)
    cp = update_checkpoint(output_dir, "6_bowtie2_run")
    sys.stdout.write("Completed check-point 6 Bowtie2 run\n\n")
    sys.stdout.flush()
  else:
    sys.stdout.write("Skipping 6_bowtie2_run because Bowtie2 wasn't in coverage_program\n\n")
    sys.stdout.flush()
    cp = update_checkpoint(output_dir, "6_bowtie2_run")
else:
  sys.stdout.write("Skipping 6_bowtie2_run because check-point wasn't 5_combined_files\n\n")
  sys.stdout.flush()

# 7. Run Minimap2
if cp == "6_bowtie2_run":
  if coverage_program in ['Minimap2', 'Both']:
    sys.stdout.write("Running Minimap2\n")
    sys.stdout.flush()
    run_minimap2(all_files, taxid, output_dir, n_proc, genome_dir)
    cp = update_checkpoint(output_dir, "7_minimap2_run")
    sys.stdout.write("Completed check-point 7 Minimap2 run\n\n")
    sys.stdout.flush()
  else:
    sys.stdout.write("Skipping 7_minimap2_run because Minimap2 wasn't in coverage_program\n\n")
    sys.stdout.flush()
    cp = update_checkpoint(output_dir, "7_minimap2_run")
else:
  sys.stdout.write("Skipping 7_minimap2_run because check-point wasn't 6_bowtie2_run\n\n")
  sys.stdout.flush()

# 8. Get the coverage and mapping of reads across the genomes
if cp == "7_minimap2_run":
  if skip_coverage:
    cp = update_checkpoint(output_dir, "8_got_coverage")
    sys.stdout.write("Skipped check-point 8 getting coverage across genomes\n\n")
    sys.stdout.flush()
  else:
    sys.stdout.write("Getting coverage across genomes\n")
    sys.stdout.flush()
    get_genome_info(taxid, output_dir, genome_dir, n_proc)
    get_coverage_across_genomes_paf(all_files, taxid, genome_dir, output_dir, n_proc, coverage_program, mapq_threshold, identity_threshold)
    cp = update_checkpoint(output_dir, "8_got_coverage")
    sys.stdout.write("Completed check-point 8 got coverage across genomes\n\n")
    sys.stdout.flush()
else:
  sys.stdout.write("Skipping 8_got_coverage because check-point wasn't 7_minimap2_run\n\n")
  sys.stdout.flush()

# 9. Collate the output
if cp == "8_got_coverage":
  sys.stdout.write("Collating all output\n")
  sys.stdout.flush()
  collate_output_paf(all_files, taxid, output_dir, kreports, samples, group_samples, skip_coverage, coverage_program, genome_dir, grouped_samples_only, no_grouped_samples)
  cp = update_checkpoint(output_dir, "9_collate_output")
  sys.stdout.write("Completed check-point 9 collating output\n\n")
  sys.stdout.flush()
else:
  sys.stdout.write("Skipping 9_collate_output because check-point wasn't 8_got_coverage\n\n")
  sys.stdout.flush()

# 10. Clean up all of the intermediate files
if cp == "9_collate_output":
  if skip_cleanup:
    cp = update_checkpoint(output_dir, "10_cleanup")
    sys.stdout.write("Skipping check-point 10 cleanup of intermediate files\n\n")
    sys.stdout.flush()
  else:
    print("Cleaning up directory\n")
    clean_up(output_dir)
    cp = update_checkpoint(output_dir, "10_cleanup")
    sys.stdout.write("Completed check-point 10 cleanup of intermediate files\n\n")
    sys.stdout.flush()
else:
  sys.stdout.write("Skipping 10_cleanup because check-point wasn't 9_collate_output\n\n")
  sys.stdout.flush()

sys.stdout.write("Finished running genome coverage checker pipeline.\n")
sys.stdout.write("Running time: --- %s seconds ---\n\n" % str(round((time.time() - start_time), 2)))
sys.stdout.flush()

f.close()
sys.stdout = o

print('See '+'Genome Coverage Checker log '+time_string+'.txt for a complete log of the Genome Coverage Checker run!\n')
