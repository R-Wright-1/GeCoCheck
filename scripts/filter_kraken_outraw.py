#!/opt/anaconda3/envs/GeCoCheck-v0.0.4/bin/python3.13

import os
import pandas as pd
import sys
import argparse
import pickle
from genome_coverage_checker.util import *
import time

start_time = time.time()

o = sys.stdout
time_string = str(time.ctime(start_time)).replace(':', '-')
print('Starting at '+time_string+'.txt\n')

parser = argparse.ArgumentParser(description='This script is to check which taxa reads have been assigned to by Kraken, pull out these reads, download reference genomes for the taxa, and map the reads to the reference genomes.')
parser.add_argument('--kraken_outraw_dir', dest='kraken_outraw_dir', default=None,
                    help="The directory containing the kraken outraw files")
parser.add_argument('--gecocheck_out_dir', dest='gecocheck_out_dir', default=None,
                    help="The directory to save the output files to")
parser.add_argument('--kraken_kreport_dir', dest='kraken_kreport_dir', default=None,
                    help="The directory containing the kraken kreport files")
parser.add_argument('--sample_metadata', dest='sample_metadata', default=None,
                    help="Location of the sample metadata file. It is expected that this will be a CSV (comma separated) file with a header and two columns. The first column contains sample names and the second contains the sample groupings")
parser.add_argument('--koutraw_filt_dir', dest='koutraw_filt_dir', default=None,
                    help="The new directory for the filtered kraken outraw files to be saved to")
parser.add_argument('--project_name', dest='project_name', default=None,
                    help="Name of this project")
parser.add_argument('--samples_grouping', dest='samples_grouping', default='overall', choices=['overall', 'per_sample', 'per_grouping'],
                    help="How to filter out the taxids that will be removed. Default is in the entire project, but other options are per sample [per_sample] or per sample grouping [per_grouping]")
parser.add_argument('--reads_or_genome_frac', dest='reads_or_genome_frac', default='reads', choices=['reads', 'genome_frac'],
                    help="Whether to filter taxids based on number of reads assigned or reference genome fraction (%) present")
parser.add_argument('--read_limit', dest='read_limit', default=None,
                    help="Keep taxids with >= this number of reads mapped to the reference genome. Default will take from what you ran GeCoCheck with (by not setting this at all).")
parser.add_argument('--genome_frac_limit', dest='genome_frac_limit', default=None,
                    help="Keep taxids with >= this reference genome % present. Can take any value between 0-100. By default this is unused and the read_limit is used instead.")
parser.add_argument('--coverage_program', dest='coverage_program', default='Bowtie2', choices=['Minimap2', 'Bowtie2'],
                    help="Which of the programs was used for getting coverage across the genome. Default is Bowtie2.")
parser.add_argument('--unchecked', dest='unchecked', default='unclassified', choices=['classified', 'unclassified'],
                    help="Whether taxids that were not checked by GeCoCheck should be classified or unclassified in the output. By default they are unclassified")
parser.add_argument('--for_stratified_out', dest='for_stratified_out', default=False, action='store_true',
                    help="Whether to write a second raw output file that doesn't have unverified taxa as unclassified but indicated that they are unverified")


args = parser.parse_args()
kraken_outraw_dir, gecocheck_out_dir, kraken_kreport_dir, sample_metadata, koutraw_filt_dir, project_name, samples_grouping, reads_or_genome_frac, coverage_program, unchecked, read_limit, genome_frac_limit, for_stratified_out = args.kraken_outraw_dir, args.gecocheck_out_dir, args.kraken_kreport_dir, args.sample_metadata, args.koutraw_filt_dir, args.project_name, args.samples_grouping, args.reads_or_genome_frac, args.coverage_program, args.unchecked, args.read_limit, args.genome_frac_limit, args.for_stratified_out

if not kraken_outraw_dir:
  sys.exit('You need to give the --kraken_outraw_dir')
if not gecocheck_out_dir:
  sys.exit('You need to give the --gecocheck_out_dir')
if not kraken_kreport_dir:
  sys.exit('You need to give the --kraken_kreport_dir')
if not sample_metadata:
  sys.exit('You need to give the --sample_metadata file name')
if not koutraw_filt_dir:
  print("You didn't give the --koutraw_filt_dir. Setting this to "+gecocheck_out_dir+"/kraken2_outraw_checked/")
  koutraw_filt_dir = gecocheck_out_dir+"/kraken2_outraw_checked/"
if not project_name:
  sys.exit('You need to give the --project_name. This should be the same as what you used for running GeCoCheck')

if reads_or_genome_frac == 'reads' and not read_limit:
  print("You have chosen to filter on the number of reads mapped to the reference genome and haven't provided a read limit. This will be taken from what you used for running GeCoCheck.")
else:
  read_limit = int(read_limit)
if reads_or_genome_frac == 'genome_frac' and not genome_frac_limit:
  sys.exit('"You have chosen to filter on the reference genome fraction (%) present so you must give a value for --genome_frac_limit. This can be any number between 0-100')
elif reads_or_genome_frac == 'genome_frac':
  genome_frac_limit = float(genome_frac_limit)
  
metadata = pd.read_csv(sample_metadata, index_col=0, header=0)
gecocheck_out = pd.read_csv(gecocheck_out_dir+'/coverage_checker_output.tsv', index_col=0, header=0, sep='\t')

if not os.path.exists(koutraw_filt_dir):
  mk = os.system('mkdir '+koutraw_filt_dir)
  print('Made a directory for the filtered kraken outraw files: '+koutraw_filt_dir)
else:
  files = os.listdir(koutraw_filt_dir)
  if len(files) != 0:
    sys.exit(koutraw_filt_dir+"isn't empty! Please try again with an empty directory")

sample_groups, groups_samples_in = {}, {}
for r in range(len(metadata.index.values)):
  group = metadata.iloc[r, 0]
  groups_samples_in[metadata.index.values[r]] = group
  if group in sample_groups:
    sample_groups[group].append(metadata.index.values[r])
  else:
    sample_groups[group] = [metadata.index.values[r]]

with open(gecocheck_out_dir+'/pickle_intermediates/args.pickle', 'rb') as f:
  cov_args = pickle.load(f)

read_lim, read_mean = cov_args[7], cov_args[8]

if not read_limit:
  read_limit = read_lim

with open(gecocheck_out_dir+'/pickle_intermediates/taxid.pickle', 'rb') as f:
    taxid_checking = pickle.load(f)

taxid_checking = set(taxid_checking)

kreports = pd.read_csv(gecocheck_out_dir+project_name+'_combined_kreport.csv', index_col=0, header=0)
kreports.index = kreports.index.map(str)
above_limit = []
for group in sample_groups:
  krep_group = kreports.copy(deep=True).loc[:, sample_groups[group]]
  krep_group = krep_group[krep_group.max(axis=1) >= read_lim]
  if read_mean != None:
    krep_group['Mean'] = krep_group.mean(axis=1)
    krep_group = krep_group[krep_group['Mean'] >= read_mean]
    krep_group = krep_group.drop('Mean', axis=1)
  above_limit += list(krep_group.index.values)

above_limit = set(above_limit)
below_limit = [tid for tid in set(list(kreports.index.values)) if tid not in above_limit]
below_limit = set(below_limit)
above_limit_not_checked = [tid for tid in above_limit if tid not in taxid_checking]

print('Found '+str(len(above_limit))+' taxids to check that were above the read limits and '+str(len(below_limit))+' taxids that were below the read limits. Of those that were above the read limit '+str(len(above_limit_not_checked))+' were not checked because they were either in the wrong domain or genomes were not found')

kreports['GeCoCheck'] = ''
kreports.loc[list(above_limit), 'GeCoCheck'] = 'Coverage checked'
kreports.loc[list(below_limit), 'GeCoCheck'] = 'Not checked: below read limit'
kreports.loc[list(above_limit_not_checked), 'GeCoCheck'] = 'Not checked: wrong domain or not found'

kreports.to_csv(gecocheck_out_dir+project_name+'_combined_kreport_filtering.csv')
  
children, parents = {}, {}
for sample_name in metadata.index.values:
  krep = kraken_kreport_dir+sample_name+'.kreport'
  if not os.path.exists(krep): krep = kraken_kreport_dir+sample_name+'_0.0.kreport'
  for row in open(krep, 'r'):
    this_row = row.split('\t')
    if this_row[3] == 'S':
      sp = this_row[4]
    elif this_row[3] in ['S1', 'S2', 'S3', 'S4']:
      parents[this_row[4]] = sp
      if sp in children:
        children[sp].append(this_row[4])
      else:
        children[sp] = [this_row[4]]

samples_keeping, samples_keeping_species = {}, {}
for sample in metadata.index.values:
  if samples_grouping == 'overall':
    this_gecocheck_out = gecocheck_out.loc[project_name, :].set_index('taxid')
  elif samples_grouping == 'per_grouping':
    this_gecocheck_out = gecocheck_out.loc[groups_samples_in[sample], :].set_index('taxid')
  else:
    this_gecocheck_out = gecocheck_out.loc[sample, :].set_index('taxid')
  if reads_or_genome_frac == 'reads':
    this_gecocheck_out = this_gecocheck_out[this_gecocheck_out[coverage_program+' reads mapped'] >= read_limit]
  elif reads_or_genome_frac == 'genome_frac':
    this_gecocheck_out = this_gecocheck_out[this_gecocheck_out[coverage_program+' genome fraction (%)'] >= genome_frac_limit]
  this_gecocheck_out.index = this_gecocheck_out.index.map(str)
  verified = set(list(this_gecocheck_out.index.values))
  #print(sample+': '+str(len(verified))+' taxids verified by GeCoCheck')
  if unchecked == 'classified':
    verified.update(above_limit_not_checked)
  samples_keeping_species[sample] = set(verified)
  strains = set()
  for tid in verified:
    if tid in children:
      strains.update(children[tid])
  verified.update(strains)
  samples_keeping[sample] = verified
  #print(sample+': '+str(len(verified))+' taxids keeping classified')
  
#check all taxids string
for sample in samples_keeping:
  for tid in samples_keeping[sample]:
    if not isinstance(tid, str):
      print(tid)

#process samples
for sample in metadata.index.values:
  in_fn = kraken_outraw_dir+sample+'.kraken'
  out_fn = koutraw_filt_dir+sample+'.kraken'
  count, changed, unclassified, classified_initially = 0, 0, 0, 0
  with open(in_fn, 'r') as infile, open(out_fn, 'w') as outfile:
    for line in infile:
      if line[0] == 'U':
        quiet_write = outfile.write(line)
        unclassified += 1
      else:
        classified_initially += 1
        working_line = line.split('\t')
        tid = working_line[2].split('(taxid ')[1].replace(')', '')
        if tid in samples_keeping[sample]:
          quiet_write = outfile.write(line)
        else:
          working_line[0] = 'U'
          working_line[2] = 'unclassified (taxid 0)'
          working_line = '\t'.join(working_line)
          quiet_write = outfile.write(working_line)
          changed += 1
      count += 1
  print(sample+': '+str(count)+' reads in file. '+str(classified_initially)+' ('+str(round((classified_initially/count)*100, 2))+'%) were classified and '+str(unclassified)+' ('+str(round((unclassified/count)*100, 2))+'%) were unclassified. An additional '+str(changed)+' ('+str(round((changed/count)*100, 2))+'%) are now unclassified.')
  if for_stratified_out:
    print('Making outraw file that can be used with workflow to generate stratified taxonomy/function output')
    out_fn = koutraw_filt_dir+sample+'_unverified.kraken'
    with open(in_fn, 'r') as infile, open(out_fn, 'w') as outfile:
      for line in infile:
        if line[0] == 'U':
          quiet_write = outfile.write(line)
        else:
          working_line = line.split('\t')
          tid = working_line[2].split('(taxid ')[1].replace(')', '')
          if tid in samples_keeping[sample]:
            quiet_write = outfile.write(line)
          else:
            working_line[2] = 'Unverified '+working_line[2]
            working_line = '\t'.join(working_line)
            quiet_write = outfile.write(working_line)

sys.stdout.write("Finished filtering reads in samples.\n")
sys.stdout.write("Running time: --- %s seconds ---\n\n" % str(round((time.time() - start_time), 2)))
