#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import os
import sys
import math

parser = argparse.ArgumentParser(description='This script is to plot the consistency of taxa after the coverage_pipeline.py script has been run.')
parser.add_argument('--taxid', dest='taxid', default=None,
                    help="Which taxonomy ID's to plot. This can be a single taxonomy ID or a list separated by commas (e.g. 980563,622,2949971)")
parser.add_argument('--top_taxa', dest='top_taxa', default=10,
                    help="How many taxa to plot. Default is the top 10. This can be any number 1-number of taxa in your samples, and if you say all then it will plot all of them")
parser.add_argument('--coverage_program', dest='coverage_program', default='Bowtie2', choices=['Minimap2', 'Bowtie2', 'Both'],
                    help="Which of the programs to use for plotting coverage across the genome. Default is Bowtie2.")
parser.add_argument('--project_folder', dest='project_folder', default=None,
                    help="The folder containing the coverage checker output. It is expected that this contains coverage_checker_output.tsv and the coverage folder at a minimum.")
parser.add_argument('--dpi', dest='dpi', default=300,
                    help="The dpi to save the figures with. Note that this may need to be reduced if plotting a large number of taxa or samples.")
parser.add_argument('--granularity', dest='granularity', default=1,
                    help="Used as the interval for plotting coverage. Default is 1. Increasing this will speed up plotting time but decrease accuracy.")
parser.add_argument('--sample_group', dest='sample_group', default=None,
                    help="Which sample group to plot.")
parser.add_argument('--read_limit', dest='read_limit', default=1,
                    help="The number of reads mapped to a position within a sample to consider it as present in that sample. Default is 1.")
# parser.add_argument('--sample_metadata', dest='sample_metadata', default=None,
#                     help="Location of the sample metadata file. It is expected that this will be a CSV (comma separated) file with a header and two columns. The first column contains sample names and the second contains the sample groupings")
parser.add_argument('--no_plots', dest='no_plots', default=False, action='store_true',
                    help="Don't make any genome consistency plots. By default these will be made for all taxid's run.")



args = parser.parse_args()
taxid, top_taxa, coverage_program, project_folder, dpi, granularity, sample_group, limit, no_plots = args.taxid, args.top_taxa, args.coverage_program, args.project_folder, args.dpi, args.granularity, args.sample_group, args.read_limit, args.no_plots
project_folder = project_folder+'/'

if not os.path.exists(project_folder+'figures/'):
  os.system('mkdir '+project_folder+'figures')

if not os.path.exists(project_folder+'consistency/'):
  os.system('mkdir '+project_folder+'consistency')
  
if top_taxa not in [None, 'all']: top_taxa = int(top_taxa)
dpi = int(dpi)
granularity = int(granularity)

if not no_plots:
  plotting = True
else:
  plotting = False

with open(project_folder+'/pickle_intermediates/group_samples.pickle', 'rb') as f:
    group_samples = pickle.load(f)
    
if sample_group not in group_samples:
  all_groups = [group for group in group_samples]
  sys.exit("The group that you have given for --sample_group is not in what you run GeCoCheck coverage_pipeline.py with. The options you have are: "+','.join(all_groups))

samples = group_samples[sample_group]
cc_out = pd.read_csv(project_folder+'/coverage_checker_output.tsv', index_col=0, header=0, sep='\t')
cc_out['taxid'] = cc_out['taxid'].map(str)
cc_out_tax = cc_out.loc[sample_group, :].set_index('taxid')
columns = ['Samples with reads > limit', 'Mean reads mapped to position', 'Median reads mapped to position', 'Positions with reads mapped > limit']
for col in columns:
  cc_out_tax[col] = ''

if taxid == None:
  if top_taxa == None:
    sys.exit('You must set one of taxid and top_taxa to run')
  taxid_list = cc_out.loc[sample_group, :]
  taxid_list = taxid_list[taxid_list[coverage_program+' reads mapped'] > 0]
  if top_taxa == 'all':
    taxid_list = list(taxid_list['taxid'].values)
    sys.stdout.write("You chose to run this with all taxa. That is "+str(len(taxid_list))+" taxa.\n")
    sys.stdout.flush()
  else:
    taxid_list = taxid_list.sort_values(by=['Kraken reads assigned'], ascending=False)
    taxid_list = taxid_list.head(int(top_taxa))
    taxid_list = list(taxid_list['taxid'].values)
elif ',' in taxid:
  taxid_list = taxid.split(',')
else:
  taxid_list = [taxid]

cc_out_tax = cc_out_tax.loc[taxid_list, :]

def get_coverage_single_sample(sample_name, taxid, length, program, all_genomes, limit):
  fn = project_folder+'coverage/'+program.replace('B', 'b').replace('M', 'm')+'_'+sample_name+'_'+taxid+'.txt'
  genome_covered = {}
  for l in range(int(length)):
    genome_covered[l] = 0
  if not os.path.exists(fn):
    print("File doesn't exist: "+fn)
    return genome_covered, all_genomes
  starts = ['']
  for row in open(fn, 'r'):
    if 'all_starting_points: ' in row:
      starts = row.split('all_starting_points: ')[1].replace('\n', '').split(',')
    elif 'all_end_points: ' in row:
      ends = row.split('all_end_points: ')[1].replace('\n', '').split(',')
    elif 'genome_identity: ' in row:
      ids = row.split('genome_identity: ')[1].replace('\n', '').split(',')
  if starts == ['']:
    return genome_covered, all_genomes
  starts = [int(n) for n in starts]
  ends = [int(n) for n in ends]
  ids = [float(n) for n in ids]
  for se in range(len(starts)):
    start, end = starts[se], ends[se]
    for g in range(start, end+1):
      try:
        genome_covered[g] += 1
      except:
        genome_covered[g] = 1
  for b in genome_covered:
    if genome_covered[b] >= limit:
      try:
        all_genomes[b] += 1
      except:
        all_genomes[b] = 1
  return genome_covered, all_genomes

def plot_genome_consistency(this_genome, name, mean, median, genome_frac, mapped_frac, save_name):
  
  fig = plt.figure(figsize=(20,3))
  ax_genome = plt.subplot2grid((1,25),(0,0), colspan=10)
  ax_cov_frac = plt.subplot2grid((1,25),(0,11), colspan=2)
  ax_map_frac = plt.subplot2grid((1,25),(0,13), colspan=2)
  ax_mean = plt.subplot2grid((1,25),(0,15), colspan=2)
  ax_median = plt.subplot2grid((1,25),(0,17), colspan=2)
  
  k, vals = list(all_genomes.keys()), list(all_genomes.values())
  groups = [vals[x:x+granularity] for x in range(0, len(vals), granularity)]
  groups_gen = [k[x:x+granularity] for x in range(0, len(k), granularity)]
  means = [sum(group)/len(group) for group in groups]
  gen_means = [sum(group)/len(group) for group in groups_gen]
  means_df = pd.DataFrame(means, index=gen_means).transpose()
  # means = list(all_genomes.values())
  # gen_means = list(all_genomes.keys())
  means_df = pd.DataFrame(means, index=gen_means).transpose()
  plt.sca(ax_genome)
  ax_genome.pcolor(means_df, vmin=0, vmax=1, cmap='plasma')
  xt = list(plt.xticks()[0][:-1])
  t = plt.xticks(xt, [round(gen_means[int(x)]/1000000, 1) for x in xt])
  yt = plt.yticks([])
  xl = plt.xlabel('Genome size (mbp)')
  ti = plt.title(name+'\n\n'+coverage_program+'\nCoverage across genome', fontweight='bold')
  
  vals, axes = [genome_frac, mapped_frac, mean, median], [ax_cov_frac, ax_map_frac, ax_mean, ax_median]
  cmaps, limits = ['RdPu', 'OrRd', 'GnBu', 'BuPu'], [[0, 100], [0, 100], [0, 1], [0, 1]]
  titles = [coverage_program+' genome\nfraction (%)', 'Genome fraction (%)\n> '+str(limit)+' reads', 'Mean reads mapped\nper genome position', 'Median reads mapped\nper genome position']
  for v in range(len(vals)):
    plt.sca(axes[v])
    m = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=limits[v][0], vmax=limits[v][1]), cmap=cmaps[v])
    ba = plt.bar([0], [1], width=1, color=m.to_rgba(vals[v]), edgecolor='k')
    xl = plt.xlim([-0.5, 0.5]), plt.ylim([0, 1]), plt.xticks([]), plt.yticks([])
    if vals[v] <= np.mean([limits[v][0], limits[v][1]]): fc = 'k'
    else: fc = 'w'
    tx = plt.text(0, 0.5, str(round(vals[v], 2)), ha='center', va='center', color=fc)
    ti = plt.title(titles[v], fontweight='bold', rotation=90)
    
  #fig.suptitle(name, fontweight='bold')
  plt.savefig(project_folder+'figures/'+save_name, dpi=dpi, bbox_inches='tight')
  
  return means_df

for taxid in taxid_list:
  print('Working on taxonomy ID: ', taxid)
  length = cc_out_tax.loc[taxid, 'Reference genome length (bp)']
  all_genomes = {}
  in_sample_count = 0
  for l in range(int(length)):
    all_genomes[l] = 0
  for sample in samples:
    # print(sample)
    sample_cc = cc_out.loc[sample, :].set_index('taxid')
    if taxid in sample_cc.index.values:
      if sample_cc.loc[taxid, coverage_program+' reads mapped'] > limit:
        in_sample_count += 1
        this_sample, all_genomes = get_coverage_single_sample(sample, taxid, length, coverage_program, all_genomes, limit)
  non_zero = [all_genomes[b]/in_sample_count for b in all_genomes if all_genomes[b] > 0]
  try:
    mean_val, median_val = np.mean(non_zero), np.median(non_zero)
  except:
    mean_val, median_val = 0, 0
  nums = [in_sample_count, mean_val, median_val, len(non_zero)]
  print('In samples: '+str(in_sample_count))
  print('Mean reads mapped: '+str(mean_val))
  print('Median reads mapped: '+str(median_val))
  print(coverage_program+' genome fraction: ', cc_out_tax.loc[taxid, coverage_program+' genome fraction (%)'], '\nProportion of genome with at least '+str(limit)+' mapped reads: ', len(non_zero)/len(all_genomes))
  for c in range(len(nums)):
    cc_out_tax.loc[taxid, columns[c]] = str(nums[c])
  cc_out_tax.loc[[taxid], :].to_csv(project_folder+'consistency/coverage_checker_output_consistency_tax'+taxid+'_limit'+str(limit)+'.tsv', sep='\t')
  name = taxid+': $'+cc_out_tax.loc[taxid, 'Species name'].split('_', 1)[0]+'$ $'+cc_out_tax.loc[taxid, 'Species name'].split('_', 1)[1]+'$\nSamples: $n$='+str(in_sample_count)+'/'+str(len(samples))
  save_name = taxid+'_'+sample_group+'_'+cc_out_tax.loc[taxid, 'Species name']+'_consistency.png'
  # plot_genome_consistency(all_genomes, name, mean_val, median_val, cc_out_tax.loc[taxid, coverage_program+' genome fraction (%)'], (len(non_zero)/len(all_genomes))*100, len(samples), save_name)
  if plotting:
    if not os.path.exists(save_name):
      plot_genome_consistency(all_genomes, name, mean_val, median_val, cc_out_tax.loc[taxid, coverage_program+' genome fraction (%)'], (len(non_zero)/len(all_genomes))*100, save_name)

cc_out_tax.to_csv(project_folder+'coverage_checker_output_consistency_'+str(limit)+'.tsv', sep='\t')


