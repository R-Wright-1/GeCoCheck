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

parser = argparse.ArgumentParser(description='This script is to plot the coverage of taxa after the coverage_pipeline.py script has been run.')
parser.add_argument('--running', dest='running', default='taxon', choices=['taxon','sample'],
                    help='Whether to run in sample oriented or taxon oriented format.')
parser.add_argument('--taxid', dest='taxid', default=None,
                    help="Which taxonomy ID's to plot. This can be a single taxonomy ID or a list separated by commas (e.g. 980563,622,2949971)")
parser.add_argument('--top_taxa', dest='top_taxa', default=30,
                    help="How many taxa to plot (only used for sample oriented format)")
parser.add_argument('--coverage_program', dest='coverage_program', default='Bowtie2', choices=['Minimap2', 'Bowtie2', 'Both'],
                    help="Which of the programs to use for plotting coverage across the genome. Default is Bowtie2.")
parser.add_argument('--sort_by', dest='sort_by', default='kraken', choices=['genome_fraction', 'kraken', 'minimap2', 'bowtie2'],
                    help="How to determine which are the top taxa. Note that if you choose minimap2/bowtie2 and ran coverage checker without a read limit then this may not be very helpful.")
parser.add_argument('--project_folder', dest='project_folder', default=None,
                    help="The folder containing the coverage checker output. It is expected that this contains coverage_checker_output.tsv and the coverage folder at a minimum.")
parser.add_argument('--dpi', dest='dpi', default=300,
                    help="The dpi to save the figures with. Note that this may need to be reduced if plotting a large number of taxa or samples.")
parser.add_argument('--granularity', dest='granularity', default=1000,
                    help="Used as the interval for plotting coverage. Default is 1000. Decreasing this will slow down plotting time but increase accuracy, and increasing it will speed it up and decrease accuracy.")
parser.add_argument('--samples', dest='samples', default=None,
                    help="Which sample(s) to plot. This can be a single sample name or a list separated by commas (e.g. Sample1,Sample2,Sample3)")

args = parser.parse_args()
running, taxid, top_taxa, sort_by, project_folder, dpi, granularity, samples, coverage_program = args.running, args.taxid, args.top_taxa, args.sort_by, args.project_folder, args.dpi, args.granularity, args.samples, args.coverage_program
project_folder = project_folder+'/'

if not os.path.exists(project_folder+'figures/'):
  os.system('mkdir '+project_folder+'figures')
  
if top_taxa != None: top_taxa = int(top_taxa)
dpi = int(dpi)
granularity = int(granularity)

if taxid == None:
  if top_taxa == None:
    sys.exit('You must set one of taxid and top_taxa to run')
  taxid_list = None
elif ',' in taxid:
  taxid_list = taxid.split(',')
else:
  taxid_list = [taxid]

if samples == 'All' or samples == 'all':
  samples = 'All'
elif samples == None:
  samples = 'All'
elif ',' in samples:
  samples = samples.split(',')
else:
  samples = [samples]

# make plot for single taxon/genome
def plot_genome_coverage(axes_genome, axes_id, sample_name, taxid, length, program):
  fn = project_folder+'coverage/'+program.replace('B', 'b').replace('M', 'm')+'_'+sample_name+'_'+taxid+'.txt'
  if not os.path.exists(fn):
    for ax in [axes_genome, axes_id]:
      plt.sca(ax)
      xl = plt.xlim([68, 102]), plt.ylim(-0.5, 0.5)
      if ax == axes_genome: text = 'No reads mapped'
      else: text = 'NA'
      tx = plt.text(0.5, 0.5, text, ha='center', va='center', transform=ax.transAxes)
      xt = plt.xticks([]), plt.yticks([])
    return 'NA', 'NA'
  starts = ['']
  for row in open(fn, 'r'):
    if 'all_starting_points: ' in row:
      starts = row.split('all_starting_points: ')[1].replace('\n', '').split(',')
    elif 'all_end_points: ' in row:
      ends = row.split('all_end_points: ')[1].replace('\n', '').split(',')
    elif 'genome_identity: ' in row:
      ids = row.split('genome_identity: ')[1].replace('\n', '').split(',')
  if starts == ['']:
    for ax in [axes_genome, axes_id]:
      plt.sca(ax)
      xl = plt.xlim([68, 102]), plt.ylim(-0.5, 0.5)
      if ax == axes_genome: text = 'No reads mapped'
      else: text = 'NA'
      tx = plt.text(0.5, 0.5, text, ha='center', va='center', transform=ax.transAxes)
      xt = plt.xticks([]), plt.yticks([])
    return 'NA', 'NA'
  starts = [int(n) for n in starts]
  ends = [int(n) for n in ends]
  ids = [float(n) for n in ids]
  plt.sca(axes_id)
  sc = plt.scatter(ids, np.random.normal(0, 0.12, len(ids)), color='#F4D03F', alpha=0.01)
  box = plt.boxplot(ids, positions=[0], widths=0.8, vert=False, showfliers=False)
  for item in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']: plt.setp(box[item], color='k')
  tx = plt.text(np.median(ids), 0.65, str(round(np.median(ids), 3)), ha='center', va='center')
  xl = plt.xlim([68, 102]), plt.xticks([]), plt.yticks([]), plt.ylim([-0.5, 0.9])
  genome_covered = {}
  for l in range(int(length)):
    genome_covered[l] = 0
  for se in range(len(starts)):
    start, end = starts[se], ends[se]
    for g in range(start, end+1):
      try:
        genome_covered[g] += 1
      except:
        genome_covered[g] = 1
  k, vals = list(genome_covered.keys()), list(genome_covered.values())
  groups = [vals[x:x+granularity] for x in range(0, len(vals), granularity)]
  groups_gen = [k[x:x+granularity] for x in range(0, len(k), granularity)]
  means = [sum(group)/len(group) for group in groups]
  gen_means = [sum(group)/len(group) for group in groups_gen]
  means_df = pd.DataFrame(means, index=gen_means).transpose()
  plt.sca(axes_genome)
  axes_genome.pcolor(means_df, vmin=0, vmax=1)
  xt = list(plt.xticks()[0][:-1])
  #t = plt.xticks(xt, [round(gen_means[int(x)]/1000000, 1) for x in xt])
  yt = plt.yticks([])
  return xt, gen_means
    

def plot_square(axes, cmap, number, min_val, max_val, rnd=None, div=None):
  plt.sca(axes)
  if not np.isnan(number):
    ba = plt.bar([0], [1], width=1, color=cmap.to_rgba(number), edgecolor='k')
    xl = plt.xlim([-0.5, 0.5]), plt.ylim([0, 1]), plt.xticks([]), plt.yticks([])
    if number <= np.mean([min_val, max_val]): fc = 'k'
    else: fc = 'w'
    if div == 'mil' and max_val > 100000:
      num = number/10000
      num = str(round(num, 2))
    else:
      if rnd == None: num = str(round(number))
      else: num = str(round(number, rnd))
    tx = plt.text(0, 0.5, num, ha='center', va='center', color=fc)
  else:
    ba = plt.bar([0], [1], width=1, color='w', edgecolor='k')
    xl = plt.xlim([-0.5, 0.5]), plt.ylim([0, 1]), plt.xticks([]), plt.yticks([])
    tx = plt.text(0, 0.5, 'NA', ha='center', va='center', color='k')
  return

# make plot for single taxon across all samples
def single_taxon_across_samples(taxid, save_name, project_folder, coverage_program, samples='All'):

  cc_out = pd.read_csv(project_folder+'coverage_checker_output.tsv', index_col=1, header=0, sep='\t')
  cc_out = cc_out.loc[int(taxid), :]
  species_name = cc_out['Species name'].values[0].replace('_', ' ')
  if coverage_program == 'Both':
    cc_out = cc_out.loc[:, ['Sample', 'Reference genome length (bp)', 'Kraken reads assigned', 'Minimap2 reads mapped', 'Minimap2 genome fraction (%)', 'Proportion kraken reads mapped with Minimap2', 'Bowtie2 reads mapped', 'Bowtie2 genome fraction (%)', 'Proportion kraken reads mapped with Bowtie2']].set_index('Sample')
    cc_out['Proportion kraken reads mapped with Minimap2'] = cc_out['Proportion kraken reads mapped with Minimap2']*100
    cc_out['Proportion kraken reads mapped with Bowtie2'] = cc_out['Proportion kraken reads mapped with Bowtie2']*100
  elif coverage_program == 'Minimap2':
    cc_out = cc_out.loc[:, ['Sample', 'Reference genome length (bp)', 'Kraken reads assigned', 'Minimap2 reads mapped', 'Minimap2 genome fraction (%)', 'Proportion kraken reads mapped with Minimap2']].set_index('Sample')
    cc_out['Proportion kraken reads mapped with Minimap2'] = cc_out['Proportion kraken reads mapped with Minimap2']*100
  elif coverage_program == 'Bowtie2':
    cc_out = cc_out.loc[:, ['Sample', 'Reference genome length (bp)', 'Kraken reads assigned', 'Bowtie2 reads mapped', 'Bowtie2 genome fraction (%)', 'Proportion kraken reads mapped with Bowtie2']].set_index('Sample')
    cc_out['Proportion kraken reads mapped with Bowtie2'] = cc_out['Proportion kraken reads mapped with Bowtie2']*100
  
  #get plotting order
  if samples == 'All':
    sample_order = sorted(list(cc_out.index.values))
  else:
    sample_order = samples
  
  l = len(sample_order)
  l += 3
  if coverage_program in ['Minimap2', 'Bowtie2']:
    fig = plt.figure(figsize=(20,l))
    programs = [coverage_program]
  elif coverage_program == 'Both':
    fig = plt.figure(figsize=(40,l))
    programs = ['Minimap2', 'Bowtie2']
  fig.suptitle(taxid+': $'+species_name.replace(' ', '$ $')+'$\n\n', fontweight='bold', fontsize=26, ha='right')
  for program in programs:
    colormaps, colnames = ['RdPu', 'GnBu', 'BuPu'], ['Kraken reads assigned', program+' genome fraction (%)', 'Proportion kraken reads mapped with '+program]
    plot_names = ['Kraken reads\nassigned', program+' genome\nfraction(%)', 'Kraken reads\nmapped by\n'+program+' (%)']
    mapping, mins, maxs = [], [], []
    for c in range(len(colormaps)):
      if colnames[c] not in cc_out.columns: continue
      all_c = list(cc_out[colnames[c]].values)
      all_c = [ac for ac in all_c if not np.isnan(ac)]
      min_c, max_c = min(all_c), max(all_c)
      m = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=min_c, vmax=max_c), cmap=colormaps[c])
      ap = mapping.append(m), mins.append(min_c), maxs.append(max_c)
      # if c == 0 and max_c > 100000:
      #   plot_names[0] = 'Kraken reads\nassigned (x10,000)'

    xticks, gen_means = [], []
    for s in range(len(sample_order)):
      if coverage_program != 'Both':
        ax_genome = plt.subplot2grid((l,25),(s+1,0), colspan=10)
        ax_identity = plt.subplot2grid((l,25),(s+1,11), colspan=3)
        ax_kraken = plt.subplot2grid((l,25),(s+1,15), colspan=2)
        ax_program_gf = plt.subplot2grid((l,25),(s+1,17), colspan=2)
        ax_program_prop = plt.subplot2grid((l,25),(s+1,19), colspan=2)
      elif program == 'Minimap2':
        ax_genome = plt.subplot2grid((l,50),(s+1,0), colspan=10)
        ax_identity = plt.subplot2grid((l,50),(s+1,21), colspan=3)
        ax_kraken = plt.subplot2grid((l,50),(s+1,28), colspan=2)
        ax_program_gf = plt.subplot2grid((l,50),(s+1,31), colspan=2)
        ax_program_prop = plt.subplot2grid((l,50),(s+1,36), colspan=2)
      elif program == 'Bowtie2':
        ax_genome = plt.subplot2grid((l,50),(s+1,10), colspan=10)
        ax_identity = plt.subplot2grid((l,50),(s+1,24), colspan=3)
        ax_kraken = plt.subplot2grid((l,50),(s+1,28), colspan=2)
        ax_program_gf = plt.subplot2grid((l,50),(s+1,33), colspan=2)
        ax_program_prop = plt.subplot2grid((l,50),(s+1,38), colspan=2)
      axes = [ax_kraken, ax_program_gf, ax_program_prop]
      rounding, div = [None, 3, 3, 3], [None, None, None, None]
      for a in range(len(axes)):
        try:
          plot_square(axes[a], mapping[a], cc_out.loc[sample_order[s], colnames[a]], mins[a], maxs[a], rnd=rounding[a], div=div[a])
        except:
          plot_square(axes[a], mapping[a], 0, mins[a], maxs[a], rnd=rounding[a], div=div[a])
        if s == 0:
          axes[a].set_title(plot_names[a], fontweight='bold', rotation=90)
      if program == 'Minimap2': ax_genome.set_ylabel(sample_order[s], fontweight='bold', rotation=0, ha='right', va='center')
      elif coverage_program != 'Both' and program == 'Bowtie2': ax_genome.set_ylabel(sample_order[s], fontweight='bold', rotation=0, ha='right', va='center')
      try:
        xt, gm = plot_genome_coverage(ax_genome, ax_identity, sample_order[s], taxid, cc_out.loc[sample_order[s], 'Reference genome length (bp)'], program)
      except:
        xt, gm = plot_genome_coverage(ax_genome, ax_identity, sample_order[s], taxid, 0, program)
      if xt != 'NA':
        xticks, gen_means = xt, gm
      if s == 0:
        ax_genome.set_title(program+'\nCoverage across genome', fontweight='bold')
        ax_identity.set_title(program+'\nIdentity to\nreference', fontweight='bold')
      if s != len(sample_order)-1:
        plt.sca(ax_genome)
        xt = plt.xticks([])
      else:
        plt.sca(ax_identity)
        xt = plt.xticks([70, 80, 90, 100])
        xl = plt.xlabel('Identity (%)')
        plt.sca(ax_genome)
        t = plt.xticks(xticks, [round(gen_means[int(x)]/1000000, 1) for x in xticks])
        xl = plt.xlabel('Genome size (mbp)')

    ad = 3
    if coverage_program != 'Both':
      ax_genome_leg = plt.subplot2grid((l,25),((s)+ad,2),colspan=6)
      ax_kraken_leg = plt.subplot2grid((l,100),((s)+ad,61), colspan=6)
      ax_program_gf_leg = plt.subplot2grid((l,100),((s)+ad,69), colspan=6)
      ax_program_prop_leg = plt.subplot2grid((l,100),((s)+ad,77), colspan=6)
      axes = [ax_genome_leg, ax_kraken_leg, ax_program_gf_leg, ax_program_prop_leg]
    elif program == 'Minimap2':
      ax_genome_leg = plt.subplot2grid((l,50),((s)+ad,7),colspan=6)
      ax_kraken_leg = plt.subplot2grid((l,200),((s)+ad,112), colspan=8)
      ax_program_gf_leg = plt.subplot2grid((l,200),((s)+ad,128), colspan=8)
      ax_program_prop_leg = plt.subplot2grid((l,200),((s)+ad,148), colspan=8)
      axes = [ax_genome_leg, ax_kraken_leg, ax_program_gf_leg, ax_program_prop_leg]
    else:
      axes = []
    mapping = [mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap='viridis')]+mapping
    colormaps = ['viridis']+colormaps
    mins, maxs = [0]+mins, [1]+maxs
    labels = ['', 'Reads', 'Genome\nfraction (%)', 'Reads\nmapped (%)']
    for a in range(len(axes)):
      cb = mpl.colorbar.ColorbarBase(axes[a], cmap=colormaps[a], norm=mpl.colors.Normalize(vmin=mins[a], vmax=maxs[a]), orientation='horizontal')
      if a == 0:
        plt.sca(axes[a])
        plt.xticks([0.1, 0.9], ['Less coverage', 'More coverage'])
      else:
        plt.sca(axes[a])
        axes[a].set_xlabel(labels[a])
        if coverage_program == 'Both': plt.xticks([mins[a], maxs[a]], ['Low', 'High'])
        elif a == 1: plt.xticks([mins[a], maxs[a]])

  plt.subplots_adjust(hspace=0.7)
  plt.savefig(save_name+' '+species_name+' '+coverage_program+'.png', bbox_inches='tight', dpi=dpi)
  return

def multiple_taxa_in_one_sample(save_name, project_folder, taxid, sample, top_taxa, sort_by, coverage_program):
  program = 'Minimap2'
  cc_out = pd.read_csv(project_folder+'coverage_checker_output.tsv', index_col=0, header=0, sep='\t')
  cc_out = cc_out.loc[sample, :].set_index('taxid')
  
  if taxid != None:
    cc_out = cc_out.loc[[int(t) for t in taxid], :]
  else:
    if sort_by == 'genome_fraction': cc_out = cc_out.sort_values(by=['Minimap2 genome fraction (%)'], ascending=False)
    elif sort_by == 'kraken': cc_out = cc_out.sort_values(by=['Kraken reads assigned'], ascending=False)
    elif sort_by in ['minimap2', 'bowtie2']:
      cc_out = cc_out.sort_values(by=['Kraken reads assigned'], ascending=False)
      if sort_by == 'minimap2': cc_out = cc_out.sort_values(by=['Proportion kraken reads mapped with Minimap2'], ascending=False)
      elif sort_by == 'bowtie2': cc_out = cc_out.sort_values(by=['Proportion kraken reads mapped with Bowtie2'], ascending=False)
    save_name += '_'+sort_by+'_top'+str(top_taxa)
    cc_out = cc_out.head(top_taxa)
  if coverage_program == 'Both':
    cc_out = cc_out.loc[:, ['Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'Minimap2 reads mapped', 'Minimap2 genome fraction (%)', 'Proportion kraken reads mapped with Minimap2', 'Bowtie2 reads mapped', 'Bowtie2 genome fraction (%)', 'Proportion kraken reads mapped with Bowtie2']]
    cc_out['Proportion kraken reads mapped with Minimap2'] = cc_out['Proportion kraken reads mapped with Minimap2']*100
    cc_out['Proportion kraken reads mapped with Bowtie2'] = cc_out['Proportion kraken reads mapped with Bowtie2']*100
  elif coverage_program == 'Minimap2':
    cc_out = cc_out.loc[:, ['Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'Minimap2 reads mapped', 'Minimap2 genome fraction (%)', 'Proportion kraken reads mapped with Minimap2']]
    cc_out['Proportion kraken reads mapped with Minimap2'] = cc_out['Proportion kraken reads mapped with Minimap2']*100
  elif coverage_program == 'Bowtie2':
    cc_out = cc_out.loc[:, ['Species name', 'Reference genome length (bp)', 'Kraken reads assigned', 'Bowtie2 reads mapped', 'Bowtie2 genome fraction (%)', 'Proportion kraken reads mapped with Bowtie2']]
    cc_out['Proportion kraken reads mapped with Bowtie2'] = cc_out['Proportion kraken reads mapped with Bowtie2']*100
  
  plot_order = list(cc_out.index.values)
  l = len(plot_order)
  l += 3
  if coverage_program in ['Minimap2', 'Bowtie2']:
    fig = plt.figure(figsize=(20,l))
    programs = [coverage_program]
  elif coverage_program == 'Both':
    fig = plt.figure(figsize=(40,l))
    programs = ['Minimap2', 'Bowtie2']
  fig.suptitle(sample, fontweight='bold', fontsize=26, ha='right')
  for program in programs:
    colormaps, colnames = ['RdPu', 'GnBu', 'BuPu'], ['Kraken reads assigned', program+' genome fraction (%)', 'Proportion kraken reads mapped with '+program]
    plot_names = ['Kraken reads\nassigned', program+' genome\nfraction(%)', 'Kraken reads\nmapped by\n'+program+' (%)']
    mapping, mins, maxs = [], [], []
    for c in range(len(colormaps)):
      if colnames[c] not in cc_out.columns: continue
      all_c = list(cc_out[colnames[c]].values)
      all_c = [ac for ac in all_c if not np.isnan(ac)]
      min_c, max_c = min(all_c), max(all_c)
      m = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=min_c, vmax=max_c), cmap=colormaps[c])
      ap = mapping.append(m), mins.append(min_c), maxs.append(max_c)
      # if c == 0 and max_c > 100000:
      #   plot_names[0] = 'Kraken reads\nassigned (x10,000)'
      
    xticks, gen_means = [], []
    for s in range(len(plot_order)):
      if coverage_program != 'Both':
        ax_genome = plt.subplot2grid((l,25),(s+1,0), colspan=10)
        ax_identity = plt.subplot2grid((l,25),(s+1,11), colspan=3)
        ax_kraken = plt.subplot2grid((l,25),(s+1,15), colspan=2)
        ax_program_gf = plt.subplot2grid((l,25),(s+1,17), colspan=2)
        ax_program_prop = plt.subplot2grid((l,25),(s+1,19), colspan=2)
      elif program == 'Minimap2':
        ax_genome = plt.subplot2grid((l,50),(s+1,0), colspan=10)
        ax_identity = plt.subplot2grid((l,50),(s+1,21), colspan=3)
        ax_kraken = plt.subplot2grid((l,50),(s+1,28), colspan=2)
        ax_program_gf = plt.subplot2grid((l,50),(s+1,31), colspan=2)
        ax_program_prop = plt.subplot2grid((l,50),(s+1,36), colspan=2)
      elif program == 'Bowtie2':
        ax_genome = plt.subplot2grid((l,50),(s+1,10), colspan=10)
        ax_identity = plt.subplot2grid((l,50),(s+1,24), colspan=3)
        ax_kraken = plt.subplot2grid((l,50),(s+1,28), colspan=2)
        ax_program_gf = plt.subplot2grid((l,50),(s+1,33), colspan=2)
        ax_program_prop = plt.subplot2grid((l,50),(s+1,38), colspan=2)
      axes = [ax_kraken, ax_program_gf, ax_program_prop]
      rounding, div = [None, 3, 3, 3], [None, None, None, None]
      for a in range(len(axes)):
        plot_square(axes[a], mapping[a], cc_out.loc[plot_order[s], colnames[a]], mins[a], maxs[a], rnd=rounding[a], div=div[a])
        if s == 0:
          axes[a].set_title(plot_names[a], fontweight='bold', rotation=90)
      if program == 'Minimap2': ax_genome.set_ylabel(str(plot_order[s])+': $'+cc_out.loc[plot_order[s], 'Species name'].replace('_', '$ $')+'$', fontweight='bold', rotation=0, ha='right', va='center')
      elif coverage_program != 'Both' and program == 'Bowtie2': ax_genome.set_ylabel(str(plot_order[s])+': $'+cc_out.loc[plot_order[s], 'Species name'].replace('_', '$ $')+'$', fontweight='bold', rotation=0, ha='right', va='center')
      xt, gm = plot_genome_coverage(ax_genome, ax_identity, sample, str(plot_order[s]), cc_out.loc[plot_order[s], 'Reference genome length (bp)'], program)
      if xt != 'NA':
        xticks, gen_means = xt, gm
      if s == 0:
        ax_genome.set_title(program+'\nCoverage across genome', fontweight='bold')
        ax_identity.set_title(program+'\nIdentity to\nreference', fontweight='bold')
      if s == len(plot_order)-1:
        plt.sca(ax_identity)
        xt = plt.xticks([70, 80, 90, 100])
        xl = plt.xlabel('Identity (%)')
        plt.sca(ax_genome)
        xl = plt.xlabel('Genome size (mbp)')
      plt.sca(ax_genome)
      if xt == 'NA':
        xt = plt.xlim([0, 1])
        rl = cc_out.loc[plot_order[s], 'Reference genome length (bp)']
        rng = math.ceil(rl/1000000)
        t = plt.xticks([(x*1000000)/rl for x in range(rng)], [str(float(x)) for x in range(rng)])
        plt.sca(ax_identity)
        plt.xticks([])
      else:
        t = plt.xticks(xticks, [round(gen_means[int(x)]/1000000, 1) for x in xticks])
  
    ad = 3
    if coverage_program != 'Both':
      ax_genome_leg = plt.subplot2grid((l,25),((s)+ad,2),colspan=6)
      ax_kraken_leg = plt.subplot2grid((l,100),((s)+ad,61), colspan=6)
      ax_program_gf_leg = plt.subplot2grid((l,100),((s)+ad,69), colspan=6)
      ax_program_prop_leg = plt.subplot2grid((l,100),((s)+ad,77), colspan=6)
      axes = [ax_genome_leg, ax_kraken_leg, ax_program_gf_leg, ax_program_prop_leg]
    elif program == 'Minimap2':
      ax_genome_leg = plt.subplot2grid((l,50),((s)+ad,7),colspan=6)
      ax_kraken_leg = plt.subplot2grid((l,200),((s)+ad,112), colspan=8)
      ax_program_gf_leg = plt.subplot2grid((l,200),((s)+ad,128), colspan=8)
      ax_program_prop_leg = plt.subplot2grid((l,200),((s)+ad,148), colspan=8)
      axes = [ax_genome_leg, ax_kraken_leg, ax_program_gf_leg, ax_program_prop_leg]
    else:
      axes = []
    mapping = [mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap='viridis')]+mapping
    colormaps = ['viridis']+colormaps
    mins, maxs = [0]+mins, [1]+maxs
    labels = ['', 'Reads', 'Genome\nfraction (%)', 'Reads\nmapped (%)']
    for a in range(len(axes)):
      cb = mpl.colorbar.ColorbarBase(axes[a], cmap=colormaps[a], norm=mpl.colors.Normalize(vmin=mins[a], vmax=maxs[a]), orientation='horizontal')
      if a == 0:
        plt.sca(axes[a])
        plt.xticks([0.1, 0.9], ['Less coverage', 'More coverage'])
      else:
        plt.sca(axes[a])
        axes[a].set_xlabel(labels[a])
        if coverage_program == 'Both': plt.xticks([mins[a], maxs[a]], ['Low', 'High'])
        elif a == 1: plt.xticks([mins[a], maxs[a]])

  plt.subplots_adjust(hspace=1.2)
  plt.savefig(save_name+' '+coverage_program+'.png', bbox_inches='tight', dpi=dpi)
  return

if running == 'taxon':
  for taxid in taxid_list:
    fig_save_name = project_folder+'figures/'+taxid
    single_taxon_across_samples(taxid, fig_save_name, project_folder, coverage_program, samples)
else:
  if samples == 'All':
    coverage_out = pd.read_csv(project_folder+'coverage_checker_output.tsv', index_col=0, header=0, sep='\t')
    samples = sorted(list(set(coverage_out.index.values)))
  for sample in samples:
    fig_save_name = project_folder+'figures/'+sample
    multiple_taxa_in_one_sample(fig_save_name, project_folder, taxid_list, sample, top_taxa, sort_by, coverage_program)


