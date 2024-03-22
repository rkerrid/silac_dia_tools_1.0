# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:27:12 2024

@author: robbi
"""
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# read in light 
path = 'G:/My Drive/Data/main experiments/20240219 baby benchmark for pydia_sis/protein_groups/'
light = pd.read_csv(f'{path}light.csv', sep=',')
light = light[['Protein.Group', 'S1_1', 'S2_1', 'S3_1']]
light.columns = ['Protein.Group', 'S1', 'S2', 'S3']

# import href for intensity col
href = pd.read_csv(f'{path}href.csv', sep=',')
href = href[['Protein.Group', 'S1_1']]
href.columns = ['Protein.Group', 'Intensity']
light = pd.merge(href, light, on='Protein.Group', how = 'inner')

# get ratios for S1/S2, S1/S3, S2/S3 (choose rep or get median value)
light['S1/S2'] = np.log2(light['S1']/light['S2'])
light['S1/S3'] = np.log2(light['S1']/light['S3'])
light['S2/S3'] = np.log2(light['S2']/light['S3'])
light['Intensity'] = np.log10(light['Intensity'])

# subset human and ecoli
human = light[light['Protein.Group'].str.contains('HUMAN_')]
ecoli = light[light['Protein.Group'].str.contains('ECOLI_')]

human['Species'] = 'Human'
ecoli['Species'] = 'E. coli'

# Combine human and ecoli into a single DataFrame
light = pd.concat([human, ecoli])

# log2 ratios between samples
human_expected = 0
ecoli_s1s2 = np.log2(1/10)
ecoli_s1s3 = np.log2(1/50)
ecoli_s2s3 = np.log2(10/50)


# plot intensity on x and ratio on y for each ratio col and eco/human subset
sns.scatterplot(light, x='Intensity', y='S1/S2', hue='Species')
plt.legend(title='Species')


# Draw horizontal lines for expected ratios
plt.hlines(human_expected, xmin=2, xmax=6, colors='blue', linestyles='dashed', label='Human Expected')
plt.hlines(ecoli_s1s2, xmin=2, xmax=6, colors='red', linestyles='dashed', label='E. coli S1/S2 Expected')

plt.title('Human and Ecoli ecpected ratios for S1:S2')
plt.show()


# plot intensity on x and ratio on y for each ratio col and eco/human subset
sns.scatterplot(light, x='Intensity', y='S1/S3', hue='Species')
plt.legend(title='Species')


# Draw horizontal lines for expected ratios
plt.hlines(human_expected, xmin=2, xmax=6, colors='blue', linestyles='dashed', label='Human Expected')
plt.hlines(ecoli_s1s3, xmin=2, xmax=6, colors='red', linestyles='dashed', label='E. coli S1/S3 Expected')

plt.title('Human and Ecoli ecpected ratios for S1:S3')
plt.show()


# plot intensity on x and ratio on y for each ratio col and eco/human subset
sns.scatterplot(light, x='Intensity', y='S2/S3', hue='Species')
plt.legend(title='Species')


# Draw horizontal lines for expected ratios
plt.hlines(human_expected, xmin=2, xmax=6, colors='blue', linestyles='dashed', label='Human Expected')
plt.hlines(ecoli_s2s3, xmin=2, xmax=6, colors='red', linestyles='dashed', label='E. coli S2/S3 Expected')

plt.title('Human and Ecoli ecpected ratios for S2:S3')
plt.show()


