#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('seaborn-white')
import matplotlib
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
matplotlib.rcParams['figure.figsize'] = (4, 4)
def read_mapstats(fname):
    with open(fname) as f:
        lines = f.readlines()
    df = pd.DataFrame([l.rstrip().split('|\t') for l in lines ]).set_index(0)
    df.index = [i.replace('  ', '') for i in df.index]
    return df


# In[3]:


import os
mapstats = []
# Genome mapped statistics
indir='/home/hsher/sara_clip/GSE146878/mapped_reads'
mapstats = [os.path.join(indir, f) for f in os.listdir(indir) if f.endswith('Log.final.out')]
indir='/home/hsher/sara_clip/202110_fmrp_adaptor/FMRP202109/results/' # new FMRP
mapstats += [os.path.join(indir, f) for f in os.listdir(indir) if f.endswith('repeat-unmapped.sorted.STARLog.final.out')]
indir='/home/hsher/sara_clip/FMRP_old/FMRP_old/results/' # old FMRP
mapstats += [os.path.join(indir, f) for f in os.listdir(indir) if f.endswith('repeat-unmapped.sorted.STARLog.final.out')]


# In[5]:


# repeat mapped statistics
indir='/home/hsher/sara_clip/202110_fmrp_adaptor/FMRP202109/results/' # new FMRP
mapstats += [os.path.join(indir, f) for f in os.listdir(indir) if 'repeat-unmapped' not in f and 'final.out' in f]
indir='/home/hsher/sara_clip/FMRP_old/FMRP_old/results/' # old FMRP
mapstats += [os.path.join(indir, f) for f in os.listdir(indir) if 'repeat-unmapped' not in f and 'final.out' in f]


# In[6]:


def MapStats():
    all_df=[]
    for f in mapstats:
        df = read_mapstats(os.path.join(indir, f))
        df.columns = [os.path.basename(f)]
        all_df.append(df)
    
    full_df = pd.concat(all_df, axis = 1)
    needed_rows = [
        "Uniquely mapped reads % ",
        " % of reads unmapped: too short ",
        " % of reads mapped to multiple loci " ,]
    mapping_stats=full_df.loc[needed_rows]
    return all_df, mapping_stats


# In[7]:


# cleaning dataframe
total = pd.concat(MapStats()[0], axis = 1).T
assert total is not None
for col in total.columns:
    if 'Number' or 'Average' in col:
        try:
            total[col] = total[col].astype(float)
        except Exception as e:
            # print(col, e)
    if '%' in col:
        total[col] = total[col].str.replace('%', '').astype(float)


# In[8]:


MapStats()[1]


# In[1]:


# generate metadata
total['experiment'] = ['old' if 'FMRP_old' in i else 'new' if 'FMRP_202109' in i else 'Organoid' for i in total.index]
total['library'] = ['CLIP' if 'CLIP' in i else 'INPUT' if 'INPUT' in i else 'IGG' for i in total.index]
total['rep'] = ['rep1' if 'CLIP1' in i or 'INPUT1' in i else 'rep2' if 'INPUT2' in i or 'CLIP2' in i else 'rep1' for i in total.index]
total['Ab'] = ['Abcam' if 'Abcam' in i else 'MBL' if 'MBL' in i else 'Other' for i in total.index]
total['RBP'] = ['FXR1' if 'FXR1' in i.split(".")[1].split("_")[0] else 'FMRP' for i in total.index]
total['Cell'] = ['SC' if 'SC' in i else 'Neu' if 'Neu' in i else 'Organoid' for i in total.index]
total['genotype'] = ['WT' if 'WT' in i else 'FXRKO' if 'FXRKO' in i else 'FMRPKO' if 'KO' in i else 'WT' for i in total.index]
total['reference'] = ['genome' if 'repeat-unmap' in i else 'repeat' if 'STAR' in i else 'genome' for i in total.index]


# # Make Table S1, statistic for eCLIP data

# In[11]:


columns = ['Number of input reads ', 'Average input read length ', ' Uniquely mapped reads number ', 'Uniquely mapped reads % ', ' % of reads mapped to multiple loci ',' % of reads mapped to too many loci ', ' % of reads unmapped: too short ',]
total.loc[total['experiment'].isin(['old', 'new'])].pivot_table(index = ['RBP','experiment', 'library', 'Ab', 'Cell', 'genotype', 'rep'], columns = 'reference', values = columns).to_csv('TableS1.csv')


# In[12]:


total.to_csv('mapping_stat.csv') 


# In[13]:



fg = sns.FacetGrid(data=total.loc[total['reference']=='genome'], hue='experiment', hue_order=['old', 'new', 'Organoid'], aspect=1.61)
fg.map(plt.scatter, 'Number of input reads ', ' Uniquely mapped reads number ').add_legend()


# In[14]:


total.loc[(total['Uniquely mapped reads % ']<40)&(total['reference']=='genome')&(total['RBP']=='FXR1')]


# In[15]:


total.plot.bar(y = columns, stacked = True)


# In[16]:


total.boxplot(by = 'experiment', column = 'Number of input reads ')
plt.ylabel('Number of input reads')
plt.tight_layout()


# In[17]:


total.boxplot(by = 'experiment', column = 'Uniquely mapped reads % ')
plt.ylabel('Uniquely mapped reads % ')
plt.tight_layout()


# In[18]:



total.boxplot(by = ['Cell','experiment', 'library'], column = 'Uniquely mapped reads % ')
plt.xticks(rotation = 90)


# In[19]:


total.boxplot(by = ['Cell', 'library'], column = 'Uniquely mapped reads % ')
plt.xticks(rotation = 90)


# In[20]:


total.boxplot(by = ['Ab', 'library'], column = 'Uniquely mapped reads % ')
plt.xticks(rotation = 90)


# In[21]:


total.boxplot(by = ['WT', 'library'], column = 'Uniquely mapped reads % ')
plt.xticks(rotation = 90)


# In[9]:


total.columns


# # Quality metrices

# In[17]:


# find all entropy
indir = '/home/hsher/sara_clip/202110_fmrp_adaptor/FMRP202109/results/'
entropies_path = [os.path.join(indir,f) for f in os.listdir(indir) if f.endswith('entropynum')]

indir='/home/hsher/sara_clip/FMRP_old/FMRP_old/results/'
entropies_path += [os.path.join(indir,f) for f in os.listdir(indir) if f.endswith('entropynum')]


# In[18]:


data = []
for file in entropies_path:
    with open(file) as f:
        entropy=float(f.readlines()[0])
    data.append([os.path.basename(file), entropy])   


# In[19]:


df = pd.DataFrame(data, columns = ['fname', 'entropy'])
df.set_index('fname', inplace = True)


# In[20]:


df['experiment'] = ['old' if 'FMRP_old' in i else 'new' if 'FMRP_202109' in i else 'Organoid' for i in df.index]
df['library'] = ['CLIP' if 'CLIP' in i else 'INPUT' if 'INPUT' in i else 'IGG' for i in df.index]
df['Ab'] = ['Abcam' if 'Abcam' in i else 'MBL' if 'MBL' in i else 'Other' for i in df.index]
df['Cell'] = ['SC' if 'SC' in i else 'Neu' if 'Neu' in i else 'Organoid' for i in df.index]
df['WT'] = ['WT' if 'WT' in i else 'FXRKO' if 'FXRKO' in i else 'FMRPKO' if 'KO' in i else 'WT' for i in df.index]


# In[21]:


df.head()


# In[22]:


df.boxplot(by = ['experiment'], column = 'entropy')


# In[23]:


df.boxplot(by = ['experiment', 'Cell'], column = 'entropy')


# In[24]:


df.boxplot(by = ['experiment', 'WT'], column = 'entropy')
plt.xticks(rotation = 90)


# In[3]:


def main(clip_data):
    res = MapStats(clip_data)
    return res

if __name__ == '__main__':
    script = sys.argv[0]
    filename = sys.argv[1 :]
    clip_data = numpy.loadtxt(filename, delimiter=',')
    main(clip_data)
        

