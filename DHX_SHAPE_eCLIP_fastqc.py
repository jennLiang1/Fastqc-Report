#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
indir='/home/hsher/scratch/20210922_dhx36shape/fshape_eclip_pipe/results/'
# see if what is left
trimmed_fastqc_txt = [os.path.join(indir,f) for f in os.listdir(indir) if f.endswith('.sorted.fq.fastqc_data.txt')]


# In[2]:


trimmed_fastqc_txt


# In[3]:


def parse_fastqc(filename):
    with open(filename) as f:
        module_lines = []
        psall = []
        all_modules = {}
        for line in f:
            
            if '>>END_MODULE' in line:
                all_modules[name]= module_lines
                psall.append([name, passfail])
                module_lines = []
                
            elif line.startswith('>>'): # start of new block
                
                name = line.split('\t')[0].replace('>>', '')
                passfail = line.split('\t')[1].rstrip()
                
            else:
                values = line.rstrip().split('\t')
                module_lines.append(values)
        
        all_modules[name]= module_lines
        psall.append([name, passfail])
        
        
        # make into dataframe
        for item in all_modules.keys():
            for i, line in enumerate(all_modules[item]):
                if '#' not in line[0]:
                    break # first line that is not column
            try:
                df = pd.DataFrame(all_modules[item][i:], columns = all_modules[item][i-1])
            except:
                df = pd.DataFrame(all_modules[item][i:])
            
            all_modules[item] = df
        
        return all_modules, psall
                


# In[4]:


all_modules,pass_fail = parse_fastqc('/home/hsher/scratch/20210922_dhx36shape/fshape_eclip_pipe/results/fSHAPE_under_eclip_pipeline.DHX36_D2_IP.umi.r1.fqTrTr.sorted.fq.fastqc_data.txt')


# In[5]:


def Fastqc_Res():
    all_data = [] # key -> value
    for file in trimmed_fastqc_txt:
        all_modules,pass_fail = parse_fastqc(file)
    
    # add filename
    #all_data[file] = pass_fail
    
    df = pd.DataFrame(pass_fail).set_index(0)
    df.columns = [file]
    
    all_data.append(df)
    return all_data


# In[29]:


pd.concat(all_data, axis = 1).shape


# In[2]:


df = pd.DataFrame(pass_fail).set_index(0)
df.columns = [file]


# In[23]:


df


# In[6]:


all_modules.keys() # dictionary key --> value


# In[8]:


pass_fail # list


# In[7]:


all_modules['Basic Statistics'] # a pandas dataframe


# In[12]:


all_modules['Per sequence GC content']


# In[10]:


all_modules['Overrepresented sequences']


# In[11]:


all_modules['Adapter Content']


# # screen adapters

# In[ ]:


from Bio import SeqIO
adaptor_path = '/projects/ps-yeolab4/software/eclip/0.7.0/examples/inputs/'
adaptors = [os.path.join(adaptor_path,f) for f in os.listdir(adaptor_path) if f.startswith('Inv')]

adaptor_seq_dict = {}
for filename in adaptors:
    sequences = []
    for record in SeqIO.parse(filename, "fasta"):
        sequences.append(str(record.seq))
    adaptor_seq_dict[os.path.basename(filename)] = sequences
    


# In[ ]:


# to one sequence
one_seq_dict = {}
for adaptor in adaptor_seq_dict.keys():
    fragments = adaptor_seq_dict[adaptor]
    
    for i,f in enumerate(fragments):
        if i == 0:
            whole_sequence = f
        else:
            whole_sequence += f[-1]
    one_seq_dict[adaptor] = whole_sequence


# In[ ]:


one_seq_dict


# In[ ]:


a['Overrepresented sequences']


# In[ ]:


Seq(unique_part[adaptor]).reverse_complement()


# In[ ]:


from Bio.Seq import Seq
df = a['Overrepresented sequences']

for adaptor in one_seq_dict.keys():
    if len(one_seq_dict[adaptor])>0:
        df[adaptor] = df['#Sequence'].str.contains(one_seq_dict[adaptor]) | df['#Sequence'].str.contains(str(Seq(one_seq_dict[adaptor]).reverse_complement()))
    


# In[ ]:


df


# In[ ]:


Seq('GGTCTACGGCCATACCACCCTGAACGCGCCCGATCTCGTCTGATC').reverse_complement()


# In[ ]:


one_seq_dict


# In[ ]:


len('NNAGGTGCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')


# In[ ]:


Seq('CGGCTATGATCTCGTATGCC').reverse_complement()


# In[ ]:


def main():
    script = sys.argv[0]
    filename = sys.argv[1]
    data = numpy.loadtxt(filename, delimiter=',')
    for row_mean in numpy.mean(data, axis=1):
        print(row_mean)
    Fastqc_Res ()        

if __name__ == '__main__':
   main()
        

