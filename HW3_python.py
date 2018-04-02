
# coding: utf-8

# In[132]:


from Bio import SeqIO 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import spline


class Kmer_spectrum:
    k = 0
    q = 0
    kmer_dict = {}
    kmers_list = []
    L = 0
    
    def __init__(self, k, q): 
        self.k = k
        self.q = q
    
    def dict_all_kmers (self, fastq):
        for record in SeqIO.parse(fastq, "fastq"):
            for index in range(len(record.seq)-self.k+1):
                flag = 0
                current_kmer = record.seq[index:(index+self.k)]
                current_quality = record.letter_annotations["phred_quality"][index:(index+self.k)]
                for val in current_quality:
                    if val < self.q:
                        flag = 1
                        break
                if flag == 0:
                    if current_kmer in self.kmer_dict:
                        self.kmer_dict[current_kmer] += 1 
                    else:
                        self.kmer_dict[current_kmer] = 1 
        return (self.kmer_dict)

    
    def list_kmers (self):
        new_dict = {}
        for key in self.kmer_dict:
            if self.kmer_dict[key] in new_dict:
                new_dict[self.kmer_dict[key]] += 1
            else:
                new_dict[self.kmer_dict[key]] = 1
        for key in new_dict:
            self.kmers_list += [[key, new_dict[key]]]
        return(self.kmers_list)
                    
    def visualise(self, xlim, ylim):
        for i in self.kmers_list:
            plt.plot(i[0], i[1], color='blue', marker='o')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel('Kmer multiplicity')
        plt.ylabel('Kmer count')
        plt.title('Kmer spectrum')
        plt.show()

    def genome_length(self, cut):
        for i in self.kmers_list:
            if i[0] >= cut:
                self.L += i[0]*i[1]/self.k
        return(self.L)

# In[133]:


with open ("/Users/yukornienko/Downloads/test_kmer.fastq") as fastq:
    kmers = Kmer_spectrum(12, 35)
    kmers_dict = kmers.dict_all_kmers(fastq)
    kmers_list = kmers.list_kmers()
    print(kmers_list)
    kmers.visualise([1, 150], [0, 100000])
    
# In[134]:

L = kmers.genome_length(30)
print(L)

# In[151]:

with open ("/Users/yukornienko/Downloads/test_kmer.fastq") as fastq:
    kmers_2 = Kmer_spectrum(14, 1)
    kmers_dict_2 = kmers_2.dict_all_kmers(fastq)
    kmers_list_2 = kmers_2.list_kmers()

    
# In[160]:
kmers_2.visualise([1, 400], [0, 150000])

# In[161]:
L = kmers_2.genome_length(40)
print(L)

