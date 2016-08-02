
# coding: utf-8

# In[36]:
from ggplot import *
import numpy as np
import pandas as pd

from generateTrees import gen_trees
from tools import flatten, timeit


# In[ ]:

c_list = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
Ne = 10000
output = np.empty((len(c_list), 3))
mean = lambda x: sum(x)/len(x)

for j, c in enumerate(c_list):
    treelist = lambda: gen_trees(n_sp_trees=5,
                                 n_gene_trees=10,
                                 n_sp=5,
                                 n_ind=20,
                                 sp_depth=Ne*c,
                                 Ne=Ne)[1]

    treelist = timeit(treelist, "making trees")
    interdist, intradist = [], []
    for k in range(len(treelist)):
        for i, tree in enumerate(treelist[k]):
            intradist.append([])
            interdist.append([])
            pdm = tree.phylogenetic_distance_matrix()
            for t1, t2 in pdm.distinct_taxon_pair_iter():
                dist = pdm.distance(t1, t2)
                if t1.label[0] == t2.label[0]:
                    intradist[i].append(dist)
                else:
                    interdist[i].append(dist)

    # inter_means = list(map(mean, interdist))
    # intra_means = list(map(mean, intradist))
    inter_mean = mean(list(flatten(*interdist)))
    intra_mean = mean(list(flatten(*intradist)))
    output[j] = [c, inter_mean, intra_mean]

cols = ["c", "interspecies mean", "intraspecies mean"]
df = pd.DataFrame(output, columns=cols)
mdf = pd.melt(df, id_vars=['c'])


# In[ ]:

title="Mean Interspecies and Intraspecies Distances\n" + \
      "5 Species Trees with 10 Gene Trees\n" + \
      "5 Species per Tree, 20 Individuals per Species"

p = ggplot(mdf, aes(x='c', y='value', color='variable', group='variable'))
p += geom_line() 
p += labs(x="c ratio", 
          y="Distance", 
          title=title)
p += scale_x_continuous(breaks=c_list, labels=c_list)
print(p)

# In[ ]:



