from __future__ import print_function
import mdtraj as md
import pandas as pd
import numpy as np
import sklearn
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
pca1 = PCA(n_components=2)
#
protein_a = md.load('mdfitdt_protein_a.xtc', top='protein_a.pdb')
protein_a.superpose(a, 0)
protein_a_reduced_cartesian = pca1.fit_transform(a.xyz.reshape(a.n_frames, a.n_atoms * 3))

protein_b = md.load('mdfitdt_protein_b.xtc', top='protein_b.pdb')
protein_b.superpose(b, 0)
protein_b_reduced_cartesian = pca1.fit_transform(b.xyz.reshape(b.n_frames, b.n_atoms * 3))

protein_c = md.load('mdfitdt_protein_c.xtc', top='protein_c.pdb')
protein_c.superpose(c, 0)
protein_c_reduced_cartesian = pca1.fit_transform(c.xyz.reshape(c.n_frames, c.n_atoms * 3))

ax=list(protein_a_reduced_cartesian[:,0])
ay=list(protein_a_reduced_cartesian[:,1])
bx=list(protein_b_reduced_cartesian[:,0])
by=list(protein_b_reduced_cartesian[:,1])
cx=list(protein_c_reduced_cartesian[:,0])
cy=list(protein_c_reduced_cartesian[:,1])

dfx = [*ax, *bx, *cx]
dfy = [*ay, *by, *cy]

at= ["protein"]*len(ax)
bt= ["protein-lig(Noncovalent)"]*len(ax)
ct= ["protein_lig(Covalent)"]*len(ax)

dft_type = [*at, *bt, *ct]
print(len(dfx))
print(len(dfy))
print(len(dft_type))
df = pd.DataFrame({'pcax':dfx, 'pcay':dfy, 'Sample_Type':dft_type})
print(df)

sns.jointplot(data=df,x="pcax",fill=True, kind="kde",y="pcay",hue="Sample_Type",alpha= 0.6)
# sns.kdeplot(data=df,x="pcax", y="pcay")
plt.show()








