#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Question: Compute the phi and psi dihedral angles for Cα[i]
# when given the coordinates of C[i-1], N[i], Cα[i],C[i] and N[i+1]
# and plot the Ramachandran plot.
# Plot Proline and Glycine residues seperately.
# Also, plot the most unstable residues with the largest 5 B-factors(temprature)

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

#Function that finds phi and psi angles given the corresponding coordinates 
#of the related atoms. 

def phi_psi(ci_1,n_i,ca_i,c_i,n_ipl1) :
    u=n_i-ci_1
    v=ca_i-n_i
    w=c_i-ca_i
    norm1= np.cross(u,v)
    norm1= norm1/np.linalg.norm(norm1)
    norm2= np.cross(v,w)
    norm2= norm2/np.linalg.norm(norm2)
    if np.dot(norm1,w)<1 : sign= -1
    else : sign=1
    phi= sign *  np.arccos(np.dot(norm1,norm2))
    u=ca_i-n_i
    v=c_i-ca_i
    w=n_ipl1-c_i
    norm1= np.cross(u,v)
    norm1= norm1/np.linalg.norm(norm1)
    norm2= np.cross(v,w)
    norm2= norm2/np.linalg.norm(norm2)
    if np.dot(norm1,w)<1 : sign= -1
    else : sign=1
    psi=   np.arccos(np.dot(norm1,norm2))
    return(phi* 180 / math.pi,psi* 180 / math.pi)

# I used the structure with pdb id 6aa3

with open('6aa3.pdb') as pdbfile:
    lines= list()
    for line in pdbfile:
        if line[:4] == 'ATOM' :
            splitted_line = [line[:7], line[7:11], line[11:16],line[16:20], 
                             line[20:22], line[22:26], line[26:38],line[38:46],
                             line[46:54], line[54:60], line[60:66],line[66:]]
            splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
            lines.append(splitted_line)

data = pd.DataFrame(lines, columns=["atom","atom_num","atom_name","res_name",
                                    "chain","res_seq_num","x_coor","y_coor",
                                    "z_coor","occup","temp","elt_sym"])
data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 'y_coor',
    'z_coor','occup', 'temp']].astype(float)
data[['atom_num', 'res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)


dihedrals= pd.DataFrame(index=range(np.min(data.loc[:,"res_seq_num"])+1,
                                    np.max(data.loc[:,"res_seq_num"])-1 ),
                columns=["phi","psi"]  )
x= pd.Series()
y=pd.Series()

#Calculating dihedrals for each residue
for i in range(np.min(data.loc[:,"res_seq_num"])+1 ,np.max(data.loc[:,"res_seq_num"])) :
    c_i=data.loc[  (data.loc[:,"res_seq_num"]==i) & (data.loc[:,"atom_name"]=="C")  ,
                 ["x_coor","y_coor","z_coor"]]  
    ci_1=data.loc[  (data.loc[:,"res_seq_num"]==(i-1)) & (data.loc[:,"atom_name"]=="C")  ,
                  ["x_coor","y_coor","z_coor"]]  
    n_i=data.loc[  (data.loc[:,"res_seq_num"]==i) & (data.loc[:,"atom_name"]=="N")  ,
                 ["x_coor","y_coor","z_coor"]]  
    ca_i=data.loc[  (data.loc[:,"res_seq_num"]==i) & (data.loc[:,"atom_name"]=="CA")  , 
                  ["x_coor","y_coor","z_coor"]]  
    n_ipl1=data.loc[  (data.loc[:,"res_seq_num"]==(i+1)) & (data.loc[:,"atom_name"]=="N"),
                    ["x_coor","y_coor","z_coor"]]     
    output =phi_psi(np.array(ci_1.iloc[0]),np.array(n_i.iloc[0]),np.array(ca_i.iloc[0]),
                    np.array(c_i.iloc[0]),np.array(n_ipl1.iloc[0]))
    dihedrals.loc[i,["phi","psi"]]=output
    x= x.append(data.loc[ (data.loc[:,"res_seq_num"]==i) & (data.loc[:,"atom_name"]=="C"), 
                         "res_name" ]) 
    y= y.append(data.loc[ (data.loc[:,"res_seq_num"]==i) & (data.loc[:,"atom_name"]=="C"),
                         "res_seq_num" ]) 

x= x.to_frame(name='res_name')
x=x.reset_index(drop=True)
y= y.to_frame(name='res_num')
y=y.reset_index(drop=True)
dihedrals= dihedrals.reset_index(drop=True)
dihedrals= dihedrals.join([x,y])

# Ramachandran plot

plt.axis('scaled')
plt.axis([-180, 180, -180, 180])
plt.plot(dihedrals.loc[:,"phi"], dihedrals.loc[:,"psi"],'ro',markersize=2)
plt.grid(True)

plt.savefig('ramachandran.png')

# Ramachandran plot of Proline residues only

pro_num = len(dihedrals.loc[  dihedrals.loc[:,"res_name"] =="PRO","phi"])

plt.axis('scaled')
plt.axis([-180, 180, -180, 180])
plt.plot(dihedrals.loc[  dihedrals.loc[:,"res_name"] =="PRO","phi"], 
         dihedrals.loc[ dihedrals.loc[:,"res_name"] =="PRO","psi"],'ro',markersize=2)
plt.grid(True)
plt.savefig('ramachandran_pro.png')

# Ramachandran plot of Glycine residues only

gly_num = len(dihedrals.loc[  dihedrals.loc[:,"res_name"] =="GLY","phi"])

plt.axis('scaled')
plt.axis([-180, 180, -180, 180])
plt.plot(dihedrals.loc[  dihedrals.loc[:,"res_name"] =="GLY","phi"], 
         dihedrals.loc[ dihedrals.loc[:,"res_name"] =="GLY","psi"],'ro',markersize=2)
plt.grid(True)
plt.savefig('ramachandran_gly.png')

#Parsing the B-factors of Cα atoms in the PDB file. Finding the most flexible 5 
#residues and plot them.

y = data.loc[   (data.loc[:,"atom_name"]=="CA")  , ["temp","res_seq_num"]  ]
y= y.sort_values(["temp","res_seq_num"],ascending=[False,True])
y=y.reset_index(drop=True)
res =y.loc[0:5,"res_seq_num"]
unstable= dihedrals.loc[  dihedrals['res_num'].isin(res)  ,:   ]

plt.axis('scaled')
plt.axis([-180, 180, -180, 180])
plt.plot(unstable.loc[:,"phi"], unstable.loc[:,"psi"],'ro',markersize=2)
plt.grid(True)
plt.savefig('ramachandran_unstable.png')

# Phi angles of prolines in our plot is around -60, which is expected because 
#prolines have a very restricted possible region. Psi and phi angles of glycine,
# on the other hand, seem more flexible.  