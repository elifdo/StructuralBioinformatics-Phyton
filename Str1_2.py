#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Question: Do the bond lengths between Cα – C, C – N and N– Cα have similar lengths or not?

import numpy as np
import pandas as pd
import math

def dist(x,y):
    return( math.sqrt( ((x[0]-y[0])**2)+ ((x[1]-y[1])**2)+ ((x[2]-y[2])**2)))

#I used 6aa3 PDB id 
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
data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 
    'y_coor','z_coor','occup', 'temp']].astype(float)
data[['atom_num','res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)

# Slicing Cα,C and N atoms
ca = (data.loc[ ( (data["atom_name"]=="CA")),
               ["atom_name","res_seq_num","x_coor","y_coor","z_coor"]  ] )
c = (data.loc[ ( (data["atom_name"]=="C")),
              ["atom_name","res_seq_num","x_coor","y_coor","z_coor"]  ] )
n = (data.loc[ ( (data["atom_name"]=="N")),
              ["atom_name","res_seq_num","x_coor","y_coor","z_coor"]  ] )
 
ca_c_dist = np.empty(len(ca.loc[:,"res_seq_num"]))
c_n_dist = np.empty(len(ca.loc[:,"res_seq_num"]))
n_ca_dist = np.empty(len(ca.loc[:,"res_seq_num"]))
res_seq = np.empty(len(ca.loc[:,"res_seq_num"]))
k=0   
    
for i in  ca.loc[:,"res_seq_num"]:
    ca_coor=ca.loc[ ca["res_seq_num"]==i   ,["x_coor","y_coor","z_coor"]]
    c_coor=c.loc[ c["res_seq_num"]==i   ,["x_coor","y_coor","z_coor"]] 
    n_coor=n.loc[ n["res_seq_num"]==i   ,["x_coor","y_coor","z_coor"]]
    ca_coor = ca_coor.iloc[0][0:3]
    c_coor = c_coor.iloc[0][0:3]
    n_coor = n_coor.iloc[0][0:3]
    ca_c_dist[k] = dist(ca_coor,c_coor )
    c_n_dist[k] = dist(c_coor,n_coor )
    n_ca_dist[k] = dist(n_coor,ca_coor )
    res_seq[k]=i
    k=k+1
    
result = pd.DataFrame(ca_c_dist, columns=["ca_c_dist"])
result =   result.join( pd.DataFrame({'c_n_dist': c_n_dist,
                                      'n_ca_dist': n_ca_dist,           
                                      'res_seq_num': res_seq }) )
 

mean =np.mean(result.loc[:,["ca_c_dist","c_n_dist","n_ca_dist"]])
std =np.std(result.loc[:,["ca_c_dist","c_n_dist","n_ca_dist"]])
mini= np.min(result.loc[:,["ca_c_dist","c_n_dist","n_ca_dist"]])
maxi=  np.max(result.loc[:,["ca_c_dist","c_n_dist","n_ca_dist"]])

result_summary =pd.DataFrame([mean,std,mini,maxi] , index=["mean","std","min","max"] )

# As can be seen, CA-C distances are around 1.52, C-N distances are around 2.45
# and N-CA distances are around 1.46. Their distances are very characteristic, 
# even their ranges don’t coincide with each other. 

# Note: Keep in mind that these statistics are only for one PDB entry. One has
# to be careful while generalizing.
