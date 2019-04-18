#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

# Question: Translate proteins such that centroid is at the origin and check 
# if the sum of the coordinates are almost zero or not
# Find an NMR structure and compare RMSD of the models with the first model
# Align the structures 1xg5 chain A, 3gy0 chain A and 1edo by using MULTIPROT 
# Write a code to calculate the RMSD value 

# Parsing atoms of the 6aa3 protein from pdb file and store as a dataframe

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

# Function that finds centroid of a protein given a dataframe of x,y,z coordinates

def centroid_finder(coordinates):
    sum = coordinates.sum(axis=0)
    centroid=sum/coordinates.shape[0]   
    return(centroid)
    
# Slicing coordinates of the protein and finding the centroid   
coordinates =data[["x_coor","y_coor","z_coor"]]
centroid = centroid_finder(coordinates)

# Subtracting centroid from each point to carry the centroid to the origin

new_coordinates= coordinates.apply(lambda x: x-centroid, axis=1) 

# Checking if the new sum of the x,y,z coordinates are almost zero
new_sum = new_coordinates.sum(axis=0)
check = (abs(new_sum)<[10**(-10),10**(-10),10**(-10)])
check

# Function finding the distance of two points with x,y,z coordinates
def dist(a,b):
    return( math.sqrt( ((a[0]-b[0])**2)+ ((a[1]-b[1])**2)+ ((a[2]-b[2])**2))  )
    
# Function finds RMSD value of two aligned proteins given coordinates of them
def rmsd(model1_coord, model2_coord) :
    square_distances= np.empty(shape=model1_coord.shape[0])
    for i in range(model1_coord.shape[0]):
        a= model1_coord.iloc[i]
        b= model2_coord.iloc[i]
        square_distances[i]=((dist(a,b))**2)
    rmsd=np.sqrt(np.mean(square_distances))
    return (rmsd)

# Parsing an NMR structure of protein 6agp, with multiple models
with open('6agp.pdb') as pdbfile:
    lines= list()
    models=list()
    i=1
    for line in pdbfile:
        if (line[:4] == 'ATOM' and line[13:15]=='CA' )  :
            splitted_line = [line[:7], line[7:11], line[11:16],line[16:20], 
                             line[20:22], line[22:26], line[26:38],line[38:46],
                             line[46:54], line[54:60], line[60:66],line[66:]]
            splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
            lines.append(splitted_line)
        if line[:6] == 'ENDMDL' :
            data = pd.DataFrame(lines, columns=["atom","atom_num","atom_name","res_name",
                                    "chain","res_seq_num","x_coor","y_coor",
                                    "z_coor","occup","temp","elt_sym"])
            data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 
                'y_coor','z_coor','occup', 'temp']].astype(float)
            data[['atom_num','res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)
            models.append(data)
            i=i+1
            lines=list()

#Finding RMSD values of the models with the first model and plotting
            
rmsd_values= np.empty(shape=len(models)-1)

for i in range(1,len(models)):
     rmsd_values[(i-1)]=rmsd(models[i][["x_coor","y_coor","z_coor"]],
               models[0][["x_coor","y_coor","z_coor"]] )

y_pos = np.arange(len(np.arange(2,21)))
plt.bar(y_pos, rmsd_values, align='center', alpha=0.5)
plt.xticks(y_pos, np.arange(2,21))
plt.ylabel('RMSD')
plt.title('Models vs RMSD values with the first model')
plt.savefig('modelvsrmsd.png')

# Functions finds RMSD for more than two models 
# Multiple RMSD calculated by mean of the pairwise rmsds 
def rmsd_mult(models):
    r=len(models)
    rmsd_values=list()   
    for i in range(1,(r)) :
        j=0
        while True :
           rmsd_values.append(rmsd(models[i][["x_coor","y_coor","z_coor"]],
                                        models[j][["x_coor","y_coor","z_coor"]]) )            
           if j==(i-1) :
               break
           j=j+1
    return(np.mean(pd.Series(rmsd_values)))  
           
# Multiple RMSD calculated by mean of the rmsds of models from artifical mean model 
def rmsd_mult2(models):
    [meanx,meany,meanz]=[pd.Series([0]* models[1].shape[0]),
                     pd.Series([0]* models[1].shape[0]),
                     pd.Series([0]* models[1].shape[0])]
    for i in range(0,(len(models))):
        meanx=meanx.add(models[i]["x_coor"])
        meany=meany.add(models[i]["y_coor"])
        meanz=meanz.add(models[i]["z_coor"])
    [meanx,meany,meanz]=[meanx/len(models),meany/len(models),meanz/len(models)]
    meanxyz= pd.DataFrame(data={'meanx': meanx, 'meany': meany,'meanz': meanz})
    rmsd_values= np.empty(shape=len(models))
    for i in range(0,len(models)):
        rmsd_values[i]=rmsd(models[i][["x_coor","y_coor","z_coor"]], meanxyz)
    return(np.mean(rmsd_values))
    
# CALCULATING RMSD FOR TWO PROTEINS 1edo and 1xg5A, alignment are done with 
#    MULTIPROT server
    
# reading aligned pdb file of two proteins 
    
with open('edoxg5.pdb') as pdbfile:
    lines= list()
    edoxg5=list()
    i=1
    for line in pdbfile:
        if (line[:4] == 'ATOM' and line[13:15]=='CA' )  :
            splitted_line = [line[:7], line[7:11], line[11:16],line[16:20], 
                             line[20:22], line[22:26], line[26:38],line[38:46],
                             line[46:54], line[54:60], line[60:66],line[66:]]
            splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
            lines.append(splitted_line)
        if line[:6] == 'ENDMDL' :
            data = pd.DataFrame(lines, columns=["atom","atom_num","atom_name",
                                                "res_name", "chain",
                                            "res_seq_num","x_coor","y_coor",
                                           "z_coor","occup","temp","elt_sym"])
            data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 
                'y_coor','z_coor','occup', 'temp']].astype(float)
            data[['atom_num','res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)
            edoxg5.append(data)
            i=i+1
            lines=list()
edo_xg5_aligned= edoxg5[0]             
xg5_edo_aligned=edoxg5[1]

#information of which residues are aligned 

with open('edoxg5.txt') as file:
    lines= list()
    for line in file:
        splitted_line = line.split()
        splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
        lines.append(splitted_line)
      
data1=list()
data2=list()
for i in range(1,int((len(lines))/2) ):
    data1.append(lines[2*i])
    data2.append(lines[(2*i)+1])
data1 = pd.DataFrame(data1, columns=["chain","aa_name","residue"])
data1[['residue']] = data1[['residue']].astype(int)
data2 = pd.DataFrame(data2, columns=["chain","aa_name","residue"])
data2[['residue']] = data2[['residue']].astype(int)

# coordinates of the aligned residues
edo_coord = edo_xg5_aligned.loc[ edo_xg5_aligned["res_seq_num"].isin(data1[['residue']].squeeze()  ),
                            ["x_coor","y_coor","z_coor","occup"]  ]
edo_coord=edo_coord.reset_index(drop=True)
xg5_coord= xg5_edo_aligned.loc[ xg5_edo_aligned["res_seq_num"].isin(data2[['residue']].squeeze()  ),
                            ["x_coor","y_coor","z_coor","occup"]  ]
xg5_coord= xg5_coord.reset_index(drop=True)

# residues with multiple entries since their occupancy values are 0.5
extras= list()

for i in range(1,xg5_coord.shape[0]):
    if(xg5_coord["occup"].iloc[i]==0.5 ) :
        if (i-1) not in extras:
            extras.append(i)

for i in range(0,(len(extras))):
    newrow=xg5_coord.loc[[extras[i], extras[i]+1 ],["x_coor","y_coor","z_coor"]].mean(axis=0)
    xg5_coord.loc[ extras[i]+1 ,["x_coor","y_coor","z_coor"]]=newrow
    
xg5_coord=xg5_coord.drop(extras ,axis=0)
xg5_coord=xg5_coord.reset_index(drop=True)


rmsd(edo_coord[["x_coor","y_coor","z_coor"]],xg5_coord[["x_coor","y_coor","z_coor"]])


# CALCULATING RMSD FOR TWO PROTEINS 1edo, 1xg5A and 3gy0A

# reading aligned pdb file of three proteins 
    
with open('1edo1xg5A3gy0A.pdb') as pdbfile:
    lines= list()
    proteins=list()
    i=1
    for line in pdbfile:
        if (line[:4] == 'ATOM' and line[13:15]=='CA' )  :
            splitted_line = [line[:7], line[7:11], line[11:16],line[16:20], 
                             line[20:22], line[22:26], line[26:38],line[38:46],
                             line[46:54], line[54:60], line[60:66],line[66:]]
            splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
            lines.append(splitted_line)
        if line[:6] == 'ENDMDL' :
            data = pd.DataFrame(lines, columns=["atom","atom_num","atom_name","res_name",
                                    "chain","res_seq_num","x_coor","y_coor",
                                    "z_coor","occup","temp","elt_sym"])
            data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 
                'y_coor','z_coor','occup', 'temp']].astype(float)
            data[['atom_num','res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)
            proteins.append(data)
            i=i+1
            lines=list()
edo_triple = proteins[0]             
xg5A_triple=proteins[1]
gy0A_triple =proteins[2]

#information of which residues are aligned 

with open('1edo1xg5A3gy0A.txt') as file:
    lines= list()
    for line in file:
        splitted_line = line.split()
        splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
        lines.append(splitted_line)
      
data1=list()
data2=list()
data3=list()
for i in range(1,int((len(lines))/3) ):
    data1.append(lines[3*i])
    data2.append(lines[(3*i)+1])
    data3.append(lines[(3*i)+2])
data1 = pd.DataFrame(data1, columns=["chain","aa_name","residue"])
data1[['residue']] = data1[['residue']].astype(int)
data2 = pd.DataFrame(data2, columns=["chain","aa_name","residue"])
data2[['residue']] = data2[['residue']].astype(int)
data3 = pd.DataFrame(data3, columns=["chain","aa_name","residue"])
data3[['residue']] = data3[['residue']].astype(int)

# coordinates of the aligned residues
edo_coord = edo_triple.loc[ edo_triple["res_seq_num"].isin(data1[['residue']].squeeze()  ),
                            ["x_coor","y_coor","z_coor","occup"]  ]
edo_coord=edo_coord.reset_index(drop=True)
xg5A_coord = xg5A_triple.loc[ xg5A_triple["res_seq_num"].isin(data2[['residue']].squeeze()  ),
                            ["x_coor","y_coor","z_coor","occup"]  ]
xg5A_coord=xg5A_coord.reset_index(drop=True)
gy0A_coord = gy0A_triple.loc[ gy0A_triple["res_seq_num"].isin(data3[['residue']].squeeze()  ),
                            ["x_coor","y_coor","z_coor","occup"]  ]
gy0A_coord=gy0A_coord.reset_index(drop=True)

# taking the mean of the ca coordinates of the same residue when occup=0.5,
# dropping old rows and add the new mean rows
extras= list()

for i in range(0,xg5A_coord.shape[0]):
    if(xg5A_coord["occup"].iloc[i]==0.5 ) :
        if (i-1) not in extras:
            extras.append(i)

for i in range(0,(len(extras))):
    newrow=xg5A_coord.loc[[extras[i], extras[i]+1 ],["x_coor","y_coor","z_coor"]].mean(axis=0)
    xg5A_coord.loc[ extras[i]+1 ,["x_coor","y_coor","z_coor"]]=newrow

xg5A_coord=xg5A_coord.drop(extras ,axis=0)
xg5A_coord=xg5A_coord.reset_index(drop=True)

extras= list()
for i in range(0,gy0A_coord.shape[0]):
    if(gy0A_coord["occup"].iloc[i]==0.5 ) :
        if (i-1) not in extras:
            extras.append(i)
for i in range(0,(len(extras))):
    newrow=gy0A_coord.loc[[extras[i], extras[i]+1 ],["x_coor","y_coor","z_coor"]].mean(axis=0)
    gy0A_coord.loc[ extras[i]+1 ,["x_coor","y_coor","z_coor"]]=newrow

gy0A_coord=gy0A_coord.drop(extras ,axis=0)
gy0A_coord=gy0A_coord.reset_index(drop=True)

models= [edo_coord,xg5A_coord,gy0A_coord]

rmsd_mult(models)
rmsd_mult2(models)