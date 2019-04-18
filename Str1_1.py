#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Question: When considering beta sheets, are consecutive residues of a
# strand closer than neighbouring residues on different strand? We will write 
# program that calculates the distance between Cbeta atoms of consecutive 
# residues in a beta strand and the Cbeta atoms of the neighboring residues 

import numpy as np
import pandas as pd
import math

def dist(x,y):
    return( math.sqrt( ((x[0]-y[0])**2)+ ((x[1]-y[1])**2)+ ((x[2]-y[2])**2))    )

# Using the protein 1bjt from Protein Data Bank, parsing the file

with open('1bjt.pdb') as pdbfile:
    lines= list()
    sheet = list()
    for line in pdbfile:
        if line[:4] == 'ATOM' :
            splitted_line = [line[:7], line[7:11], line[11:16],line[16:20], 
                             line[20:22], line[22:26], line[26:38],line[38:46],
                             line[46:54], line[54:60], line[60:66],line[66:]]
            splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
            lines.append(splitted_line)
        if line[:5] == 'SHEET' :
            splitted_line2 = [line[:6], line[7:10], line[11:14],line[14:16],
                              line[17:20], line[21], line[22:26],line[26], 
                              line[28:31], line[32], line[33:37],line[37], 
                              line[38:40],line[41:45], line[45:48], line[49],
                              line[50:54], line[54], line[56:60],line[60:63], 
                              line[64], line[65:69],line[69]]
            splitted_line2 =[ y.replace(' ', '')  for y in splitted_line2 ]
            sheet.append(splitted_line2)

# Protein data
data = pd.DataFrame(lines, columns=["atom","atom_num","atom_name","res_name",
                                    "chain","res_seq_num","x_coor","y_coor",
                                    "z_coor","occup","temp","elt_sym"])
data[['x_coor', 'y_coor','z_coor','occup', 'temp']] = data[['x_coor', 'y_coor',
    'z_coor','occup', 'temp']].astype(float)
data[['atom_num', 'res_seq_num']] = data[['atom_num', 'res_seq_num']].astype(int)

# Sheet data
sheet_data = pd.DataFrame(sheet, columns=["sheet","str_num","sheet_id",
                                          "num_str","ini_res_name",
                                          "init_chain_id","init_seq_num",
                                          "init_inst","ter_res_name",
                                          "ter_chain_id","ter_seq_num",
                                          "ter_ins","strand_sense","curr_atom",
                                          "curr_res","curr_chain","curr_res_seq",
                                          "curr_ins","prev_atom","prev_res",
                                          "pre_chain","prev_res_seq","prev_ins"])
sheet_data[['str_num','num_str','init_seq_num',
            'ter_seq_num','strand_sense']] = sheet_data[['str_num','num_str',
                                         'init_seq_num','ter_seq_num',
                                         'strand_sense']].astype('int64' )
sheet_data['prev_res_seq'] = pd.to_numeric(sheet_data['prev_res_seq'], errors='coerce')
sheet_data['curr_res_seq'] = pd.to_numeric(sheet_data['curr_res_seq'], errors='coerce')

# Slicing information we need, since first residue in the strand doesnt have previous
#residue information, we discard them 

u = sheet_data.loc[sheet_data["str_num"]>1, ["curr_res_seq","prev_res_seq",
                   "sheet_id","strand_sense"] ]
u= u.reset_index()   

#Distances with consecutive residues
cons_dist = np.empty(u.shape[0])

for x in range(0, u.shape[0]):
    # Slicing x,y,z coordinates of the current residues
    coord = (data.loc[ ( (data["res_seq_num"]==u["curr_res_seq"][x]) & (data["atom_name"]=="CB")),
                      ["res_seq_num","x_coor","y_coor","z_coor"]  ] )
    # Slicing x,y,z coordinates of the consecutive residues
    coord2 = (data.loc[ ( (data["res_seq_num"]==u["curr_res_seq"][x]+1) & (data["atom_name"]=="CB")),
                       ["res_seq_num","x_coor","y_coor","z_coor"]  ] )
    #If Cbeta atoms exists for both residues
    if(coord.shape == (1,4) and coord2.shape == (1,4) ) : 
        point1 = coord.iloc[0][1:4]
        point2 = coord2.iloc[0][1:4]
        cons_dist[x] = dist(point1,point2)
    else:
        cons_dist[x] = None
        
#Distances with neighboring residues               
neighbor_dist = np.empty(u.shape[0])

for x in range(0, u.shape[0]):
    coord = (data.loc[ ( (data["res_seq_num"]==u["curr_res_seq"][x]) &  (data["atom_name"]=="CB")),
                      ["res_seq_num","x_coor","y_coor","z_coor"]  ] )
    # Slicing x,y,z coordinates of the neighbouring residues on a different strand
    # corresponding neighbor information is being kept on prev_res_seq vrb. in sheetdata 
    coord2 = (data.loc[ ( (data["res_seq_num"]==u["prev_res_seq"][x]) &  (data["atom_name"]=="CB")),
                       ["res_seq_num","x_coor","y_coor","z_coor"]  ] )
    if(coord.shape == (1,4) and coord2.shape == (1,4) ) : 
        point1 = coord.iloc[0][1:4]
        point2 = coord2.iloc[0][1:4]
        neighbor_dist[x] = dist(point1,point2)
    else:
        neighbor_dist[x] = None

u =    u.join(pd.DataFrame(cons_dist,columns=['cons_dist']))
u =    u.join(pd.DataFrame(np.array(u.loc[:,"curr_res_seq"]+1),columns=['consequtive_res']))
result =    u.join(pd.DataFrame(neighbor_dist,columns=['neighbor_dist']))
result = result[['curr_res_seq','consequtive_res','prev_res_seq','sheet_id','strand_sense',
                 'cons_dist','neighbor_dist']]
# There are NA values for the distances if there is no Cbeta atoms for the 
# corresponding residue. As can be seen in the last three columns of the "result"
# data frame, for parallel beta sheets, consecutive residues of a
# strand closer than neighbouring residues on different strand. However, for 
# antiparallel beta sheets, it is the opposite.

# Note: Keep in mind that these statistics are only for one PDB entry. One has to
# be careful while generalizing.
