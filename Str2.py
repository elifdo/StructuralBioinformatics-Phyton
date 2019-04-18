# -*- coding: utf-8 -*-

# Question: Find the contact matrix of a protein. For two residues to be contacting 
# distance  between Cα atoms should be less than 7 Å. 
# Find the relation between coordination number and B-factor.
# Find the relation between coordination number and Relative Surface Area.
# Is there any difference between core and surface in terms of the physicochemical
# properties of the aminoacids (small, large, hydrophobic, polar, charged etc.)

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

# Function calculating distance of two points given their x,y,z coordinates
def dist(x,y):
    return( math.sqrt( ((x[0]-y[0])**2)+ ((x[1]-y[1])**2)+ ((x[2]-y[2])**2))  )

# Reading data from PDB id 6aa3

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
data[['atom_num', 'res_seq_num']]=data[['atom_num', 'res_seq_num']].astype(int)

# Slicing Ca atoms od the protein
ca_data= data.loc[ ( (data["atom_name"]=="CA")),["res_seq_num","x_coor","y_coor","z_coor"]  ]

dist_mat = pd.DataFrame( index= range(min(ca_data["res_seq_num"]),max(ca_data["res_seq_num"])+1 ))

# Filling distance values of the residues in a matrix
for i in  ca_data["res_seq_num"]:
    for j in ca_data["res_seq_num"]:
        coord1= ca_data.loc[ ca_data["res_seq_num"]==i  ,["x_coor","y_coor","z_coor"  ]]
        coord2= ca_data.loc[ ca_data["res_seq_num"]==j  ,["x_coor","y_coor","z_coor"  ]]
        dist_mat.loc[i,j]=dist(coord1.iloc[0][0:3],coord2.iloc[0][0:3]  )

contact_mat = pd.DataFrame( index= range(min(ca_data["res_seq_num"]),max(ca_data["res_seq_num"])+1 ))

# Function converts the distance into a binary variable, two residues are
#in contact if distance is less than 7 or the other residue is not itself.
def convert(x):
    if(x<7 and x!=0) :return(1)
    else :return(0)

#Turning distance matrix into contact matrix
contact_mat= dist_mat.applymap(lambda x: convert(x))

# Turning consecutive residues from the matrix from contacting to 
# non-contacting to be able to see sheet structures more clearly.
num= ca_data["res_seq_num"].tolist()

for i in num:
    if((i-1) not in  num): contact_mat.loc[i,i+1]=0
    elif((i+1) not in  num): contact_mat.loc[i-1,i]=0
    else:
        contact_mat.loc[i,i-1]=0
        contact_mat.loc[i,i+1]=0

# Contact map: lines parallel to the diagonal and very near are helices,
# lines parallel to the diagonal and far are parallel strands
# the lines which are perpendicular to the diagonal are antiparallel strands

dotplot=plt.imshow(np.array( contact_mat  ))
plt.savefig('contact.png')

# Plot showing the relation between the coordination number and 
# temperature(B-factor) 

coord_num = pd.DataFrame([0 for y in range(contact_mat.shape[0]) ], 
                          index= ca_data["res_seq_num"] )

for i in num: 
    coord_num.loc[i]= sum(contact_mat.loc[ i,: ] )

temp = data.loc[ ( (data["atom_name"]=="CA") ),["res_seq_num","temp"]  ]

coordandtemp = pd.merge(temp, coord_num, on='res_seq_num')
coordandtemp=  coordandtemp.rename(columns={0: "coord_num"})

x= np.array(coordandtemp.loc[:,"temp"])
y= np.array(coordandtemp.loc[:,"coord_num"])
plt.scatter(x, y,  alpha=0.5)
plt.title('B-factor vs Coordination number')
plt.xlabel('B-factor')
plt.ylabel('Coordination number')
plt.savefig('bfactorvscoordnum.png')

# relative surface area data taken from GETAREA server as RTF file

with open('rsadata.rtf') as rsa:
    lines= list()
    for line in rsa:
        splitted_line = [line[:8], line[8:12], line[12:20],line[20:29], 
                         line[29:38], line[38:47], line[47:56],line[56:62]]
        splitted_line =[ x.replace(' ', '')  for x in splitted_line ]
        lines.append(splitted_line)
        
rsadata = pd.DataFrame(lines, columns=["res_name","res_num","total","apolar",
                                       "backbone","sidechain","ratio","in/out"])           
rsadata[['ratio']] = rsadata[['ratio']].astype(float)
rsadata[[ 'res_num']] = rsadata[['res_num']].astype(int)
rsa = rsadata[["res_name","res_num","ratio"]]     

# Plot of coordination number versus RSA value
x= np.array(coordandtemp.loc[:,"coord_num"])
y= np.array(rsa.loc[:,"ratio"]  )
plt.scatter(x, y,  alpha=0.5)
plt.title('Coordination number vs RSA')
plt.xlabel('Coordination number')
plt.ylabel('RSA')
plt.savefig('coordnumvsrsa.png')

# Setting residues with RSA value less than 0.05 as in(core region residue) and
# others as out (surface region residue)  
for i in range(0, len(rsadata.index)) :
    if (rsadata.loc[i,"ratio"]<5 ) :rsadata.loc[i,"in/out"]="in"
    else :rsadata.loc[i,"in/out"] ="out"
rsa = rsadata[["res_name","res_num","ratio","in/out"]]

rsa_in = rsadata.loc[rsadata["in/out"]=="in"]
rsa_out = rsadata.loc[rsadata["in/out"]=="out"]

aminoacids= set(["SER","ASN","GLN","ASP","GLU","ARG","LYS","HIS","THR","TYR",
                 "TRP","PHE","LEU","ILE","VAL","ALA","CYS","GLY","MET","PRO"])


# Funcion for automated ploting of the percentages of aminoacids in core and 
# surface region according to the classification as first and second group 
def comparison(first,name1, second, name2):
    in_first = rsa_in.loc[rsa_in["res_name"].isin(first)]
    in_second =rsa_in.loc[rsa_in["res_name"].isin(second)]
    out_first = rsa_out.loc[rsa_out["res_name"].isin(first)]
    out_second =rsa_out.loc[rsa_out["res_name"].isin(second)]
    in_first_per = in_first.shape[0]/(in_first.shape[0]+in_second.shape[0])
    in_second_per = in_second.shape[0]/(in_first.shape[0]+in_second.shape[0])
    out_first_per = out_first.shape[0]/(out_first.shape[0]+out_second.shape[0])
    out_second_per = out_second.shape[0]/(out_first.shape[0]+out_second.shape[0])
    n_groups = 2
    first = (in_first_per, out_first_per)
    second= (in_second_per, out_second_per)
    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.8
    rects1 = plt.bar(index, first, bar_width,
                     alpha=opacity,
                     color='b',
                     label=name1)
    rects2 = plt.bar(index + bar_width, second, bar_width,
                     alpha=opacity,
                     color='g',
                     label=name2)
    plt.ylabel('Percantages')
    plt.title('Percentages of '+ name1 +' and '+name2+' residues')
    plt.xticks(index + bar_width/2, ('Core', 'Surface'))
    plt.legend()
    plt.tight_layout()
    plt.savefig(name1+'_'+name2 +'_comp.png')
    

######### Polar Apolar Comparison
polar =set(["SER","ASN","GLN","ASP","GLU","ARG","LYS","HIS","THR","TYR","TRP"])
apolar = aminoacids.difference(polar)
comparison(polar,"polar",apolar,"apolar")

########## Hydrophobic Hdyrophilic Comparison

hydrophilic=set(["SER","ASN","GLN","ASP","GLU","ARG","PRO"])
hydrophobic  = aminoacids.difference(hydrophilic)
comparison(hydrophobic,"hydrophobic",hydrophilic,"hydrophilic")


########## Charged Notcharged Comparison

charged =set(["HIS","LYS","ASP","GLU","ARG"])
notcharged = aminoacids.difference(charged)
comparison(charged,"charged",notcharged,"notcharged")

########## Small Large Comparison

small =set(["CYS","GLY","ASP","ASN","ALA","PRO","THR","SER"])
large = aminoacids.difference(small)
comparison(small,"small",large,"large")