#!/usr/bin/env python
#This script calculate FDR for the P values

import sys,re,os,numpy

def myorder(p,reverse):
	p_dict={};
	#print('len(p)');print(len(p));
	for index in range(len(p)):
		if not(p[index] in p_dict):
			p_dict[p[index]]=[];
		p_dict[p[index]].append(index+1);
	#print('p_dict');print(p_dict);
	o_index=sorted(p_dict,reverse=reverse);o=[];
	#print('o_index');print(o_index);
	for index in o_index:
		for this_p in p_dict[index]:
			o.append(this_p);
	return(o);

def mycummin(p):
	res=[];
	if len(p) == 0:
		return res;

	current_min = p[0];
	for value in p:
		if value < current_min:
			current_min = value;

		res.append(current_min);

	return(res);

def myFDR(p):
	lp=len(p);
	i=range(lp,0,-1);
	o=myorder(p,True);
	ro=myorder(o,False);
	p_new=[];
	for index in range(len(o)):
		p_new.append(p[o[index]-1]*lp/i[index]);
	p_mycummin=mycummin(p_new);
	p_mycummin_new=[];
	for index in range(len(p_mycummin)):
		p_mycummin_new.append(min(p_mycummin[index],1));
	res=[];
	for index in range(len(p_mycummin_new)):
		res.append(p_mycummin_new[ro[index]-1]);
	return(res);

#input P values
ifile=open(sys.argv[1]);title=ifile.readline();
ilines=ifile.readlines();
P=[];
for i in ilines:
	element=re.findall('[^\t\n]+',i);
	P.append(float(element[-1]));

#Convert FDR
FDR=myFDR(P)

#output FDR
ofile=open(sys.argv[2],'w');ofile.write(title[:-1]+'\tFDR\n');
for i in range(len(ilines)):
	ofile.write(ilines[i][:-1]+'\t'+str(FDR[i])+'\n');
ofile.close();
