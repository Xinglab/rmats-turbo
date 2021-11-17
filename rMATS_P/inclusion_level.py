#!/usr/bin/env python
#This script calculates the inclusion level of RNA-Seq reads

import sys

#Inc.txt
ofile=open(sys.argv[2],'w');
ofile.write('IncLevel1\tIncLevel2\tIncLevelDifference\n');

def vec2float(vec):
	res=[];
	for i in vec:
		res.append(float(i));
	return(res);

def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
	psi=[];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	for i in range(len(inc)):
		if (float(inc[i])+float(skp[i]))==0:
			psi.append("NA");
		else:
			psi.append(str(round(float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length),3)));
	return(psi);

def parse_comma_separated(s):
	if s == '':
		return list()

	return s.split(',')

#Data.txt
ifile=open(sys.argv[1]);
ifile.readline();
ilines=ifile.readlines();
for i in ilines:
	element=i.strip('\n').split('\t')
	inc1=parse_comma_separated(element[0])
	skp1=parse_comma_separated(element[1])
	inc2=parse_comma_separated(element[2])
	skp2=parse_comma_separated(element[3])
	inc1=vec2float(inc1);
	skp1=vec2float(skp1);
	inc2=vec2float(inc2);
	skp2=vec2float(skp2);
	effective_inclusion_length=int(element[4]);
	effective_skipping_length=int(element[5]);
	psi1=vec2psi(inc1,skp1,effective_inclusion_length,effective_skipping_length);
	psi2=vec2psi(inc2,skp2,effective_inclusion_length,effective_skipping_length);
	inc1='';inc2='';sum1=0;sum2=0;count1=0;count2=0;
	for j in psi1:
		inc1+=j+',';
		if j!="NA":
			sum1+=float(j);
			count1+=1;
	for j in psi2:
		inc2+=j+',';
		if j!="NA":
			sum2+=float(j);
			count2+=1;
	if (count1 == 0) or (count2 == 0):
		diff="NA";
	else:
		diff=str(round(sum1/count1-sum2/count2,3));
	ofile.write(inc1[:-1]+'\t'+inc2[:-1]+'\t'+diff+'\n');
ofile.close();
