import re,os,sys,warnings,numpy,scipy,math,itertools;

from scipy import stats;
from numpy import *;
from multiprocessing import Pool;
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log;

numpy.random.seed(1231);
warnings.filterwarnings('ignore');

#ReadLength
read_length=int(sys.argv[3]);

#JunctionLength
junction_length=int(sys.argv[4]);

#MultiProcessor
MultiProcessor=1;
if len(sys.argv)>=6:
	MultiProcessor=int(sys.argv[5]);

#splicing difference cutoff
cutoff=0.1;
if len(sys.argv)>=7:
	cutoff=float(sys.argv[6]);


#binomial MLE optimization functions
def logit(x):
	if x<0.01:
		x=0.01;
	if x>0.99:
		x=0.99;
	return(log(x/(1-x)));

def myfunc_multivar(x,*args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	#print('psi1');print(psi1);print('psi2');print(psi2);
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=pow(logit(i)-logit(x[0]),2);
	sum1=sum1/var1/2;
	for i in psi2:
		sum2+=pow(logit(i)-logit(x[1]),2);
	sum2=sum2/var2/2;
	return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x[0]),2)+pow(stats.norm.ppf(x[1]),2)-2*rho*stats.norm.ppf(x[0])*stats.norm.ppf(x[1])));

def myfunc_multivar_der(x,*args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=-2*(logit(i)-logit(x[0]))/x[0]/(1-x[0]);
	sum1=sum1/var1/2;
	res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[0])-2*rho*stats.norm.ppf(x[1]))/stats.norm.pdf(stats.norm.ppf(x[0]));
	#print('1');print(x[1]);print(res1);print(stats.norm.pdf(stats.norm.ppf(x[0])));
	for i in psi2:
		sum2+=-2*(logit(i)-logit(x[1]))/x[1]/(1-x[1]);
	sum2=sum2/var2/2;
	res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[1])-2*rho*stats.norm.ppf(x[0]))/stats.norm.pdf(stats.norm.ppf(x[1]));
	return(numpy.array([res1,res2]));

def myfunc_1(x, *args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=pow(logit(i)-logit(x+cutoff),2);
	sum1=sum1/var1/2;
	for i in psi2:
		sum2+=pow(logit(i)-logit(x),2);
	sum2=sum2/var2/2;
	return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)));

def myfunc_der_1(x, *args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff);
	sum1=sum1/var1/2;
	res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff));
	for i in psi2:
		sum2+=-2*(logit(i)-logit(x))/x/(1-x);
	sum2=sum2/var2/2;
	res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x));
	return(res1+res2);

def myfunc_2(x, *args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=pow(logit(i)-logit(x),2);
	sum1=sum1/var1/2;
	for i in psi2:
		sum2+=pow(logit(i)-logit(x+cutoff),2);
	sum2=sum2/var2/2;
	return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)));

def myfunc_der_2(x, *args):
	psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
	sum1=0;sum2=0;
	for i in psi1:
		sum1+=-2*(logit(i)-logit(x))/(x)/(1-x);
	sum1=sum1/var1/2;
	res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff));
	for i in psi2:
		sum2+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff);
	sum2=sum2/var2/2;
	res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x));
	return(res1+res2);

def myfunc_marginal_2der(x, args):
	I=args[0];S=args[1];beta=args[2];var=args[3];
	inclusion_length=args[4];
	skipping_length=args[5]
	#print('test x');print(x);print(var);
	temp1=1/pow(x,2)/pow(1-x,2)*(((2*x-1)*(logit(beta)-logit(x))-1)/var-1);#temp1=(2*x-1)/pow(x,2)/pow(1-x,2)*((logit(beta)-logit(x)-1/(2*x-1))/var+1);
	temp2=I*skipping_length*((2*inclusion_length+skipping_length)*x+skipping_length*(1-x))/pow(x,2)/pow(inclusion_length*x+skipping_length*(1-x),2);
	temp3=S*inclusion_length*((inclusion_length+2*skipping_length)*(1-x)+inclusion_length*x)/pow(1-x,2)/pow(inclusion_length*x+skipping_length*(1-x),2);
	#print('test');print(beta);print(x);print(var);print(temp1);print(temp2);print(temp3);
	return(temp1-temp2-temp3);
	
def myfunc_marginal(x, *args):
	beta=x;
	I=args[0];S=args[1];psi=args[2];var=args[3];
	inclusion_length=args[4];skipping_length=args[5];
	sum=0;
	for i in range(len(psi)):
		new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]));
		f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var-log(psi[i])-log(1-psi[i]);
		f1_2der=abs(myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length]));
		sum+=(0.5*log(abs(f1_2der)+0.00001)-f1);
	return(sum);

def myfunc_marginal_der(x, *args):
	beta=x;
	I=args[0];S=args[1];psi=args[2];var=args[3];
	inclusion_length=args[4];skipping_length=args[5];
	sum=0;
	for i in range(len(psi)):
		new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]));
		f1_3der=1*(2*psi[i]-1)/pow(psi[i],2)/pow(1-psi[i],2)/beta/(1-beta)/var;
		f1_2der=myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length]);
		#print('test2');print(f1_2der);
		f1_der=(logit(psi[i])-logit(beta))/beta/(1-beta)/var;
		f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var+log(psi[i])+log(1-psi[i]);
		sum+=(0.5*f1_3der/(f1_2der)-f1_der);
	return(sum);

def myfunc_marginal_1(x, *args):
	beta2=x;beta1=x+cutoff;
	I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
	I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
	inclusion_length=args[8];skipping_length=args[9];
	return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_1_der(x, *args):
	beta2=x;beta1=x+cutoff;
	I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
	I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
	inclusion_length=args[8];skipping_length=args[9];
	return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_2(x, *args):
	beta2=x+cutoff;beta1=x;
	I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
	I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
	inclusion_length=args[8];skipping_length=args[9];
	return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_2_der(x, *args):
	beta2=x+cutoff;beta1=x;
	I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
	I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
	inclusion_length=args[8];skipping_length=args[9];
	return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));
	
def myfunc_individual(x,*args):
	I=args[0];S=args[1];beta=args[2];var=args[3];
	inclusion_length=args[4];
	skipping_length=args[5]
	new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
	return(-1*(I*log(new_psi)+S*log(1-new_psi)-(logit(x)-logit(beta))*(logit(x)-logit(beta))/2/var-log(x)-log(1-x)));

def myfunc_individual_der(x,*args):
	I=args[0];S=args[1];beta=args[2];var=args[3];
	inclusion_length=args[4];
	skipping_length=args[5];
	new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
	new_psi_der=inclusion_length*skipping_length/pow(inclusion_length*x+skipping_length*(1-x),2);
	return(-1*(I/new_psi*new_psi_der-S/(1-new_psi)*new_psi_der-(logit(x)-logit(beta))/var/x/(1-x)-1/x+1/(1-x) ));

def myfunc_likelihood(x, args):
	I=args[0];S=args[1];beta=args[2];var=args[3];sum=0;N=I+S;
	#return(-1*(-log(sqrt((I+S)*x*(1-x)))-(I-(I+S)*x)*(I-(I+S)*x)/2/((I+S)*x*(1-x))-log(sqrt(var))-(x-beta)*(x-beta)/2/var));
	#print('debug');print(N);print(var);print(x);print(beta);
	if N==0:
		return(0);
	return(-0.5*((I-N*x)*(I-N*x)/(N*x)+(S-N*(1-x))*(S-N*(1-x))/(N*(1-x)))-log(sqrt(var))-(x-beta)*(x-beta)/2/var);

def MLE_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	beta_0=sum(psi1)/len(psi1);
	beta_1=sum(psi2)/len(psi2);
	var1=10*scipy.var(numpy.array(psi1)-beta_0);
	var2=10*scipy.var(numpy.array(psi2)-beta_1);
	if var1<=0.01:
		var1=0.01;
	if var2<=0.01:
		var2=0.01;
	print('var1');print(var1);print('var2');print(var2);
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		print('var1');print(var1);print('var2');print(var2);
		#if abs(sum(psi1)/len(psi1)-sum(psi2)/len(psi2))>cutoff:
		if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
		else:#minize psi1 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
		print('constrain_1xopt');print('theta');print(theta1);print(theta2);print(xopt);
		#else:
		#	theta1=sum(psi1)/len(psi1);theta2=sum(psi2)/len(psi2);
		beta_0=theta1;beta_1=theta2;
		#iteration of psi
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		print('constrain_2xopt');
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
		print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0]);
		psi1=new_psi1;psi2=new_psi2;
		print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum);
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	#print('constrain');print(theta1);print(theta2);print(psi1);print(psi2);print(current_sum);print(likelihood_sum);
	#print(xopt);
	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
	#return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	beta_0=sum(psi1)/len(psi1);
	beta_1=sum(psi2)/len(psi2);
	var1=10*scipy.var(numpy.array(psi1)-beta_0);
	var2=10*scipy.var(numpy.array(psi2)-beta_1);
	if var1<=0.01:
		var1=0.01;
	if var2<=0.01:
		var2=0.01;
	print('var1');print(var1);print('var2');print(var2);
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
		beta_0=float(xopt[0][0]);
		beta_1=float(xopt[0][1]);
		print('unconstrain_1xopt');print(xopt);
		print('theta');print(beta_0);print(beta_1);print('theta_end');
		#iteration of psi
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
		print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0]);
		psi1=new_psi1;psi2=new_psi2;print
		print('count');print(count);('previous_sum');print(previous_sum);print('current_sum');print(current_sum);
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	#print('unconstrain');print(beta_0);print(beta_0+beta_1);print(psi1);print(psi2);print(current_sum);print(likelihood_sum);
	#print(xopt);
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);
	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
	#return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	#initial value
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	beta_0=sum(psi1)/len(psi1);
	beta_1=sum(psi2)/len(psi2);
	var1=10*scipy.var(numpy.array(psi1)-beta_0);
	var2=10*scipy.var(numpy.array(psi2)-beta_1);
	if var1<=0.01:
		var1=0.01;
	if var2<=0.01:
		var2=0.01;
	print('var1');print(var1);print('var2');print(var2);
	
	#MLE of the full likelihood
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
		beta_0=float(xopt[0][0]);
		beta_1=float(xopt[0][1]);
		print('unconstrain_MLE_xopt');print(xopt);
		print('theta');print(beta_0);print(beta_1);print('theta_end');
		#iteration of psi
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);
			#likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);
			#likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
		print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0]);
		psi1=new_psi1;psi2=new_psi2;
		print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum);
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	#print('unconstrain');print(beta_0);print(beta_0+beta_1);print(psi1);print(psi2);print(current_sum);print(likelihood_sum);
	print('unconstrain_psi_MLE');print(psi1);print(psi2);
	print('unconstrain_beta_MLE');print(beta_0);print(beta_1);
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);
	
	#MLE of the marginal likelihood
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		print('unconstrain_MLE_marginal_value_1_der');print(myfunc_marginal_der(beta_0,i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length));
		print('unconstrain_MLE_marginal_value_2_der');print(myfunc_marginal_der(beta_1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length))
		xopt1=fmin_l_bfgs_b(myfunc_marginal,[beta_0],myfunc_marginal_der,args=[i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
		xopt2=fmin_l_bfgs_b(myfunc_marginal,[beta_1],myfunc_marginal_der,args=[i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
		beta_0=float(xopt1[0][0]);
		beta_1=float(xopt2[0][0]);
		print('unconstrain_1xopt');print(xopt1);
		print('unconstrain_2xopt');print(xopt2);
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);
		psi1=new_psi1;psi2=new_psi2;
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);
	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
	#initial value
	psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	beta_0=sum(psi1)/len(psi1);
	beta_1=sum(psi2)/len(psi2);
	var1=10*scipy.var(numpy.array(psi1)-beta_0);
	var2=10*scipy.var(numpy.array(psi2)-beta_1);
	if var1<=0.01:
		var1=0.01;
	if var2<=0.01:
		var2=0.01;
	print('var1');print(var1);print('var2');print(var2);
	
	#MLE of the full likelihood
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		#iteration of beta
		beta_0=sum(psi1)/len(psi1);
		beta_1=sum(psi2)/len(psi2);
		print('var1');print(var1);print('var2');print(var2);
		#if abs(sum(psi1)/len(psi1)-sum(psi2)/len(psi2))>cutoff:
		if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
		else:#minize psi1 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
		print('constrain_1xopt');print('theta');print(theta1);print(theta2);print(xopt);
		#else:
		#	theta1=sum(psi1)/len(psi1);theta2=sum(psi2)/len(psi2);
		beta_0=theta1;beta_1=theta2;
		#iteration of psi
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		print('constrain_2xopt');
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
			#likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
		print('new_psi[0]');print(new_psi1[0]);print(new_psi2[0]);
		psi1=new_psi1;psi2=new_psi2;
		print('count');print(count);print('previous_sum');print(previous_sum);print('current_sum');print(current_sum);
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);
	#print('constrain');print(theta1);print(theta2);print(psi1);print(psi2);print(current_sum);print(likelihood_sum);
	#print(xopt);
	
	#MLE of the marginal likelihood
	iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
	while((iter_cutoff>0.01)&(count<=iter_maxrun)):
		count+=1;
		if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_marginal_1,[beta_1],myfunc_marginal_1_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			beta_1 = max(min(float(xopt[0]),1-cutoff),0);beta_0=beta_1+cutoff;
		else:#minize psi1 if this is the case
			xopt = fmin_l_bfgs_b(myfunc_marginal_2,[beta_0],myfunc_marginal_2_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
			beta_0 = max(min(float(xopt[0]),1-cutoff),0);beta_1=beta_0+cutoff;
		print('constrain_xopt');print(xopt);
		new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
		for i in range(len(psi1)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
		for i in range(len(psi2)):
			xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
			new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);print(xopt);
		psi1=new_psi1;psi2=new_psi2;
		if count>1:
			iter_cutoff=abs(previous_sum-current_sum);
		previous_sum=current_sum;
	if count>iter_maxrun:
		return([current_sum,[psi1,psi2,0,0,var1,var2]]);

	return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
	#return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

	
#Random Sampling Function
def likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id):
	print('testing'+str(id));
	if flag==0:
		return([1,1]);
	else:
		res=MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
		if abs(res[1][2]-res[1][3])<=cutoff:
			print('1<=cutoff');print(res);print((res[1][2]-res[1][3]));
			return([1,res[0]]);
		else:
			res_constrain=MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
			print('2>cutoff');print(res);print(res_constrain);
			print(abs(res_constrain[0]-res[0]));print('2end');
			#print(abs(res_constrain[2]-res[2]));print('2end_marginal');			
			if len(i1)<=3:
				return([1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1)]);
			else:
				return([1-scipy.stats.chi2.cdf(2*(abs(res_constrain[0]-res[0])),1)]);

#MultiProcessorFunction
def MultiProcessorPool(n_original_diff):
	i1=n_original_diff[0];i2=n_original_diff[1];s1=n_original_diff[2];s2=n_original_diff[3];
	effective_inclusion_length=n_original_diff[4];effective_skipping_length=n_original_diff[5];
	flag=n_original_diff[6];id=n_original_diff[7];
	P=likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id);
	return(P);

#Function for vector handling
def vec2float(vec):
	res=[];
	for i in vec:
		res.append(float(i));
	return(res);

#add 1 in both inclusion and skipping counts for robustness in small sample size
def vecAddOne(vec):
	res=[];
	for i in vec:
		res.append(i);
	return(res);

def vecprod(vec):
	res=1;
	for i in vec:
		res=res*i;
	return(res);

def vecadd(vec1,vec2):
	res=[];
	for i in range(len(vec1)):
		res.append(vec1[i]+vec2[i]);
	return(res);
	
def vec2remove0psi(inc,skp):
	res1=[];res2=[];
	for i in range(len(inc)):
		if (inc[i]!=0) | (skp[i]!=0):
			res1.append(inc[i]);res2.append(skp[i]);
	return([res1,res2]);

def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
	psi=[];
	inclusion_length=effective_inclusion_length;
	skipping_length=effective_skipping_length;
	for i in range(len(inc)):
		if (float(inc[i])+float(skp[i]))==0:
			psi.append(0.5);
		else:
			psi.append(float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length));
	return(psi);

def vec210(vec):
	res=[];
	for i in vec:
		if i>0:
			res.append(1);
		else:
			res.append(-1);
	return(res);

def myttest(vec1,vec2):
	if (len(vec1)==1) & (len(vec2)==1):
		res=stats.ttest_ind([vec1[0],vec1[0]],[vec2[0],vec2[0]]);
	else:
		res=stats.ttest_ind(vec1,vec2);
	return(res);

ifile=open(sys.argv[1]);
title=ifile.readline();
#analyze the title of the inputed data file to find the information of how much simulation are involved
#the min simulated round is 10, each time it increases by 10 times
element=re.findall('[^ \t\n]+',title);
ofile=open(sys.argv[2]+'/rMATS_Result_P.txt','w');
ofile.write(title[:-1]+'\tPValue'+'\n');

list_n_original_diff=[];probability=[];psi_list_1=[];psi_list_2=[];
ilines=ifile.readlines();
for i in range(len(ilines)):
	if 'NA' in ilines[i]:
		list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,element[0]]);
		continue;
	element=re.findall('[^ \t\n]+',ilines[i]);
	inc1=re.findall('[^,]+',element[1]);skp1=re.findall('[^,]+',element[2]);inc2=re.findall('[^,]+',element[3]);skp2=re.findall('[^,]+',element[4]);
	effective_inclusion_length=int(element[5]);
	effective_skipping_length=int(element[6]);
	inc1=vec2float(inc1);skp1=vec2float(skp1);inc2=vec2float(inc2);skp2=vec2float(skp2);
	if (vecprod(vecadd(inc1,skp1))==0) | (vecprod(vecadd(inc2,skp2))==0):
		list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,0,element[0]]);
	else:
		inc1=vecAddOne(inc1);skp1=vecAddOne(skp1);inc2=vecAddOne(inc2);skp2=vecAddOne(skp2);
		psi_list_1.append(sum(inc1)/(sum(inc1)+sum(skp1)));
		psi_list_2.append(sum(inc2)/(sum(inc2)+sum(skp2)));
		#temp1=vec2remove0psi(inc1,skp1);temp2=vec2remove0psi(inc2,skp2);
		#inc1=temp1[0];skp1=temp1[1];inc2=temp2[0];skp2=temp2[1];
		list_n_original_diff.append([inc1,inc2,skp1,skp2,effective_inclusion_length,effective_skipping_length,1,element[0]]);
	#if i>2:
	#	break;

rho=0.9
rho=stats.pearsonr(numpy.array(psi_list_1),numpy.array(psi_list_2));rho=rho[0];
if rho>0.9:
	rho=0.9;
rho=0.9
print('rho');print(rho);

if MultiProcessor>1:
	pool=Pool(processes=MultiProcessor);
	probability=pool.map(MultiProcessorPool,list_n_original_diff);
else:
	for i in range(len(list_n_original_diff)):
		probability.append(MultiProcessorPool(list_n_original_diff[i]));
#print(probability);
index=0;
for i in range(len(ilines)):
    element=re.findall('[^ \t\n]+',ilines[i]);
    print(probability);
    ofile.write(ilines[i][:-1]+'\t'+str(probability[i][0])+'\n');
ofile.close();
