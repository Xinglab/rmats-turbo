
# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab')


# In[2]:

import re


# In[3]:

all_prob1 = []
with open('rMATS_Result_P_backup.txt', 'r') as inputf:
    all_lines = inputf.readlines()
    for line in all_lines[1:]:
        element=re.findall('[^ \t\n]+', line)[-1]
        all_prob1.append(float(element))


# In[4]:

all_prob1 = np.array(all_prob1)


# In[5]:

# all_prob = all_prob[all_prob > 0]
# all_prob = all_prob[all_prob < 1]


# In[6]:

all_prob1 = np.log10(all_prob1+1)


# In[7]:

all_prob1


# In[8]:

fig1 = plt.figure(dpi=100)
ax1 = fig1.add_subplot(1, 1, 1)
ax1.scatter(range(0, all_prob1.shape[0]), all_prob1, s=8, c='b', marker='o', label='first')


# In[9]:

all_prob2 = []
with open('output.txt', 'r') as inputf:
    all_lines = inputf.readlines()
    for line in all_lines[1:]:
        element=re.findall('[^ \t\n]+', line)[-1]
        all_prob2.append(float(element))


# In[10]:

all_prob2 = np.array(all_prob2)


# In[11]:

# all_prob = all_prob[all_prob > 0]
# all_prob = all_prob[all_prob < 1]


# In[12]:

all_prob2 = np.log10(all_prob2+1)


# In[13]:

all_prob2


# In[14]:

fig2 = plt.figure(dpi=100)
ax2 = fig2.add_subplot(1, 1, 1)
ax2.scatter(range(0, all_prob2.shape[0]), all_prob2, s=8, c='r', marker='o', label='second')


# In[15]:

fig2.dpi


# In[16]:

all_prob1 = []
with open('rMATS_Result_P_backup.txt', 'r') as inputf:
    all_lines = inputf.readlines()
    for line in all_lines[1:]:
        element=re.findall('[^ \t\n]+', line)[-1]
        all_prob1.append(float(element))
all_prob1 = np.array(all_prob1)
all_prob2 = []
with open('rMATS_Result_P.txt', 'r') as inputf:
    all_lines = inputf.readlines()
    for line in all_lines[1:]:
        element=re.findall('[^ \t\n]+', line)[-1]
        all_prob2.append(float(element))
all_prob2 = np.array(all_prob2)
mean((all_prob1-all_prob2)**2)


# In[ ]:



