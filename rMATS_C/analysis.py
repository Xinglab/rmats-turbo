#!/usr/bin/env python
# -*- coding: utf-8 -*-


##
# @file analysis.py
# @brief 
# @author Zhijie Xie
# @version 1.0
# @date 2015-10-28




import numpy as np
import matplotlib.pyplot as plt
import re
from scipy.stats import pearsonr
from optparse import OptionParser
from os import listdir


VERSION = 'v1.0.0'
USAGE = "usage: %prog [options] arg1 arg2"


def get_options():
    """TODO: Docstring for get_options.
    :returns: TODO

    """
    parser = OptionParser(usage=USAGE, version=VERSION)

    parser.add_option('--pyf', action='store', type='string',
                      help='Raw input file.', dest='pyf')
    parser.add_option('--cf', action='store', type='string',
                      help='Raw input file.', dest='cf')
    parser.add_option('--tn', action='store', type='string',
                      help='test name.', dest='tn')

    return parser.parse_args()


def analyze_py_c(tn, pyf, cf):
    """TODO: Docstring for analyze_py_c.
    :returns: TODO

    """
    all_prob1 = []
    with open(pyf, 'r') as inputf:
        all_lines = inputf.readlines()
        for line in all_lines[1:]:
            element=re.findall('[^ \t\n]+', line)[-1]
            all_prob1.append(float(element))
    all_prob1 = np.array(all_prob1)
    all_prob2 = []
    with open(cf, 'r') as inputf:
        all_lines = inputf.readlines()
        for line in all_lines[1:]:
            element=re.findall('[^ \t\n]+', line)[-1]
            all_prob2.append(float(element))
    all_prob2 = np.array(all_prob2)

    # pearsonr of raw p-value

    print 'pearsonr of raw p-value: %e, %e' % (pearsonr(all_prob1, all_prob2))

    # The first value in the following tuple is the Pearson’s correlation coefficient,
    # and the second one is a P-value roughly indicates the probability of an uncorrelated system
    # producing datasets that have a Pearson correlation at least as extreme as the one
    # computed from these datasets. The p-values are not entirely reliable but are probably
    # reasonable for datasets larger than 500 or so.

    # So, we can refer from the following tuple that our results is correlative with yours, because the p-value is 0.0.

    # -log10(p-value) and set log10(0) to 10.

    all_prob1[all_prob1 <= 10-10] = 10
    all_prob1[all_prob1 != 10] = -np.log10(all_prob1[all_prob1 != 10])
    all_prob2[all_prob2 <= 10-10] = 10
    all_prob2[all_prob2 != 10] = -np.log10(all_prob2[all_prob2 != 10])

    # Drawing the picture.

    fig1 = plt.figure(dpi=100)
    ax1 = fig1.add_subplot(1, 1, 1)
    ax1.set_xlabel(u'Py', fontsize=14)
    ax1.set_ylabel(u'C', fontsize=14)
    ax1.set_title(tn,fontsize=18)
    ax1.scatter(all_prob1, all_prob2, s=8, c='b', marker='o', label='first')
    ax1.plot(ax1.get_xlim(), ax1.get_ylim(), ls="--", c=".3")
    plt.savefig(tn, dpi=600)
    # plt.show()

    # pearsonr of -log10(p-value)

    print 'pearsonr of -log10(p-value): %e, %e' % (pearsonr(all_prob1, all_prob2))


def diff():
    """TODO: Docstring for diff.
    :returns: TODO

    """
    ids = []
    tmp1 = listdir('.')
    for i in tmp1:
        if i[-2:] == 'py':
            b = [j for j in listdir(i) if j[-3:] == 'txt'][0]
            with open(i+'/'+b, 'r') as inputf:
                all_lines = inputf.readlines()
                for line in all_lines[1:]:
                    element=re.findall('[^ \t\n]+', line)[0]
                    ids.append(float(element))
    ids = np.array(ids)
    all_prob1 = []
    tmp1 = listdir('.')
    for i in tmp1:
        if i[-2:] == 'py':
            b = [j for j in listdir(i) if j[-3:] == 'txt'][0]
            with open(i+'/'+b, 'r') as inputf:
                all_lines = inputf.readlines()
                for line in all_lines[1:]:
                    element=re.findall('[^ \t\n]+', line)[-1]
                    all_prob1.append(float(element))
    all_prob1 = np.array(all_prob1)
    all_prob2 = []
    tmp1 = listdir('.')
    for i in tmp1:
        if i[-1:] == 'c':
            b = [j for j in listdir(i) if j[-3:] == 'txt'][0]
            with open(i+'/'+b, 'r') as inputf:
                all_lines = inputf.readlines()
                for line in all_lines[1:]:
                    element=re.findall('[^ \t\n]+', line)[-1]
                    all_prob2.append(float(element))
    all_prob2 = np.array(all_prob2)
    content1 = []
    tmp1 = listdir('.')
    for i in tmp1:
        if i[-2:] == 'py':
            b = [j for j in listdir(i) if j[-3:] == 'txt'][0]
            with open(i+'/'+b, 'r') as inputf:
                all_lines = inputf.readlines()
                for line in all_lines[1:]:
                    element=re.findall('[^ \t\n]+', line)
                    content1.append(element)
    content2 = []
    tmp1 = listdir('.')
    for i in tmp1:
        if i[-1:] == 'c':
            b = [j for j in listdir(i) if j[-3:] == 'txt'][0]
            with open(i+'/'+b, 'r') as inputf:
                all_lines = inputf.readlines()
                for line in all_lines[1:]:
                    element=re.findall('[^ \t\n]+', line)
                    content2.append(element)
    a = all_prob1 - all_prob2
    idx = np.where(a>0.01)
    with open('large_margin.txt', 'w') as tmp:
        for i in idx[0]:
            tmp.write(str(content1[int(i)]))
            tmp.write(str(content2[int(i)]))
            tmp.write('\n')
        idx = np.where(a>0.001)
    with open('special.txt', 'w') as tmp:
        for i in idx[0]:
            tmp.write('%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % (int(content1[i][0]), str(content1[i][1]), str(content1[i][2]), str(content1[i][3]), str(content1[i][4]), str(content1[i][5]), str(content1[i][6])))


def duplicate():
    """TODO: Docstring for duplicate.
    :returns: TODO

    """
    content = []
    with open('/u/home/s/shiehshi/test_output/cluster_testing/SE.JC.input.txt', 'r') as inputf:
        all_lines = inputf.readlines()
        for line in all_lines[1:]:
            element=re.findall('[^ \t\n]+', line)
            element[1] = ','.join([element[1] for i in range(33)])
            element[2] = ','.join([element[2] for i in range(33)])
            element[3] = ','.join([element[3] for i in range(33)])
            element[4] = ','.join([element[4] for i in range(33)])
            # element[1] = element[1]
            # element[2] = element[2]
            # element[3] = element[3]
            # element[4] = element[4]
            content.append(element)

    content = ['%d\t%s\t%s\t%s\t%s\t%s\t%s\n' % (int(line[0]), str(line[1]), str(line[2]), str(line[3]), str(line[4]), str(line[5]), str(line[6])) for line in content]
    with open('dummy.txt', 'w') as tmp:
        tmp.write(all_lines[0])
        for i in range(5000):
            tmp.writelines(content)



def compare():
    """TODO: Docstring for compare.
    :returns: TODO

    """
    all_prob1 = []
    with open('qsub_run_1_1_0001_py/PC3E_py_1_1_0001_output.txt', 'r') as inputf:
        all_lines = inputf.readlines()
        for line in all_lines[1:]:
            element=re.findall('[^ \t\n]+', line)[-1]
            all_prob1.append(float(element))
    all_prob1 = np.array(all_prob1)
    all_prob2 = []
    with open('qsub_run_1_1_0001_c/PC3E_c_1_1_0001_output.txt', 'r') as inputf:
        all_lines = inputf.readlines()
        for line in all_lines[1:]:
            element=re.findall('[^ \t\n]+', line)[-1]
            all_prob2.append(float(element))
    all_prob2 = np.array(all_prob2)


def main():
    """TODO: Docstring for main.
    :returns: TODO

    """
    (options, args) = get_options()
    analyze_py_c(options.tn, options.pyf, options.cf)


if __name__ == "__main__":
    main()
