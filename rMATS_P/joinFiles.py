#
## this program joins two files based on the column
#

### import necessary libraries
import re,os,sys,logging,time,datetime;

### checking out the number of arguments
if (len(sys.argv)<6): 
  print('Not enough arguments!!');
  print ('It takes at least 5 arguments.');
  print ('Usage:\n\tProgramName inFile_1 inFile_2 col_1 col_2 outFile');
  print ('Example\n\tProgramName allEvents exisingEvents 0 0 out.joined.txt');
  sys.exit();

def listToString(x):
  rVal = '';
  for a in x:
    rVal += a+' ';
  return rVal;


### setting up the logging format 
#logging.basicConfig(level=logging.DEBUG,
#                    format='%(asctime)s %(message)s',
#                    filename='log.joinFiles.'+  str(datetime.datetime.now()),
#                    filemode='w')

##### Getting Start Time ######
#logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

###
iFile_1 = open(sys.argv[1]); ## input file 1
iFile_2 = open(sys.argv[2]); ## input file 2
col_1 = int(sys.argv[3]); ## join on this column on file 1
col_2 = int(sys.argv[4]); ## join on this column on file 2
oFile = open(sys.argv[5], 'w'); ## out fastq file

fd={}; 

c=0;
dupKey = 0;

header=iFile_1.readline().strip();
oFile.write(header+'\t');

for line in iFile_1: ## for each line
  ele = line.strip().split('\t');
  key = ele[col_1];
  if key in fd: ## should not happen
    dupKey += 1;
#    logging.debug("duplicate key: %s" % line.strip());
  else:
    fd[key] = line.strip();
#logging.debug("Finished populating first input file dictionary with %d keys and %d duplicates" % (len(fd), dupKey));

header = iFile_2.readline().strip();
oFile.write(header + '\n');

numFound =0;
notFound =0;
for line in iFile_2: ## for each line
  ele = line.strip().split('\t');
  key = ele[col_2];
  if key in fd: ## let's write it out
    numFound += 1;
    oFile.write(fd[key] + '\t' + line.strip() + '\n');
  else: ## not found
    notFound += 1;
#logging.debug("Finished accesing 2nd file. Found %d lines and did not find %d lines" % (numFound, notFound));

iFile_1.close();
iFile_2.close();
oFile.close();

#############
## calculate total running time
#############
#logging.debug("Program ended");
currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
#logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));

sys.exit(0);
