
m sys import argv
import getopt
import re
import numpy as np

proj_path = '/u/biancad/COS424'
data_path = 'data'
train_file = 'intersected_final_chr1_cutoff_20_train_revised.bed'
sample_file = 'intersected_final_chr1_cutoff_20_sample.bed'
test_file = 'intersected_final_chr1_cutoff_20_test.bed'
num_interval = '5'
interval_size = '100'

try:
  opts, remains = getopt.getopt(argv[1:],"n:i:f:",["--num=","--intsize=","--file="])
except getopt.GetoptError:
  print 'error'
  sys.exit(2)
  
for opt, arg in opts:
  if opt in ("-n", "--num"):
num_interval = arg
  elif opt in ("-f", "--file"):
train_file = arg
  elif opt in ("-i","--intsize"):
    interval_size = arg
# first, collect the features
print 'reading train data'

inputfile = proj_path + '/' + data_path + '/' + train_file
raw_data = open(inputfile)

data = {}
label = {}
for line in raw_data:
  row = line.split()
  start_loc = int(row[1])
  data[start_loc] = np.loadtxt(row[4:-1])
  label[start_loc] = int(row[-1])
raw_data.close()

print 'buidling features'

features = {}
for i in data:
  features[i] = data[i]
  dat_size = data[i].shape
  
  for j in range(1,(int(num_interval)-1)/2+1):
    start_loc = i + (j-1)*int(interval_size)
    # pos offset
    end_int = i + j*int(interval_size)
    sum = np.zeros(dat_size)
    count = np.zeros(dat_size)
    k = start_loc + 2
    while(k <= end_int):
 if k in data:
   checknan = ~np.array(np.isnan(data[k])) #skip nan, so need to keep count correctly
   data[k][np.isnan(data[k])] = 0 #change nan to zero to skip it in the average
   sum += data[k]
   count += checknan
 k += 2
    average = np.zeros(dat_size)
    count[count == 0] = 1 #clear the case where count is zero
    average = sum/count
    features[i] = np.concatenate((features[i],average),axis=0)
    # neg offset
    end_int = i - j*int(interval_size)
    sum = np.zeros(dat_size)
    count = np.zeros(dat_size)
    start_loc = i - (j-1)*int(interval_size)
    k = start_loc - 2
  
    while(k >= end_int):
 if k in data:
   checknan = ~np.array(np.isnan(data[k])) #skip nan, so need to keep count correctly
   data[k][np.isnan(data[k])] = 0 #change nan to zero to skip it in the average
   sum += data[k]
   count += checknan
 k -= 2
    average = np.zeros(dat_size)
    count[count == 0] = 1 #clear the case where count is zero 
    average = sum/count
    features[i] = np.concatenate((features[i],average),axis=0)
# now, collect the teacher
print 'reading sample/test data'

inputfile = proj_path + '/' + data_path + '/' + sample_file
raw_data = open(inputfile)

teacher = {}
for line in raw_data:
  row = line.split()
  if(int(row[-1]) == 1):
    loc = int(row[1])
    teacher[loc] = float(row[4])
raw_data.close()
inputfile = proj_path + '/' + data_path + '/' + test_file
raw_data = open(inputfile)
for line in raw_data:
  row = line.split()
  if(int(row[-1]) == 0):
    loc = int(row[1])
    teacher[loc] = float(row[4])
raw_data.close()

# write to files 
print 'writing to files'

output_Xtrain = './cos424hw2_Xtrain_' + num_interval
output_Xtest = './cos424hw2_Xtest_' + num_interval
output_Ytrain = './cos424hw2_Ytrain_' + num_interval
output_Ytest = './cos424hw2_Ytest_' + num_interval
output_Ltrain = './cos424hw2_locTrain_' + num_interval
output_Ltest = './cos424hw2_locTest_' + num_interval

outXtrain = open(output_Xtrain,'w')
outYtrain = open(output_Ytrain,'w')
outXtest = open(output_Xtest,'w')
outYtest = open(output_Ytest,'w')
outLtrain = open(output_Ltrain,'w')
outLtest = open(output_Ltest,'w')

for i in features:
  
  if(label[i] == 0):
    # test data
np.savetxt(outXtest,features[i],newline=' ')
outXtest.write("\n")
outYtest.write(str(teacher[i]))
outYtest.write("\n")
outLtest.write(str(i))
outLtest.write("\n")
  elif(label[i] == 1):
    #train data
np.savetxt(outXtrain,features[i],newline=' ')
outXtrain.write("\n")
outYtrain.write(str(teacher[i]))
outYtrain.write("\n")
outLtrain.write(str(i))
outLtrain.write("\n")
  else:
    print "something odd just happened."
outXtest.close()
outYtest.close()
outXtrain.close()
outYtrain.close()
outLtrain.close()
outLtest.close()
