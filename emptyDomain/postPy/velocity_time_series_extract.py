import cloud
import numpy as np
import os

case_dir = '../'
case_name = 'testCase'

ic = cloud.Cloud(case_dir+case_name+'/postProcessing/internalCloudU/', 'U', 5e-4, 6)

# cloudData = np.unique(ic.cloudData, axis=1)
UTimeSeries = ic.cloudData[:,80:90,0];
VTimeSeries = ic.cloudData[:,80:90,1];
WTimeSeries = ic.cloudData[:,80:90,2];

if not os.path.exists('../results/'):
	os.makedirs('../results/')

np.savetxt("../results/UTimeSeries_"+case_name+".csv", UTimeSeries, delimiter=',', fmt='%8.4f')
np.savetxt("../results/VTimeSeries_"+case_name+".csv", VTimeSeries, delimiter=',', fmt='%8.4f')
np.savetxt("../results/WTimeSeries_"+case_name+".csv", WTimeSeries, delimiter=',', fmt='%8.4f')
