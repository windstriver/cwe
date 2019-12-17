import numpy as np
import pandas as pd
import os

case_dir = '../'
case_name = 'testCase-RW3'

# 6-12 s mean
UMean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UMean.xy')
UPrime2Mean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UPrime2Mean.xy')

UMean_dat = np.unique(UMean_dat, axis=0)
UPrime2Mean_dat = np.unique(UPrime2Mean_dat, axis=0)

index = ['z/H = {:3.1f}'.format(n) for n in np.linspace(0.1, 4.0, num=20, endpoint=True, dtype=float)]
columns = ['x/H = {:d}'.format(n) for n in np.linspace(1, 15, num=15, endpoint=True, dtype=int)]

UMean_12_df = pd.DataFrame(UMean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
uuMean_12_df = pd.DataFrame(UPrime2Mean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
vvMean_12_df = pd.DataFrame(UPrime2Mean_dat[:300,6].reshape((20,15), order='F'), index=index, columns = columns)
wwMean_12_df = pd.DataFrame(UPrime2Mean_dat[:300,8].reshape((20,15), order='F'), index=index, columns = columns)

# 12-18 s mean
UMean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UMean.xy')
UPrime2Mean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UPrime2Mean.xy')

UMean_dat = np.unique(UMean_dat, axis=0)
UPrime2Mean_dat = np.unique(UPrime2Mean_dat, axis=0)

index = ['z/H = {:3.1f}'.format(n) for n in np.linspace(0.1, 4.0, num=20, endpoint=True, dtype=float)]
columns = ['x/H = {:d}'.format(n) for n in np.linspace(1, 15, num=15, endpoint=True, dtype=int)]

UMean_18_df = pd.DataFrame(UMean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
uuMean_18_df = pd.DataFrame(UPrime2Mean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
vvMean_18_df = pd.DataFrame(UPrime2Mean_dat[:300,6].reshape((20,15), order='F'), index=index, columns = columns)
wwMean_18_df = pd.DataFrame(UPrime2Mean_dat[:300,8].reshape((20,15), order='F'), index=index, columns = columns)

# 18-24 s mean
UMean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UMean.xy')
UPrime2Mean_dat = np.genfromtxt(case_dir+case_name+'/postProcessing/internalCloudUMean/12/cloud_UPrime2Mean.xy')

UMean_dat = np.unique(UMean_dat, axis=0)
UPrime2Mean_dat = np.unique(UPrime2Mean_dat, axis=0)

index = ['z/H = {:3.1f}'.format(n) for n in np.linspace(0.1, 4.0, num=20, endpoint=True, dtype=float)]
columns = ['x/H = {:d}'.format(n) for n in np.linspace(1, 15, num=15, endpoint=True, dtype=int)]

UMean_24_df = pd.DataFrame(UMean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
uuMean_24_df = pd.DataFrame(UPrime2Mean_dat[:300,3].reshape((20,15), order='F'), index=index, columns = columns)
vvMean_24_df = pd.DataFrame(UPrime2Mean_dat[:300,6].reshape((20,15), order='F'), index=index, columns = columns)
wwMean_24_df = pd.DataFrame(UPrime2Mean_dat[:300,8].reshape((20,15), order='F'), index=index, columns = columns)

# 6-24 s mean
UMean_df = (UMean_12_df + UMean_18_df + UMean_24_df) / 3;
uuMean_df = (uuMean_12_df + uuMean_18_df + uuMean_24_df) / 3;
vvMean_df = (vvMean_12_df + vvMean_18_df + vvMean_24_df) / 3;
wwMean_df = (wwMean_12_df + wwMean_18_df + wwMean_24_df) / 3;

Iu_df = pd.DataFrame(np.sqrt(uuMean_df.values)/UMean_df.values, index=index, columns = columns)
Iv_df = pd.DataFrame(np.sqrt(vvMean_df.values)/UMean_df.values, index=index, columns = columns)
Iw_df = pd.DataFrame(np.sqrt(wwMean_df.values)/UMean_df.values, index=index, columns = columns)

if not os.path.exists('../results/'):
	os.makedirs('../results/')

UMean_df.to_csv('../results/U_' + case_name + '.csv', index=False, header=False, float_format='%6.2f')
Iu_df.to_csv('../results/Iu_' + case_name + '.csv', index=False, header=False, float_format='%5.3f')
Iv_df.to_csv('../results/Iv_' + case_name + '.csv', index=False, header=False, float_format='%5.3f')
Iw_df.to_csv('../results/Iw_' + case_name + '.csv', index=False, header=False, float_format='%5.3f')
