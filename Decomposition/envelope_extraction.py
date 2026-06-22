# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

df = pd.read_csv("Example2.csv",header=0) #ECOS PPI, CPI dataset

df.head()

def Envelope_extract(data):
  idx_max, _ = find_peaks(data)
  idx_min, _2 = find_peaks(-1*data)

  idx_max = np.unique(np.concatenate([[0],idx_max,[len(data) - 1]],axis=0))
  idx_min = np.unique(np.concatenate([[0],idx_min,[len(data) - 1]],axis=0))

  lc_max = data[idx_max]
  lc_min = data[idx_min]


  t_lc_max= idx_max#time stamp
  x_lc_max = np.array(lc_max).reshape(-1) #value
  cs_lc_max = CubicSpline(t_lc_max,x_lc_max)

  t_lc_min= idx_min#time stamp
  x_lc_min = np.array(lc_min).reshape(-1) #value
  cs_lc_min = CubicSpline(t_lc_min,x_lc_min)

  ts = np.linspace(0, len(np.array(data).reshape(-1))-1, len(np.array(data)))

  u_envelope = cs_lc_max(ts, nu=0)
  l_envelope = cs_lc_min(ts, nu=0)
  num_peaks = len(lc_max)+len(lc_min)
  return u_envelope, l_envelope, num_peaks

_1 , _2, _3 = Envelope_extract(df['PPI'])

ts = np.linspace(0, len(np.array(df['PPI']).reshape(-1))-1, len(np.array(df['PPI'])))

sns.set_style("darkgrid")
plt.plot(ts,df['PPI'],color="black",label="original")
plt.plot(ts,_1,color="red",linestyle="--",label="upper_envelope")
plt.plot(ts,_2,color="blue",linestyle="--",label="lower_envelope")
plt.legend();
