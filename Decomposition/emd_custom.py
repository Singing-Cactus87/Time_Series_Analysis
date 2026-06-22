# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks

df = pd.read_csv("Example2.csv",header=0) #ECOS PPI, CPI dataset

df.head()

class EMD_custom():
  def __init__(self,sd_thres):
    super(EMD_custom,self).__init__()
    self.sd_thres = sd_thres

  def emd_compute(self,x): #x: dataframe

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

    current_dt = x
    sd = 10; C = []; k= 0
    _1, _2 = find_peaks(current_dt)
    __1, _2 = find_peaks((-1)*current_dt)
    num = len(_1)+len(__1)
    while num > 2:
      U_,L_,pk_num = Envelope_extract(current_dt)
      M_ = np.mean([U_,L_],axis=0)
      H_ = current_dt - M_; H_ = pd.Series(H_); h_current = H_
      u1,l1,pk_n = Envelope_extract(h_current)
      num_UL = pk_n-4
      SD = 100; zero_cross = len(np.where(np.diff(np.sign(np.array(h_current).reshape(-1))))[0])
      while ((SD > self.sd_thres) | (np.abs(zero_cross-num_UL) > 1)):
        h_new = h_current-np.mean([u1,l1],axis=0)
        SD = np.sum(np.abs(np.array(h_current-h_new).reshape(-1))**2/(np.array(h_current).reshape(-1)**2+1e-6))
        h_current = h_new
        u1,l1,pk_n = Envelope_extract(h_current)
        num_UL = pk_n-4
        zero_cross = len(np.where(np.diff(np.sign(np.array(h_current))))[0])
      C_ = h_current; C.append(np.array(C_).reshape(-1))
      current_dt = current_dt - C_
      _1, _2 = find_peaks(current_dt)
      __1, _2 = find_peaks((-1)*current_dt)
      num = len(_1)+len(__1)
      k+= 1; print(k)

    residue = current_dt

    return C, residue

start_emd = EMD_custom(sd_thres=0.2)

C, res = start_emd.emd_compute(df['PPI'])

len(C)

fig, axs = plt.subplots(2, 2, figsize=(8, 6))
axs = axs.flatten()

for a in range(4):
  ax = axs[a]
  if a != 3:
    ax.plot(C[a],color="darkcyan",label=f"IMF_{a}")
    ax.set_title(f"IMF_{a}")
    ax.legend()
  else:
    ax.plot(res,color="darkcyan",label="Residue")
    ax.set_title("Residue")
    ax.legend()

