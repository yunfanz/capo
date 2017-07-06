import numpy as np, aipy as ap

#Generates baseline informations for Adrian's chromaticity calculations. 

cal='psa6622_v003'
fqs = np.linspace(0.145, 0.155, 10)
aa = ap.cal.get_aa(cal, fqs)

t0 = 2456681.5
t1 = 2456681.5054
bl1 = (96,102)
bl2 = (96,94)
# bl1 = (0,26)
# bl2 = (0,38)

aa.set_jultime(t0)
bl1v = aa.get_baseline(bl1[0],bl1[1],'z')*ap.const.len_ns/100
bl2v = aa.get_baseline(bl2[0],bl2[1],'z')*ap.const.len_ns/100

m = np.linalg.pinv(aa.eq2top_m)
bl1eq = np.dot(m, bl1v)
bl2eq = np.dot(m, bl2v)

aa.set_jultime(t1)
m = np.linalg.pinv(aa.eq2top_m)
bl1top = np.dot(aa.eq2top_m, bl1eq)
bl2top = np.dot(aa.eq2top_m, bl2eq)

print bl1v
print bl2v
print bl1top
print bl2top