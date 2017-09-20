import numpy as np
import periodictable as pt
pt._load_covalent_radius()
cov_rad = pt.covalent_radius.Cordero.split('\n')
import matplotlib.pyplot as plt
plt.ion()
# cov_rad2 = [x.split() for x in cov_rad ]
cov_rad2 = []

for el in cov_rad:
    el_sp = el.split()
    inner_list = []
    append_bool = True
    for i,el_sp_el in enumerate(el_sp):
        try:
            res = float(el_sp_el)
        except:
            res = el_sp_el
            if i == 0:
                append_bool = False
                break
        if i==2 or i == 0:
            inner_list.append(res)
    if append_bool:
        cov_rad2.append(inner_list)

cov_rad_arr = np.array(cov_rad2)
np.savetxt('cov_radii.dat',cov_rad_arr[:,1])
plt.plot(cov_rad_arr[:,0],cov_rad_arr[:,1])