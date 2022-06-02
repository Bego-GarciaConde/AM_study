# -*- coding: utf-8 -*-
#
"""
Created on 29/04/2022
@author: B. Garcia-Conde
"""
import gc
import yt
from yt import YTArray


import numpy as np
from yt.units import G
import array
import pandas as pd

import gc


from config import *
from snapshot_definition import *

with open(path_save_data+ 'angle_difference_am.csv','w') as csvfile:
    csvfile.write("")
   # np.savetxt(csvfile, data, fmt = "%.3f", newline=" ")
csvfile.close()



def iniciar_angular_momentum (component, r_max):
 #if iniciar == 1:
    dato_vacio = []
    df_vacio = pd.DataFrame(dato_vacio, columns=['Snapshot', "Lookback", "am_x", "am_y", "am_z"])
    df_vacio.to_csv(path_save_data+ f'angular_momentum_{component}_{r_max}.csv', sep = ",")

def angular_momentum_vector_csv (X, VX, Y, VY, Z, VZ, mass):  
    Lx = Y*VZ-Z*VY
    Ly = Z*VX-X*VZ
    Lz = X*VY-Y*VX
    #finding median L vector only with selected particles
    mLx,mLy,mLz= np.average(Lx, weight=mass),np.average(Ly, weight=mass),np.average(Lz, weight = mass)

    m=np.sqrt(mLx*mLx+mLy*mLy+mLz*mLz)
#        print('0',np.round((mLx,mLy,mLz),2))
    mLx,mLy,mLz=-mLx/m,-mLy/m,-mLz/m #normalization of the median L vector         
    v3=[mLx,mLy,mLz]
    return v3

def angular_momentum_comparison (snapshot, snapshot_pre):
        #csvfile.close()
        print(f"Calculating angular momentum for {snapshot.name}")
        # try:
        #     snapshot.read_angular_momentum()
        # except:
        snapshot.load_dm()
        snapshot.load_stars()
        df_dm = snapshot.dm 
        df_stars = snapshot.stars
        df_dm["R_sph"] = np.sqrt(df_dm["X"]**2 + df_dm["Y"]**2 + df_dm["Z"]**2)
        df_stars["R_sph"] = np.sqrt(df_stars["X"]**2 + df_stars["Y"]**2 + df_stars["Z"]**2)
        df_dm = df_dm[(df_dm["R_sph"]<r_max)]
        df_stars = df_stars[(df_stars["R_sph"]<r_max)]
        

        v_dm = angular_momentum_vector_csv(df_dm["X"], df_dm["VX"],df_dm["Y"], df_dm["VY"],
                                                  df_dm["Z"], df_dm["VZ"], df_dm["Mass"] )
        v_stars = angular_momentum_vector_csv(df_stars["X"], df_stars["VX"],df_stars["Y"], df_stars["VY"],
                                                      df_stars["Z"], df_stars["VZ"], df_stars["Mass"]  )
      #   snapshot_pre.angular_momentum()
        time = snapshot_pre.lb - snapshot.lb

        snapshot.save_angular_momentum(component="dm_am_total")
        snapshot.save_angular_momentum(component="dm_am_inner")
        snapshot.save_angular_momentum(component="stars_am")
        snapshot.save_angular_momentum(component = "stars_am_inner")
        #snapshot.plot_angular_momentum_quiver()

        angle_reference_theta = [0,0,1]
        angle_reference_phi = [1,1,0]
        # angle_am_stars.append(angle_between_vectors(angle_reference, snapshot.stars_am))
        # angle_am_dm_inner.append(angle_between_vectors(angle_reference, snapshot.dm_am_inner))
        # rotation_matrix = calculate_rotation_matrix(angular_momentum_ref)
        # vector_am_stars = apply_transformation_matrix(rotation_matrix,snapshot.stars_am[0], snapshot.stars_am[1], snapshot.stars_am[2])
        # vector_am_dm_inner = apply_transformation_matrix(rotation_matrix,snapshot.dm_am_inner[0], snapshot.dm_am_inner[1], snapshot.dm_am_inner[2])
        # vector_am_dm_total = apply_transformation_matrix(rotation_matrix,snapshot.dm_am_total[0], snapshot.dm_am_total[1], snapshot.dm_am_total[2])
        
        # a_stars = angle_between_vectors(angular_momentum_ref, snapshot.stars_am)
        # a_dm_inner = angle_between_vectors(angular_momentum_ref_dm, snapshot.dm_am_inner)
        # a_difference_total = angle_between_vectors(angular_momentum_ref_dm, snapshot.dm_am_total)

    #     a_stars_theta = angle_between_vectors(angle_reference_theta, snapshot.stars_am)
    #     a_dm_inner_theta = angle_between_vectors(angle_reference_theta, snapshot.dm_am_inner)
    #     a_difference_total_theta = angle_between_vectors(angle_reference_theta, snapshot.dm_am_total)

    #     a_stars_phi = angle_between_vectors(angle_reference_phi, snapshot.stars_am)
    #     a_dm_inner_phi = angle_between_vectors(angle_reference_phi, snapshot.dm_am_inner)
    #     a_difference_total_phi = angle_between_vectors(angle_reference_phi, snapshot.dm_am_total)

    #     a_difference = angle_between_vectors(snapshot.stars_am, snapshot.dm_am_inner)
    #     a_stars_diff = angle_between_vectors(snapshot.stars_am, snapshot.stars_am_inner)
        
    #     a_diff_dm_pre = angle_between_vectors(snapshot_pre.dm_am_inner, snapshot.dm_am_inner)
    #     a_diff_stars_pre = angle_between_vectors(snapshot_pre.stars_am, snapshot.stars_am)

    #     data = np.array([snapshot.name, snapshot.lb,a_difference, a_stars_diff, a_stars_phi, a_stars_theta, 
    #                   a_dm_inner_phi, a_dm_inner_theta, a_difference_total_phi,a_difference_total_theta, a_diff_dm_pre, a_diff_stars_pre])
    #    # data = np.array([snapshot.name, snapshot.lb,a_stars, a_dm_inner, a_difference, a_difference_total, a_stars_diff])
    #     with open(path_save_data+ 'angle_difference_am.csv','a+') as csvfile:
    #         csvfile.write("\n")
    #         np.savetxt(csvfile, data, fmt = "%.3f", newline=" ")
    #     csvfile.close()

