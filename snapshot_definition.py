
from os import path
from config import *
from functions_tools import *
import yt
from yt import YTArray

import numpy as np
from yt.units import G
import array
import pandas as pd

global datos_edades
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
#global angular_momentum_ref
#angular_momentum_ref = YTArray([2.63425158e+29, 4.18786719e+28, -1.09858375e+29], "cm**2/s")


def cartesian_to_cylindrical (df):
    df["Phi"] =  np.mod(np.arctan2(df["Y"],df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    df["Vphi"] = (df["X"]*df["VY"] - df["Y"]*df["VX"])/df["R"] #todo revisar signos de phi y vphi
    df["VR"] = (df["X"]*df["VX"] + df["Y"]*df["VY"])/df["R"]
    return df
from os import path
from config import *
import yt
from yt import YTArray

import numpy as np
from yt.units import G
import array
import pandas as pd

#global datos_edades
datos_edades = pd.read_csv(path_datos + "edades.csv", sep = ",",index_col = 0)
#global angular_momentum_ref
#angular_momentum_ref = YTArray([2.63425158e+29, 4.18786719e+28, -1.09858375e+29], "cm**2/s")

def cartesian_to_cylindrical (df):
    df["Phi"] =  np.mod(np.arctan2(df["Y"],df["X"]), 2*np.pi)
    df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
    df["Vphi"] = (df["X"]*df["VY"] - df["Y"]*df["VX"])/df["R"] #todo revisar signos de phi y vphi
    df["VR"] = (df["X"]*df["VX"] + df["Y"]*df["VY"])/df["R"]
    return df


class Snapshot:
    def __init__(self, name):
        self.name = name
        self.path_snapshot = None
        self.lb = None
        self.center = None
        self.Rvir = None
        self.ds = None 
        self.dm = None
        self.gas = None
        self.stars = None
        self.disk = None


        def read_lb():
            self.lb = datos_edades.loc[datos_edades['Snapshot'] == self.name, 'Lookback'].iloc[0]

        def read_center_Rvir ():
            centro = np.loadtxt(path_datos +f'center_{self.name}.txt')
            center = YTArray([centro[0], centro[1], centro[2]], "cm")
        #   cx,cy,cz = center[0].in_units("cm"), center[1].in_units("cm"),  center[2].in_units("cm")
            Rvir = YTArray(centro[3], "kpc")
            self.center = center
            self.Rvir = Rvir

        def find_path_for_yt():
            # name = snapshots_analysis[i]
            if self.name < 425:
                path_snapshot = "/media/temp1/bego/GARROTXA_ART/"
            elif (self.name >= 425)&(name < 600):
                path_snapshot = "/srv/cab1/garrotxa/GARROTXA_ART/MW_003/RUN2.2/"
            elif (self.name >=600 )&(name < 800):
                path_snapshot = "/home/Garrotxa_ART/New_Run/"
            elif (self.name >= 800) & (name < 900) :
                path_snapshot = "/media/temp/bego/New_Resim/"
            elif self.name >= 900 :
                path_snapshot = "/media/temp1/GARROTXA_ART/MW_003/RUN2.2/"
            self.path_snapshot = path_snapshot
        
        print(f"Initializing snapshot {name}")
        find_path_for_yt()
        read_lb()
        print(f"Lookback time: {self.lb} Gyr")
        read_center_Rvir()
        
    def load_stars (self):
        self.stars = pd.read_csv(path_csv + f"{self.name}_stars_Rvir.csv",sep = ",")
        self.stars = cartesian_to_cylindrical(self.stars)

    def load_dm (self):
        self.dm = pd.read_csv(path_csv + f"{self.name}_dm_Rvir.csv",sep = ",")
        self.dm = cartesian_to_cylindrical(self.dm)

    def load_gas (self):
        self.gas = pd.read_csv(path_csv + f"Gas_{self.name}.csv",sep = ",")
        self.gas = cartesian_to_cylindrical(self.gas)

    def load_disk(self):
        self.disk = pd.read_csv(path_disk + f"Stars_disco_{self.name}.csv")

    def filter_disk_particles(self):
        dfA = self.stars[self.stars['ID'].isin(self.disk["ID"])]
        df = dfA[(dfA['R']< 25)].copy()
        df["Phi"] = np.mod(np.arctan2(df["Y"], df["X"]), 2*np.pi)
        df["R"] = np.sqrt(df["X"]**2 + df["Y"]**2)
        return df



    def angular_momentum (self):
        
        self.dm_am_total = np.array(calculate_angular_momentum(self.ds, self.center, 120, self.Rvir, "darkmatter"))
        self.dm_am_inner = np.array(calculate_angular_momentum(self.ds, self.center, 20, self.Rvir, "darkmatter"))
        self.stars_am = np.array(calculate_angular_momentum(self.ds, self.center, 20, self.Rvir, "young_stars"))
        self.stars_am_inner = np.array(calculate_angular_momentum(self.ds, self.center, 6, self.Rvir, "young_stars"))

    def save_angular_momentum (self, component):
        print(path_save_data  + f"angular_momentum_{component}.csv")
        am_table = pd.read_csv(path_save_data  + f"angular_momentum_{component}.csv" , index_col = 0, sep = ",")
        print("Saving angular momentum")
        am_values = getattr(self, f"{component}" )
        print(am_values)
        new_row ={'Snapshot':self.name, "Lookback": self.lb, "am_x": am_values[0], "am_y":am_values[1], "am_z":am_values[2]}
        am_table = am_table.append(new_row, ignore_index = True)
        am_table = am_table.sort_values(by=["Snapshot"], ascending = True)
        am_table.to_csv(path_save_data  +f"angular_momentum_{component}.csv", sep = ",")

    def read_angular_momentum(self):
        
        am_table_stars = pd.read_csv(path_save_data + f"angular_momentum_stars_am.csv" , index_col = 0, sep = ",")
        am_table_stars_inner = pd.read_csv(path_save_data + f"angular_momentum_stars_am_inner.csv" , index_col = 0, sep = ",")
        am_table_dm_inner = pd.read_csv(path_save_data + f"angular_momentum_dm_am_inner.csv" , index_col = 0, sep = ",")
        am_table_dm_total = pd.read_csv(path_save_data + f"angular_momentum_dm_am_total.csv" , index_col = 0, sep = ",")
        
        self.dm_am_inner = am_table_dm_inner.loc[am_table_dm_inner['Snapshot'] == self.name, ['am_x', 'am_y', 'am_z']].iloc[0]
        self.dm_am_total = am_table_dm_total.loc[am_table_dm_total['Snapshot'] == self.name, ['am_x', 'am_y', 'am_z']].iloc[0]
        self.stars_am = am_table_stars.loc[am_table_stars['Snapshot'] == self.name, ['am_x', 'am_y', 'am_z']].iloc[0]
        self.stars_am_inner = am_table_stars_inner.loc[am_table_stars_inner['Snapshot'] == self.name, ['am_x', 'am_y', 'am_z']].iloc[0]

    def plot_angular_momentum_quiver (self, angle):
        matriz_rotacion = calculate_rotation_matrix(angular_momentum_ref)
        am_stars_rotated = apply_transformation_matrix(matriz_rotacion, self.stars_am[0],self.stars_am[1],self.stars_am[2])
        am_dm_total_rotated = apply_transformation_matrix(matriz_rotacion, self.dm_am_total[0],self.dm_am_total[1],self.dm_am_total[2])
        am_dm_inner_rotated = apply_transformation_matrix(matriz_rotacion, self.dm_am_inner[0],self.dm_am_inner[1],self.dm_am_inner[2])
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        ax.quiver(0, 0, 0, am_stars_rotated[0], am_stars_rotated[1],am_stars_rotated[2], length=0.1, 
        normalize=True, color = "red", label = "Stars")
        ax.quiver(0, 0, 0, am_dm_total_rotated[0], am_dm_total_rotated[1], am_dm_total_rotated[2], length=0.1, 
        normalize=True, color = "blue", label = "Dm total")
        ax.quiver(0, 0, 0, am_dm_inner_rotated[0], am_dm_inner_rotated[1], am_dm_inner_rotated[2], length=0.1,
         normalize=True, color = "black",  label = "Dm inner")
        ax.view_init(elev=25., azim=120)
        ancho = 0.08
        ax.axes.set_xlim3d(left=-ancho, right=ancho) 
        ax.axes.set_ylim3d(bottom=-ancho, top=ancho) 
        ax.axes.set_zlim3d(bottom=ancho-ancho/2, top=-ancho)
        ax.legend()
        ax.xaxis.set_tick_params(labelsize=4)
        ax.yaxis.set_tick_params(labelsize=4)
        ax.zaxis.set_tick_params(labelsize=4)
        ax.set_title (f"Lookback time {self.lb:.2f} Gyr")
        ax.view_init(elev=50., azim=0)
        ax.text(0.05, 0, 0, "", color='red')
        plt.savefig(path_save_angular_momentum_figs + f"am_vectors_{self.name}.png")







