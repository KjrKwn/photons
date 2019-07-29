import numpy as np
import pandas as pd
import copy

__all__ = ["Mesh", "read_mesh_data"]


class Mesh(object):
    """
    Type containing mesh data replaced from SPH data
    
    members
    =======
    radius :
    thetas :
    phi    :
    mass   : mass coordinate
    data   : [N_r, N_theta, N_phi, 62]
    
    data
    ====
    0  x
    1  y
    2  z
    3  r_out
    4  r_plusHalf
    5  v_r_from_fit
    6  v_x
    7  v_y
    8  v_z
    9  v_r
    10 v_theta
    11 v_phi
    12 density
    13 volume
    14 mass in the mesh
    15-27 X_alpha[i]
    28-57 X_stable[i]
    58-61 X_active[i]
    """
    
    def __init__(self):
        self.radius = None
        self.thetas = None
        self.phis = None
        self.data = None
        
    def __getitem__(self, key):
        new_mesh = copy.deepcopy(self)
        if type(key) != tuple:
            # only 1st index is given
            new_mesh.radius = np.array([self.radius[key]]).flatten()
        else:
            N = len(key)
            new_mesh.radius = np.array([self.radius[key[0]]]).flatten()
            if (N > 1):
                new_mesh.thetas = np.array([self.thetas[key[1]]]).flatten()
                if (N > 2):
                    new_mesh.phis = np.array([self.phis[key[2]]]).flatten()
        flatten_data = np.array(self.data[key])
        new_mesh.data = flatten_data.reshape(new_mesh.times.size, new_mesh.thetas.size, new_mesh.phis.size, flatten_data.size / (new_mesh.times.size * new_mesh.thetas.size * new_mesh.phis.size))
        return new_mesh

        
def read_mesh_data(filepath):
    """
    read mesh.dat file and return Mesh class object
    """
    N_col = 66
    colarray = np.linspace(1, N_col, num=N_col).astype(int).astype(str)
    colarray.astype(int).astype(str)
    col = np.array(['c'] * N_col)
    colarray2 = np.core.defchararray.add(col, colarray)

    df_ = pd.read_csv(filepath, delim_whitespace=True, header=None, names=colarray2, skiprows=9)
    
    mesh = Mesh()
    mesh.radius = np.unique(df_["c2"])
    mesh.thetas = np.unique(df_["c3"])
    mesh.phis   = np.unique(df_["c4"])
    mesh.data   = df_.to_numpy()[:,4:].reshape(mesh.radius.size, mesh.thetas.size, mesh.phis.size, 62)
    mesh.mass   = np.sum(mesh.data[:,:,:,14], axis=(1,2)).cumsum()
    return mesh
    
