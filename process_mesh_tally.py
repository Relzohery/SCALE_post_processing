import numpy as np

def readTally(fileName, npg):
    lines = open(fileName).readlines()
    tallyName       = lines[0]
    Nx, Ny, Nz      = tuple(int(i) for i in lines[1].split())
    total           = Nx*Ny*Nz
    meshes          = []
    data            = np.zeros((total, npg))
    data_sigma      = data.copy()
    for no, line in enumerate(lines[2:]):
        if no < Nx + Ny + Nz + 3:
            meshes.append(float(line))
        for group in range(1, npg+1):
            if f', group {group}' in line:
                print(line, f'geoup {group}', no, no+2*total)
                data_1 = [float(i) for i in lines[no+4: (no+4)+2*total]]
                data[:, group-1]       = [float(i) for i in lines[no+4: (no+4)+total]]
                data_sigma[:, group-1] = [float(i) for i in lines[(no+4)+total: (no+4)+2*total]]

    meshX = meshes[:Nx]
    meshY = meshes[Nx+1: Nx+Ny+1]
    meshZ = meshes[Nx+Ny+2: -1]
  
    out = {'Nx': Nx, 
          'Ny': Ny, 
          'Nz': Nz, 
          'data': data, 
          'data_sigma':data_sigma, 
          'meshX': meshX, 
          'meshY': meshY,
          'meshZ': meshZ, 
          'tallyName': tallyName
        }
    
    return out