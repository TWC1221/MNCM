import os
import sys
import numpy as np

### Run code for a given set of inputs
# Inputs
PolyOrders = np.arange(1,5)
Refinement = np.arange(2,12)

Disc = sys.argv[1]
Corr = 'None'
Curvilinear = True
MaxIter = 1
Preconditioner = 0
Source = 'Eigenvalue'

problem = 'htr'
NumMaterials = 3
NumGroups = 7

MeshPath = f'../PhD_Meshes/2D/{problem}/{problem}'
MaterialPath = 'benchmarks/htr.dat'
OutVTKPath = 'output_files/TEST.vtk'

FluxInt = False
FluxAvg = False
Cond = True
keff = False
tol = 1e-10

for i in range(len(PolyOrders)):
    OutDataPath = f'../Paper/NESEM_Diffusion/{problem}/{Disc}_{PolyOrders[i]}_{Corr}'

    if Curvilinear:
        OutDataPath += '_curvilinear'

    if FluxInt:
        with open(OutDataPath + '_int_flux.txt', 'w') as f:
            f.close()
    if FluxAvg:
        with open(OutDataPath + '_mean_flux.txt', 'w') as f:
            f.close()
    if Cond:
        with open(OutDataPath + '_cond_number.txt', 'w') as f:
            f.close()
    if keff:
        with open(OutDataPath + '_keff.txt', 'w') as f:
            f.close()

    for j in range(len(Refinement)):

        InputFile = f'run_{Disc}_{PolyOrders[i]}.txt'
        with open(InputFile, 'w') as f:
            
            print(f'InputFile: {InputFile}')
            if Disc in ['FEM', 'SEM']:
                if Curvilinear:
                    temp_mesh_path = MeshPath + f'_{PolyOrders[i]}_{Refinement[j]}_{Disc}_{Corr}_curvilinear.asmg'
                else:
                    temp_mesh_path = MeshPath + f'_{PolyOrders[i]}_{Refinement[j]}_{Disc}_{Corr}.asmg'
            else:
                if Curvilinear:
                    temp_mesh_path = MeshPath + f'_{PolyOrders[i]}_{Refinement[j]}_{Disc}_curvilinear.asmg'
                else:
                    temp_mesh_path = MeshPath + f'_{PolyOrders[i]}_{Refinement[j]}_{Disc}.asmg'

            print(f'Meshpath: {temp_mesh_path}')

            f.write(f'MeshPath        {temp_mesh_path}' + '\n')
            f.write(f'MatPath         {MaterialPath}' + '\n')
            f.write(f'OutVTKPath      {OutVTKPath}' + '\n')
            f.write(f'OutDataPath     {OutDataPath}' + '\n')
            f.write(f'OutVTK          False' + '\n')
            f.write(f'OutData         False' + '\n')
            f.write(f'OutMaterial     False' + '\n')

            if FluxInt:
                f.write(f'OutFluxInt      True' + '\n')
            else:
                f.write(f'OutFluxInt      False' + '\n')

            if FluxAvg:
                f.write(f'OutFluxAvg      True' + '\n')
            else:
                f.write(f'OutFluxAvg      False' + '\n')

            if keff:
                f.write(f'CalcKeff        True' + '\n')
            else:
                f.write(f'CalcKeff        False' + '\n')

            f.write(f'CalcDiffCoeff   False' + '\n')
            f.write(f'Method          Diffusion' + '\n')
            f.write(f'CoordSystem     0' + '\n')
            f.write(f'NoCores         {10}' + '\n')
            f.write(f'Preconditioner  {Preconditioner}' + '\n')
            f.write(f'Source          {Source}' + '\n')
            f.write(f'Groups          {NumGroups}' + '\n')
            f.write(f'Materials       {NumMaterials}' + '\n')
            f.write(f'SN              {2}' + '\n')
            f.write(f'Anisotropic     0' + '\n')
            f.write(f'Adjoint         False' + '\n')
            f.write(f'Tolerance       {tol}' + '\n')
            f.write(f'MaxIterations   {MaxIter}' + '\n')

            if Cond:
                f.write(f'CalcCondition   True' + '\n')
            else:
                f.write(f'CalcCondition   False' + '\n')

            f.write(f'TerminalOutput  True' + '\n')

            f.close()

        os.system(f'./output {InputFile}')
        os.remove(InputFile)        