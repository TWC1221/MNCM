import os 
import numpy as np

### Run code for a given set of inputs
# Inputs
PolyOrders = np.arange(1,5)
Refinement = np.arange(1,11)

Discs = ['FEM']
Corr = 'None'
Curvilinear = False
MaxIter = 1
Preconditioner = 0
Source = 'Fixed'

problem = 'two_region_mms'
NumMaterials = 2
NumGroups = 1

MeshPath = f'../PhD_Meshes/2D/{problem}/{problem}'
MaterialPath = 'benchmarks/two_mat.dat'
OutVTKPath = 'output_files/TEST.vtk'

FluxInt = True
Cond = False

if FluxInt and Cond:
    raise Warning('Cannot output both flux integral and condition number currently...')

for Disc in Discs:
    for i in range(len(PolyOrders)):

        if FluxInt:
            OutDataPath = f'/home/ciaran/Documents/Paper/NESEM_Diffusion/two_region_mms/{Disc}_{PolyOrders[i]}_fluxint_{Corr}'
            if Curvilinear:
                OutDataPath+=f'_curvilinear.txt'
            else:
                OutDataPath+=f'.txt'

            with open(OutDataPath, 'w') as f:
                f.write(f'poly, ref, fluxes, num_nodes' +'\n')
                f.close()

        elif Cond:
            OutDataPath = f'/home/ciaran/Documents/Paper/NESEM_Diffusion/two_region_mms/{Disc}_{PolyOrders[i]}_cond_{Corr}'
            if Curvilinear:
                OutDataPath+=f'_curvilinear.txt'
            else:
                OutDataPath+=f'.txt'

            with open(OutDataPath, 'w') as f:
                f.write('poly, ref, lambda_min, lambda_max, cond_number, num_nodes, num_elements' +'\n')
                f.close()

        for j in range(len(Refinement)):

            # if Curvilinear:
            #     OutDataPath = f'/home/ciaran/Documents/Paper/NESEM_Diffusion/two_region_mms/L2_files/{Disc}_{PolyOrders[i]}_{Refinement[j]}_{Corr}_curvilinear.txt'
            # else:
            #     OutDataPath = f'/home/ciaran/Documents/Paper/NESEM_Diffusion/two_region_mms/L2_files/{Disc}_{PolyOrders[i]}_{Refinement[j]}_{Corr}.txt'

            # with open(OutDataPath, 'w') as f:
            #     f.write('poly, ref, lambda_min, lambda_max, cond_number, num_nodes, num_elements' +'\n')
            #     f.close()

            with open('run.txt', 'w') as f:
                
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
                f.write(f'CalcDiffCoeff   True' + '\n')
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
                f.write(f'Tolerance       1e-8' + '\n')
                f.write(f'MaxIterations   {MaxIter}' + '\n')

                if Cond:
                    f.write(f'CalcCondition   True' + '\n')
                else:
                    f.write(f'CalcCondition   False' + '\n')

                f.write(f'TerminalOutput  True' + '\n')

                f.close()

            os.system('./output')