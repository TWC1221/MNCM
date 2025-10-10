import os 
import numpy as np

### Run code for a given set of inputs
# Inputs
# PolyOrders = [2, 4, 6, 8]
Refinement = [1]
PolyOrders = [1]
Discs = ['FEM']
Corr = 'None'
Curvilinear = True

NCores = 10

MaxIter = 1
Preconditioner = 0
Source = 'Eigenvalue'
Groups = 7

MeshPath = f'../PhD_Meshes/2D/C5G7/C5G7'
MaterialPath = 'materials/C5G7.dat'
OutVTKPath = 'output_files/TEST.vtk'
OutDataPath = '../outputs/TEST_data'

for Disc in Discs:
    for i in range(len(PolyOrders)):
        OutDataPath = f'../outputs/C5G7_{Disc}_{str(PolyOrders[i])}'
        
        for j in range(len(Refinement)):
            with open('run.txt', 'w') as f:
                
                if Disc in ['FEM', 'SEM']:
                    temp_mesh_path = MeshPath + f'_{str(PolyOrders[i])}_{str(Refinement[j])}_{Disc}_{Corr}'
                else:
                    temp_mesh_path = MeshPath + f'_{str(PolyOrders[i])}_{str(Refinement[j])}_{Disc}'

                if Curvilinear:
                    temp_mesh_path += '_Curvilinear.asmg'
                else:
                    temp_mesh_path += '.asmg'

                f.write(f'MeshPath        {temp_mesh_path}' + '\n')
                f.write(f'MatPath         {MaterialPath}' + '\n')
                f.write(f'OutVTKPath      {OutVTKPath}' + '\n')
                f.write(f'OutDataPath     {OutDataPath}' + '\n')
                f.write(f'OutVTK          False' + '\n')
                f.write(f'RUNANALYSIS     True' + '\n')
                f.write(f'OutData         False' + '\n')
                f.write(f'OutMaterial     False' + '\n')
                f.write(f'OutFluxInt.     True' + '\n')
                f.write(f'OutKEFF         True' + '\n')
                f.write(f'CalcDiffCoeff   True' + '\n')
                f.write(f'Method          Diffusion' + '\n')
                f.write(f'CoordSystem     0' + '\n')
                f.write(f'Preconditioner  {Preconditioner}' + '\n')
                f.write(f'Source          {Source}' + '\n')
                f.write(f'Groups          {Groups}' + '\n')
                f.write(f'SN              {2}' + '\n')
                f.write(f'Anisotropic     0' + '\n')
                f.write(f'Adjoint         False' + '\n')
                f.write(f'Tolerance       1e-8' + '\n')
                f.write(f'MaxIterations   {MaxIter}' + '\n')
                f.write(f'CalcCondition   False' + '\n')
                f.write(f'TerminalOutput  True' + '\n')

                f.close()

            os.system(f'./akatosh run.txt omp={NCores}')