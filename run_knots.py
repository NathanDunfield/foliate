#! /bin/env sage-python
#
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=3000
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=slurm_out/%j
#SBATCH --error=slurm_error/%j

import taskdb2.worker
import edge_orient, util, peripheral

def search_for_persist(task):
    for isosig in util.cusped_isosigs(task['name']):
        M = peripheral.Triangulation(isosig)
        try:
            for eo in edge_orient.edge_orientations(M):
                if eo.gives_foliation():
                    task['foliar_tri'] = isosig
                    task['foliar_orient'] = repr(eo.signs).replace(' ', '')
                    task['done'] = True
                    return
        except AssertionError:
            return 
            

task1 = {'name':'K10a12'}

taskdb2.worker.run_function('KnotsInS3NonAlt15', 'task_persist', search_for_persist)
