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
import snappy, edge_orient, util, peripheral

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

def add_alex(task):
    M = snappy.Triangulation(task['name'])
    alex = M.alexander_polynomial()
    task['alex'] = repr(alex).replace(' ', '')
    coeffs = alex.coefficients()
    if (coeffs != (len(coeffs)//2)*[1, -1] + [1] or coeffs == [1] or
        (coeffs == [1, -1, 1] and len(M.link()) != 3)):
        task['floer_simple'] = False
    elif :
        task[
    task['done'] = True

task1 = {'name':'K10a12'}
task2 = {'name':'K7a7'}
#search_for_persist(task1)


taskdb2.worker.run_function('KnotsInS3', 'task_alex', add_alex)
#taskdb2.worker.run_function('KnotsInS3NonAlt16', 'task_persist', search_for_persist)
