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
import snappy, main
import snappy.snap.t3mlite as t3m
import edge_orient

def search_for_taut(task):
    M = snappy.Manifold(task['name'])
    fol = main.first_foliation(M, 1000, 1000)
    if fol is not None:
        task['taut'] = True
        task['laminar_tri'] = fol.mcomplex.name
        task['done'] = True

def count_taut(task):
    N = t3m.Mcomplex(str(task['laminar_tri']))
    laminar_orients = [eo for eo in edge_orient.edge_orientations(N) if eo.gives_foliation()]
    task['laminar_orients'] = repr([eo.signs for eo in laminar_orients]).replace(' ', '')
    task['taut_euler_0'] = repr([1 if eo.euler_class_vanishes() else 0 for eo in laminar_orients]).replace(' ', '')
    task['done'] = True

task1 = {'name':'m003(-1, 3)', 'laminar_tri':'jLLvMQQcdfigihghihsafroggnw'}
task2 = {'name':'o9_34819(5, 1)', 'laminar_tri':'uLALvvLPMQvAQQccbbeilkjpknmqoprrtsstqqnnbmxeonkvtngpfrkkk'}
    
taskdb2.worker.run_function('QHSpheres', 'task_fol', search_for_taut)
