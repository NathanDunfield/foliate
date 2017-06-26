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
import edge_orient, link

def compute_euler(task):
    N = t3m.Mcomplex(str(task['laminar_tri']))
    L = link.LinkSphere(N)
    laminar_orients = [edge_orient.EdgeOrientation(N, L, signs)
                       for signs in eval(task['laminar_orients'])]
    task['taut_euler_0'] = repr([1 if eo.euler_class_vanishes() else 0
                                 for eo in laminar_orients]).replace(' ', '')
    task['done'] = True


task1 = {'name':'m376(-2, 3)',
         'laminar_tri':'oLLvLMLPQQcacgikmkimjlnnnmnkjaaagnnnwkwkw',
         'laminar_orients':'[[1,1,-1,1,1,-1,-1,1,1,1,1,-1,1,1,-1]]'}
    
taskdb2.worker.run_function('QHSpheres', 'task_fix_euler', compute_euler)
