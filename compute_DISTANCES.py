import pytraj as pt
import numpy as np

 
def CH_distances(traj, time_averaged=True, save_to_file=True, plot=True):
    """ return a python array of all the distances between the Carbon-sp3 holding 
    the fluorinated structure (AtomName: 'C') and the 8 hydroxyl groups represented by
    alchoholic hydrogen (AtomType: 'ho') as a function of simulation timesteps """
    print('Start CH_distances')
    
    my_hydrogens_names= [atoms.name for atoms in traj.top.atoms if atoms.resid==0 if atoms.type=='ho']
    
    distances = []

    if traj.ndendr > 1:
        for dendrons in range(1,traj.ndendr+1):
            print("I'm in dendron"+str(dendrons))
            mymask=[':'+str(dendrons)+'@C'+ ' :' + str(dendrons) + '@'+ hydrogens for hydrogens in my_hydrogens_names]
            distances.append(pt.distance(traj, mymask))
    else:
        my_hydrogen_names = [atoms.name for atoms in traj.top.atoms if atoms.resid== 0 if atoms.type=='ho']
        distances = {}
        for hydrogen_i in my_hydrogen_names:
            certain_CH_couple = ':1@C :1@' + hydrogen_i
            distances[hydrogen_i] = pt.distance(traj, certain_CH_couple)
    return my_hydrogens_names, distances





##########################################################################################################

def FF_distances(traj):
    """ return a python array of all the INTER-molecular distances between
    fluorine atoms in the dendrons present in the simulation box"""
    if traj.ndendr < 2:
        return 
    #fluorine_dendron1=traj.fluorine_atomnames(0)
    #fluorine_dendron2=traj.fluorine_atomnames(1)
    fluorine_dendron1=([atom for atom in traj.top.atoms if atom.type=='f' if atom.resid==0])
    fluorine_dendron2=([atom for atom in traj.top.atoms if atom.type=='f' if atom.resid==1])
    mymask = [':'+str(fluorines1.resid+1)+'@'+str(fluorines1.name)+' :'+str(fluorines2.resid+1)+'@'+str(fluorines2.name) \
    for fluorines1 in fluorine_dendron1 for fluorines2 in fluorine_dendron2]
    distances=pt.distance(traj, mymask)
    return distances
