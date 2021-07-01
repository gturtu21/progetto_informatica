import os
import sys
import numpy as np
import pytraj as pt
import matplotlib.pyplot as plt
import pickle 
############
try:
    sys.path.remove('/usr/local/amber16/lib/python2.7/site-packages')

except ValueError:
    print('The path variable was already correct!!')
    pass
################################# IMPORT THE TRAJECTORY ##################################

class my_traj(pt.TrajectoryIterator):
    """ cosolvent should be 'ETN' or 'TFN' or None """
    from mytrajectories import get_trajectory_files
    all_trajectories= get_trajectory_files()
    def __init__(self, trajectory_index, cosolvent, number_of_dendrons, timerange=[(0,-2)]):
        self.filepath, self.filetraj, self.filetop = self.all_trajectories[trajectory_index] 
        self.index = trajectory_index+1
        self.cosolvent = cosolvent
        self.ndendr = number_of_dendrons
        self.time = timerange
        """ the following command it's used to instantiate a pt.TrajectoryIterator object:
        https://amber-md.github.io/pytraj/latest/_api/pytraj.trajectory_iterator.html"""
        super().__init__(self.filepath+self.filetraj,self.filepath+self.filetop)
    
    #################

    def get_OH_hydrogens(self):
        """ return a list of the atom names of hydrogens in the -OH groups
        present in the dendrons, this method is called to compute distances """
        return [atoms.name for atoms in self.top.atoms if atoms.resid==0 if atoms.type=='ho']
    

    def fluorine_atomnames(self,resnumber):
        return ([atom.name for atom in self.top.atoms if atom.type=='f' if atom.resid==resnumber])

    def average_structure(self):
        """ Compute the Root Mean Square Displacement of each dendron using
            the starting point as reference structure """
        #return pt.mean_structure(self,mask=':'+str(dendron_id), frame_indices=None)
        #return pt.rmsf(self,mask=':'+str(dendron_id), frame_indices=None)
        results=[]
        for dendrons in range(1,self.ndendr+1):
            results.append(pt.rmsd(self,mask=':'+str(dendrons), frame_indices=None))
        with open('rmsd_dendrons_traj_'+str(self.index)+'.pickle','wb') as handle:
            pickle.dump(results,handle,protocol=pickle.HIGHEST_PROTOCOL)
    #### AVAILABLE ANALYSIS IMPLEMENTED FOR THIS CLASS 



    def radius_of_gyration(self):
        from compute_RDF import all_rdf
        results = all_rdf(self)
        with open('novolume_rdf_traj'+str(self.index)+'.pickle', 'wb') as handle:
            pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 
    def extract_dihedral(self):
        from compute_ANGLES import triazol_dihedral
        results = triazol_dihedral(self)
        np.savetxt('dihedrals_media_traj'+str(self.index)+'.txt',results[1],fmt='%.3f')
        dihedrals = np.transpose(results[0])
        np.savetxt('dihedrals_traj'+str(self.index)+'.txt', dihedrals ,fmt='%.3f')

    def extract_CCC_angle(self):
        from compute ANGLES import CCC_angle
        results = CCC_angle(self)
        np.savetxt('angleCCC_media_traj'+str(self.index)+'.txt',results[1],fmt='%.3f')
        angles = np.transpose(results[0])
        np.savetxt('angleCCC_traj'+str(self.index)+'.txt', angles ,fmt='%.3f')

    def extract_rdf(self, couples='all'):
        from compute_RDF import all_rdf
        if couples=='all':
            results=all_rdf(self)
            with open('novolume_rdf_traj'+str(self.index)+'.pickle', 'wb') as handle:
                pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 

    def extract_ff_distances(self):
        from compute_DISTANCES import FF_distances
        results = FF_distances(self)
        try:
            with open('FF_distances_traj'+str(self.index)+'.pickle','wb') as handle:
                pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
            np.savetxt('distances_FF_traj'+str(self.index)+'.txt', results)
        except ValueError:
            print('This method is available just for trajectory with more than one dendron')
        #with open('FF_distances_traj'+str(self.index)+'.pickle','wb') as handle:
        #    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
        return 

    def extract_c_ho_distances(self):
        from compute_DISTANCES import CH_distances
        results = CH_distances(self)
        with open('CH_distances_traj'+str(self.index)+'.pickle','wb') as handle:
            pickle.dump(results,handle,protocol=pickle.HIGHEST_PROTOCOL)
        return 

    ####### THIS METHOD PERFORMS ALL THE ANALYSIS LISTED ABOVE

    def all_analysis(self):
        import time
        start = time.time()
        self.radius_of_gyration()
        print('ROG analysis: Done')
        self.extract_dihedral()
        print('Dihedral analysis: Done')
        self.extract_rdf()
        print('RDF analysis: Done')
        self.extract_ff_distances()
        print('FF distances: Done')
        self.extract_c_ho_distances()
        print('C-ho distances: Done')
        end = time.time()
        print('Total time of execution: ', end-start ,'minutes.')
        return 

