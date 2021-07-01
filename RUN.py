#    The function returns the list containing all the path to possible trajectories : 
#    - 0: Two fdg3 in ETN 5%;            # instantiate as my_traj(0, 'ETN', 2)
#    - 1: One fdg3 in ETN 5% ;           # instantiate as my_traj(1, 'ETN', 1)
#    - 2: One fdg3 in WATER;             # instantiate as my_traj(2,  None, 1)
#    - 3: Two fdg3 in TFN (0-100 ns) ;   # instantiate as my_traj(3, 'TFN', 2)
#    - 4: Two fdg3 in TFN (100-200 ns) ; # instantiate as my_traj(4, 'TFN', 2)
#    - 5: Two fdg3 in TFN (0-200 ns);    # instantiate as my_traj(5, 'TFN', 2)
#    - 6: Two fdg3 in ETN (0-200 ns);    # instantiate as my_traj(6, 'ETN', 2)
#    - 7: One fdg3 in TFN 5%        ;    # instantiate as my_traj(7, 'TFN', 1)

import pytraj as pt
import matplotlib.pyplot as plt
from main import my_traj


mynewtraj=my_traj(0,'ETN', 2)

mynewtraj.all_analysis()
#mynewtraj.extract_c_ho_distances()

