###################################################################
#dump_pickle.py#
#################################################################################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
import propigation	
import read_rays
import radar_math
import radar_to_cart
import mathematics
#import grad_conj_solver_plus
import read_sounding
import met
from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random
import pickle_zip


def dump_pickle(pickle_name):
	locs={'Berrimah':[-12.457, 130.925], 'C-POL':[-12.2492,  131.0444], 'Berrimah_deal':[-12.457, 130.925], 'C-POL_deal':[-12.2492,  131.0444]}
	basedir='/bm/gkeep/scollis/deal_ber/'
	radar=pickle_zip.load(basedir+pickle_name)
	pres.dump_radar(radar, locs[radar[0]['radar_name']])

if __name__ == "__main__":
	t0=systime()
	print "the uber cool test"
	print sys.argv
	#save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	#simple_reconstruction_3d_pytest(sys.argv[1],sys.argv[2], sys.argv[3])
	#recon(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	#test_pert_winds()
	#test_sounding()
	dump_pickle(sys.argv[1])
	print "Finished running runtime=",systime()-t0, "Seconds"

