###################################################################
#ascii_scan_to_pickle.py                                          #
###################################################################
#Module for bom_mds                                               #
#Save a scan variable and GZIP it and GUNZIPING and loading       #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from matplotlib import use as use_interface
# chose a non-GUI backend
use_interface( 'Agg' )
import propigation	
import read_rays
import radar_math
import radar_to_cart
import read_sounding
import met
from time import time as systime
import save_scan

def lass_ascii_to_pgz(dirname, **kwargs):
	base_path=kwargs.get('base_path','/bm/gscratch/scollis/gunn_pt/') 
	debug=kwargs.get('debug', False)
	if debug: print "Loading "+base_path+dirname
	lass_scan=read_rays.construct_lassen_scan(path=base_path+dirname)
	save_scan.pzip_radar(lass_scan, debug=debug)

def cpoluf_ascii_to_pgz(dirname, **kwargs):
	base_path=kwargs.get('base_path','/bm/gscratch/scollis/cpol/') 
	debug=kwargs.get('debug', False)
	if debug: print "Loading "+base_path+dirname
	uf_scan=read_rays.construct_uf_scan(path=base_path+dirname, radar_name='C-POL')
	if 'deal' in dirname: uf_scan[0]['radar_name']=uf_scan[0]['radar_name']+'_deal'
	save_scan.pzip_radar(uf_scan, debug=debug)


def uf_ascii_to_pgz(dirname, **kwargs):
	base_path=kwargs.get('base_path','/bm/gscratch/scollis/berrimah/') 
	debug=kwargs.get('debug', False)
	if debug: print "Loading "+base_path+dirname
	uf_scan=read_rays.construct_uf_scan(path=base_path+dirname)
	if 'deal' in dirname: uf_scan[0]['radar_name']=uf_scan[0]['radar_name']+'_deal'
	save_scan.pzip_radar(uf_scan, debug=debug)

if __name__ == "__main__":
	if len(sys.argv)==5:
		ber_path='/bm/gscratch/scollis/berrimah/'
		files_ber=os.listdir(ber_path)
		files_ber.sort()
		st_ber=files_ber.index(sys.argv[1])
		end_ber=files_ber.index(sys.argv[2])
		print files_ber[st_ber:end_ber]
		for dirname in files_ber[st_ber:end_ber]: uf_ascii_to_pgz(dirname+'/', debug=True)
		
		gp_path='/bm/gscratch/scollis/gunn_pt/'
		files_gp=os.listdir(gp_path)
		files_gp.sort()
		st_gp=files_gp.index(sys.argv[3])
		end_gp=files_gp.index(sys.argv[4])
		print files_gp[st_gp:end_gp]
		for dirname in files_gp[st_gp:end_gp]: lass_ascii_to_pgz(dirname+'/', debug=True)
	elif len(sys.argv)==2:
		dirname=sys.argv[1]
		if '_' in dirname:
			uf_ascii_to_pgz(dirname+'/', debug=True)
		else:
			cpoluf_ascii_to_pgz(dirname+'/', debug=True)
			#lass_ascii_to_pgz(dirname+'/', debug=True)