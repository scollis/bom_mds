###################################################################
#save_scan.py                                                     #
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
#import shelve
import pickle_zip

def pzip_radar(radar, **kwargs):
	debug=kwargs.get('debug', False)
	date_obj=radar[0]['date']
	date_dic={'y':date_obj.year, 'm':date_obj.month, 'd':date_obj.day, 'H':date_obj.hour, 'M':date_obj.minute}
	dfname=radar[0]['radar_name']+"_%(y)04d%(m)02d%(d)02d_%(H)02d%(M)02d.pickle.gz" %date_dic
	fname=kwargs.get('fname', dfname)
	path=kwargs.get('path', '/bm/gdata/scollis/radar_shelves/')
	#dbm=shelve.open(path+fname)
	#dbm['scan']=radar
	#dbm.close
	if debug: print "Saving Gzipped Pickle file to "+path+fname
	pickle_zip.save(radar, path+fname)

