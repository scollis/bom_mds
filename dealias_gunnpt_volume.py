#################################################################
#Dealias a Gunn Pt  Lassen file using FourDD_lassen
#################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
sys.path.append(os.getenv('HOME')+'/bom_mds')
import ascii_scan_to_pickle
import read_sounding
from pylab import date2num, num2date, datestr2num
from numpy import array, argsort
from numpy import where as nwhere
from time import time as systime
import parse_ini

def dealias_single_volume(filename, **kwargs):	
	pattern=kwargs.get('pattern', 'Gunn_pt_')
	deal_add=kwargs.get('deal_add', '_deal')
	lassen_path=kwargs.get('lassen_path', '/bm/gscratch/scollis/lassen_cpol/')
	uf_path=kwargs.get('lassen_path', '/bm/gscratch/scollis/uf_cpol/')
	deal_files=os.listdir(uf_path)
	raw_files=os.listdir(lassen_path)
	#deal_files=[]
	#raw_files=[]
	#for file in uf_files:
	#	if 'deal' in file: 
	#		deal_files.append(file)
	#	else:
	#		raw_files.append(file)
	#
	if not(filename in raw_files):
		raise IOError, 'File not there: '+filename
	deal_files.sort()
	raw_files.sort()
	this_date_str=filename[len(pattern):-11]
	#this_date_str=this_date_str[0:this_date_str.find(' ')+5]
	sonde_name=read_sounding.make_deal_sonde(this_date_str)
	deal_date=sonde_name[0:sonde_name.find('_')]
	date_num_list=[]
	for file in deal_files:
		deal_str=file[len(pattern):-8]
		#print deal_str
		#deal_str=deal_str[0:deal_str.find(' ')+5]
		#print deal_str
		date_num_list.append(datestr2num(deal_str))
	offset=array(date_num_list)-datestr2num(this_date_str)
	idec_where=nwhere(offset < 0.0)[0]
	offset_least=100.0
	if len(idec_where)!=0: offset_least=offset[idec_where[-1]]
	print "**********************"
	if abs(offset_least) < 30.0/(60.0*24.0):
		deal_fname=deal_files[idec_where[-1]]
		print "Found a previous de-aliased file ", deal_fname
	else:
		print "No de-aliased file found within 30 minutes, only using sounding"
		deal_fname='dummy'
	print "********************"
	outfile=filename[0:-11]+deal_add+".uf"
	cwd=os.getcwd()
	os.chdir('/flurry/home/scollis/bom_mds/dealias/')
	execbl='./FourDD_lassen '
	command=execbl+uf_path+deal_fname+' '+lassen_path+filename+' '+uf_path+outfile+' '+ sonde_name+' '+deal_date+' '+'1 1 1 1'
	print command
	os.system(command)
	os.chdir(cwd)	
	return outfile

def dealias_multi(first, last, **kwargs):
	pattern=kwargs.get('pattern', 'Gunn_pt_')
	deal_add=kwargs.get('deal_add', '_deal')
	lassen_path=kwargs.get('lassen_path', '/bm/gscratch/scollis/lassen_cpol/')
	uf_path=kwargs.get('lassen_path', '/bm/gscratch/scollis/deal_cpol/')
	deal_files=os.listdir(uf_path)
	raw_files=os.listdir(lassen_path)
	if not(first in raw_files) or not(last in raw_files):
		raise IOError, 'File not there'
	deal_files.sort()
	raw_files.sort()
	fi=raw_files.index(first)
	li=raw_files.index(last)
	sublist=raw_files[fi:li]
	cwd=os.getcwd()
	for file in sublist:
		outfile=dealias_single_volume(file)
		bits=outfile[len(pattern):-3]
		com2="./cpol_uf_to_ascii.sh "+outfile
		os.chdir("/flurry/home/scollis/bom_mds/uf2ascii/")
		os.system(com2)
		com3="python /flurry/home/scollis/bom_mds/ascii_scan_to_pickle.py " +bits
		#os.system(com3)
		ascii_scan_to_pickle.cpoluf_ascii_to_pgz(bits+'/', debug=True)
	os.chdir(cwd)



#dealias_single_volume('BerrimaVol20060117_054004_deal.uf')
#dealias_single_volume('BerrimaVol20060117_004004.uf')
if __name__ == "__main__":
	t0=systime()
	print "DEALIASING!"
	print sys.argv
	#save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	#simple_reconstruction_3d_pytest(sys.argv[1],sys.argv[2], sys.argv[3])
	#recon(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	#test_pert_winds()
	#test_sounding()
	dealias_multi(sys.argv[1], sys.argv[2])
	print "Finished running runtime=",systime()-t0, "Seconds"

