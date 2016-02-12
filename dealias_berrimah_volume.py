#################################################################
#Dealias a berrimah UF file using FourDD_berrmiah
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
#import ascii_scan_to_pickle

#/bm/gscratch/scollis/uf_ber/BerrimaVol20060117_035004.uf
def dealias_single_volume(filename, **kwargs):
	pattern=kwargs.get('pattern', 'BerrimaVol')
	deal_add=kwargs.get('deal_add', '_deal')
	uf_path=kwargs.get('uf_path', '/bm/gscratch/scollis/uf_ber/')
	uf_files=os.listdir(uf_path)
	
	deal_files=[]
	raw_files=[]
	for file in uf_files:
		if 'deal' in file: 
			deal_files.append(file)
		else:
			raw_files.append(file)
	
	if not(filename in raw_files):
		raise IOError, 'File not there'
	deal_files.sort()
	raw_files.sort()
	this_date_str=filename[len(pattern):-3].replace('_', ' ')
	this_date_str=this_date_str[0:this_date_str.find(' ')+5]
	sonde_name=read_sounding.make_deal_sonde(this_date_str)
	deal_date=sonde_name[0:sonde_name.find('_')]
	
	date_num_list=[]
	for file in deal_files:
		deal_str=file[len(pattern):-8].replace('_', ' ')
		deal_str=deal_str[0:deal_str.find(' ')+5]
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
	outfile=filename[0:-3]+deal_add+".uf"
	cwd=os.getcwd()
	os.chdir('/flurry/home/scollis/bom_mds/dealias/')
	execbl='./FourDD_berrimah '
	command=execbl+uf_path+deal_fname+' '+uf_path+filename+' '+uf_path+outfile+' '+ sonde_name+' '+deal_date+' '+'1 1 1 1'
	print command
	os.system(command)
	os.chdir(cwd)
	return outfile

def dealias_multi(first, last, **kwargs):
	uf_path=kwargs.get('uf_path', '/bm/gscratch/scollis/uf_ber/')
	uf_files=os.listdir(uf_path)
	raw_files=[]
	deal_files=[]
	for file in uf_files:
		if 'deal' in file: 
			deal_files.append(file)
		else:
			raw_files.append(file)
	if not(first in raw_files) or not(last in raw_files):
		raise IOError, 'File not there'
	deal_files.sort()
	raw_files.sort()
	fi=raw_files.index(first)
	li=raw_files.index(last)
	sublist=raw_files[fi:li]
	pattern=kwargs.get('pattern', 'BerrimaVol')
	cwd=os.getcwd()
	for file in sublist:
		outfile=dealias_single_volume(file)
		bits=outfile[len(pattern):-3]
		com2="./uf_to_ascii.sh "+outfile
		os.chdir("/flurry/home/scollis/bom_mds/uf2ascii/ber/")
		os.system(com2)
		com3="python /flurry/home/scollis/bom_mds/ascii_scan_to_pickle.py " +bits
		#os.system(com3)
		ascii_scan_to_pickle.uf_ascii_to_pgz(bits+'/', debug=True)
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

