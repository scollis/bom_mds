#################################################################
#Dealias a berrimah UF file using FourDD_berrmiah
#################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
#sys.path.append(os.getenv('HOME')+'/bom_mds')
#import ascii_scan_to_pickle
import read_sounding
from pylab import date2num, num2date, datestr2num
from numpy import array, argsort
from numpy import where as nwhere
from time import time as systime
#import ascii_scan_to_pickle

def dealias_arb(filename, format, path,deal_path, prefix):
	print "radar type: "+format
	if format=="lassen":
		fn=dealias_cpol_volume(filename, raw_path=path,deal_path=deal_path, pattern=prefix)
	elif format=="uf":
		fn=dealias_berrimah_volume(filename,  raw_path=path,deal_path=deal_path, pattern=prefix)
	return fn

def dealias_cpol_volume(filename, **kwargs):	
	pattern=kwargs.get('pattern', 'Gunn_pt_')
	deal_add=kwargs.get('deal_add', '_deal')
	raw_path=kwargs.get('raw_path', '/data/lassen_cpol/')
	deal_path=kwargs.get('deal_path', '/data/uf_cpol/')
	deal_files=os.listdir(deal_path)
	raw_files=os.listdir(raw_path)
        print deal_files
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
		print deal_str
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
	os.chdir('/home/scollis/bom_mds/dealias/')
	execbl='./FourDD_lassen '
	command=execbl+deal_path+deal_fname+' '+raw_path+filename+' '+deal_path+outfile+' '+ sonde_name+' '+deal_date+' '+'0 1 1 1'
	print command
	os.system(command)
	os.chdir(cwd)	
	return outfile

#/bm/gscratch/scollis/uf_ber/BerrimaVol20060117_035004.uf
def dealias_berrimah_volume(filename, **kwargs):
	pattern=kwargs.get('pattern', 'BerrimaVol')
	deal_add=kwargs.get('deal_add', '_deal')
	raw_path=kwargs.get('raw_path', '/data/uf_ber/')
	deal_path=kwargs.get('deal_path', '/data/deal_ber/')
	
	deal_files=os.listdir(deal_path)
	raw_files=os.listdir(raw_path)
	#deal_files=[]
	#raw_files=[]
	#for file in uf_files:
	#	if 'deal' in file: 
	#		deal_files.append(file)
	#	else:
	#		raw_files.append(file)
	#
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
	os.chdir('/home/scollis/bom_mds/dealias/')
	execbl='./FourDD_berrimah '
	command=execbl+deal_path+deal_fname+' '+raw_path+filename+' '+deal_path+outfile+' '+ sonde_name+' '+deal_date+' '+'0 1 1 1'
	print command
	os.system(command)
	os.chdir(cwd)
	return outfile
