# -*- coding: utf-8 -*-
###################################################################
#netcdf_utis.py: Utilities for saving radar data to NetCDF        #
###################################################################
#Module for bom_mds                                               #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
from Scientific.IO.NetCDF import *
from numpy import array
from pylab import date2num, num2date
def make_var(ncf, var_name, typecode, dim_tuple, vals):
	myvar=ncf.createVariable(var_name, typecode, dim_tuple)
	#print myvar.shape
	#print vals[50,50,0], ' ', vals[50,51,0]
	#vals=zeros(myvar.shape, dtype=float)
	for i in range(myvar.shape[0]):
		for j in range(myvar.shape[1]):
			for k in range(myvar.shape[2]):
				myvar[i,j,k]=float(vals[i,j,k])
	return myvar



def save_data_cube(radar1, radar2, ncf_fname):
	ncf=NetCDFFile(ncf_fname, 'w')
	ncf.createDimension('nx', len(radar1['xar']))
	ncf.createDimension('ny', len(radar1['yar']))
	ncf.createDimension('nl', len(radar1['levs']))
	ncf.createDimension('one', 1)
	ncf.createDimension('two', 2)
	avail_vars_radar1=set(radar1.keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'KD', 'i_comp', 'j_comp', 'k_comp', 'u_array', 'v_array', 'w_array']) 
	avail_vars_radar2=set(radar2.keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW','KD', 'i_comp', 'j_comp', 'k_comp']) 
	#testme=make_var(ncf, var, 'f',  ('nx', 'ny', 'nl'), gp[var])
	ncf_varlist_radar1=[make_var(ncf, var+'_radar1', 'f',  ('nx', 'ny', 'nl'), radar1[var]) for var in avail_vars_radar1]
	ncf_varlist_radar2=[make_var(ncf, var+'_radar2', 'f',  ('nx', 'ny', 'nl'), radar2[var]) for var in avail_vars_radar2]
	xvar=ncf.createVariable('xar', 'f', ('one', 'nx'))
	for i in range(len(radar1['xar'])): xvar[0,i]=float(radar1['xar'][i])
	yvar=ncf.createVariable('yar', 'f', ('one', 'ny'))
	for j in range(len(radar1['yar'])): yvar[0,j]=float(radar1['yar'][j])
	#yvar.assignValue(array([gp['yar']]))
	lvar=ncf.createVariable('levs', 'f', ('one', 'nl'))
	for k in range(len(radar1['levs'])): lvar[0,k]=float(radar1['levs'][k])
	#lvar.assignValue(array([gp['levs']]))
	rad1_locvar=ncf.createVariable('radar1_loc', 'f', ('one','two'))
	rad2_locvar=ncf.createVariable('radar2_loc', 'f', ('one','two'))
	rad1_disvar=ncf.createVariable('radar1_dis', 'f', ('one','two'))
	rad2_disvar=ncf.createVariable('radar2_dis', 'f', ('one','two'))
	rad1_datvar=ncf.createVariable('radar1_date', 'd', ('one','one'))
	rad2_datvar=ncf.createVariable('radar2_date', 'd', ('one','one'))
	rad1_locvar[0,0]=float(radar1['radar_loc'][0])
	rad1_locvar[0,1]=float(radar1['radar_loc'][1])
	rad2_locvar[0,0]=float(radar2['radar_loc'][0])
	rad2_locvar[0,1]=float(radar2['radar_loc'][1])
	rad1_disvar[0,0]=float(radar1['displacement'][0])
	rad1_disvar[0,1]=float(radar1['displacement'][1])
	rad2_disvar[0,0]=float(radar2['displacement'][0])
	rad2_disvar[0,1]=float(radar2['displacement'][1])
	rad1_datvar[0,0]=float(date2num(radar1['date']))
	rad2_datvar[0,0]=float(date2num(radar2['date']))
	setattr(ncf, 'radar1_name', radar1['radar_name'])
	setattr(ncf, 'radar2_name', radar2['radar_name'])
	
	#ncf_varlist_gp=dict([(var+'_gp',ncf.createVariable(var+'_gp', 'f', ('nx', 'ny', 'nl'))) for var in avail_vars_gp])
	#for var in avail_vars_gp:
	#	for i in range(ncf_varlist_gp[var+'_gp'].shape[0]):
	#		for j in range(ncf_varlist_gp[var+'_gp'].shape[1]):
	#			for k in range(ncf_varlist_gp[var+'_gp'].shape[2]):
	#				ncf_varlist_gp[var+'_gp'][i,j,k]=gp[var][i,j,k]
	ncf.close()
	

def load_cube(fname):#weird edits added for changed behaviour
	ncf=NetCDFFile(fname, 'r')
	print "Ok1"
	xar=ncf.variables['xar'][0]
	#print mxar.shape
	#xar=mxar[0]
	print "plew"
	yar=array(ncf.variables['yar'][0])
	print "Ok2"
	levs=array(ncf.variables['levs'][0])
	print "Ok3"
	parms=['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW','KD', 'i_comp', 'j_comp', 'k_comp', 'v_array', 'u_array', 'w_array']
	radar1_poss_parms=[par+'_radar1' for par in parms]
	radar2_poss_parms=[par+'_radar2' for par in parms]
	radar1_parms=set(radar1_poss_parms) & set(ncf.variables.keys())
	radar2_parms=set(radar2_poss_parms) & set(ncf.variables.keys())
	print "Doing radar1"
	radar1_dict=dict([(par[0:-7],array(ncf.variables[par].getValue())) for par in radar1_parms])
	radar1_dict.update({'xar':xar, 'yar':yar, 'levs':levs, 'radar_loc':ncf.variables['radar1_loc'][0], 'displacement':ncf.variables['radar1_dis'][0], 'date':num2date(ncf.variables['radar1_date'][0,0]), 'radar_name':getattr(ncf, 'radar1_name')})
	print "doing radar2"
	radar2_dict=dict([(par[0:-7],array(ncf.variables[par].getValue())) for par in radar2_parms])
	radar2_dict.update({'xar':xar, 'yar':yar, 'levs':levs, 'radar_loc':ncf.variables['radar2_loc'][0], 'displacement':ncf.variables['radar2_dis'][0], 'date':num2date(ncf.variables['radar2_date'][0,0]), 'radar_name':getattr(ncf, 'radar2_name')})
	ncf.close()
	return radar1_dict, radar2_dict







