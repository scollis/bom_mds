import sys
from os import getenv
sys.path.append(getenv('HOME')+'/bom_mds/modules')
from numpy import sqrt, array, argsort, abs, zeros, float, pi, linspace
from pylab import meshgrid, num2date, date2num

def diff_dates(d1,d2):
	return num2date(date2num(d1)-date2num(d2))

def date_linspace(d1,d2,N):
	return num2date(linspace(date2num(d1), date2num(d2), N))

