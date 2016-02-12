#Utilities to read in an ascii file as prouced by a Lassen to ascii converter 
#Scott Collis
#CAWCR April 2008

#Initial idea is to open the whole file and then dump each radial scan to file and use the  built in tool to read and
#parse that file, except that idea sucked... I found a better way using stringobj.split()
import sys
sys.path.append('/flurry/home/scollis/pylibs/lib64/python2.4/site-packages')
from numpy import *
from scipy import *
#import matplotlib
from pylab import num2date, datestr2num
#sys.path.append('/home/scott/python/text_formatter/modules/')
#sys.path.append('/home/scott/fortran/read_cappi/')
sys.path.append('/flurry/home/scollis/python/radar/')
from radar_math import *
import propigation
import mathematics
#from mathematics import *
#from read_cappi import read_cappi
#from matplotlib.toolkits.basemap import Basemap #map plotting things for matplotlib
from os import listdir, getenv

def isvalue(instr):
	striped_str1=instr.replace('.','')
	striped_str=striped_str1.replace('-','')
	return striped_str.isdigit()



def stripz(nstr):
	slen=len(nstr)
	for i in range(slen-1):
		if nstr[0]=='0': nstr=nstr[1:slen-i]
		#print nstr
	return nstr

def evalz(nstr):
	#print nstr
	if isvalue(nstr):
		rv=eval(stripz(nstr))
	else:
		rv='ER'
	
	return rv


def parse_info_line(info_line):
	#date_obj=num2date(datestr2num(info_line[(info_line.find('date/time=')+len('date/time=')):-1]))
	#n_rays=int(info_line[(info_line.find('number_of_rays=')+len('number_of_rays=')):(info_line.find('date/time=')-1)])
	ilsplit=info_line.split()
	date_obj=num2date(datestr2num(ilsplit[-2]+' '+ilsplit[-1]))
	sweep_number=int(ilsplit[ilsplit.index("sweep_no=")+1])
	nrays=int(ilsplit[ilsplit.index("number_of_rays=")+1])
	return {'date':date_obj, 'sweep number':sweep_number, 'rays':nrays}

def parse_info_line_uf(info_line, ascii_list):
	#date_obj=num2date(datestr2num(info_line[(info_line.find('date/time=')+len('date/time=')):-1]))
	#n_rays=int(info_line[(info_line.find('number_of_rays=')+len('number_of_rays=')):(info_line.find('date/time=')-1)])
	ilsplit=info_line.split()
	if ilsplit[6][0]=='2':
		d1=ilsplit[6][2:4]
		d2=ilsplit[6][5:7]
		d3=ilsplit[6][8:10]
	else:
		d1=ilsplit[6][0:2]
		d2=ilsplit[6][3:5]
		d3=ilsplit[6][6:8]
	date_obj=num2date(datestr2num(d3+'/'+d2+'/'+d1 + ' '+ilsplit[7]))
	sweep_number=int(ilsplit[ilsplit.index("Scan")+2])
	ngates=int(ilsplit[ilsplit.index("ngates")+2])
	nrays=(len(ascii_list)-1)/(ngates+2+1) #Number of lines - one blank line at top / 511
	return {'date':date_obj, 'sweep number':sweep_number, 'rays':nrays}



def split_cpol_rays(info_line, ascii_list):
	#bound it
	nrays=parse_info_line(info_line)['rays']
	data=ascii_list[3:-2]
	block_size=len(data)/nrays
	blocks=[data[(n*(block_size+1)):(block_size+n*(block_size+1))] for n in range(nrays)]
	return blocks


def split_uf_rays(info_line, ascii_list):
	#bound it
	nrays=parse_info_line_uf(info_line, ascii_list)['rays']
	data=ascii_list[2:-3]
	block_size=len(data)/nrays
	blocks=[data[(n*(block_size+1)):(block_size+n*(block_size+1))] for n in range(nrays)]
	return blocks


def parse_block(block):
	preamble=dict([(block[0].split()[2*n].replace('=',''), float(block[0].split()[2*n+1])) for n in range(4)])
	headers=block[1].split()
	data_dict=dict([(headers[i],array([float(block[j+2].split()[i]) for j in range(len(block)-2)])) for i in range(len(headers))])
	data_dict.update(preamble)
	return data_dict

def parse_block_uf(block):
	ilsplit=block[0].split()
	extract=['ngates', 'Azmth', 'Elev']
	preamble=dict([(item, float(ilsplit[ilsplit.index(item)+2])) for item in extract])
	headers=block[1].split()
	data_dict=dict([(headers[i],array([float(block[j+2].split()[i]) for j in range(len(block)-2)])) for i in range(len(headers))])
	data_dict.update(preamble)
	return data_dict


def read_ufascii(**kwargs):
	radar_name=kwargs.get('radar_name', 'Berrimah')
	fname=kwargs.get('fname', '/bm/gdata/scollis/berrimah/20060122_035004/Berrimah_20060122_035004_01.ascii')
	sweep_file=open(fname, 'r')
	sweep_ascii=sweep_file.readlines()
	sweep_file.close()
	sweep_ascii.append('\n')
	sweep_ascii.append('\n')
	sweep_ascii.append('\n')
	info_line=sweep_ascii[2]
	info_dict=parse_info_line_uf(info_line, sweep_ascii)
	blocks=split_uf_rays(info_line, sweep_ascii)
	#assuming all rays have the same number of gates
	dict0=parse_block_uf(blocks[0])
	gated_vars=[ 'VR', 'CZ', 'VE', 'TI', 'range', 'ZD', 'SW', 'VR', 'PH', 'RH']
	sweep_dict=dict([(key,zeros([len(blocks),len(dict0[key])], dtype=float)) for key in set(dict0.keys())&set(gated_vars) ])
	sweep_dict.update(dict([(key,zeros([len(blocks)],dtype=float)) for key in set(dict0.keys())-set(gated_vars) ]))
	ray_num=1
	for block in blocks:
		#print ray_num
		#print block[-1]
		bdict=parse_block_uf(block)
		for key in set(bdict.keys())&set(gated_vars):
			#print key
			sweep_dict[key][ray_num-1, :]=bdict[key]
		for key in set(bdict.keys())-set(gated_vars):
			#print key
			sweep_dict[key][ray_num-1]=bdict[key]
		ray_num=ray_num+1
	sweep_dict.update(info_dict)
	sweep_dict.update({'radar_name':radar_name})
	return sweep_dict



def read_cpol(**kwargs):
	radar_name=kwargs.get('radar_name', 'C-POL')
	fname=kwargs.get('fname', '/data/cpol_ppi/CPOL_20060122_1540_01.ascii')
	sweep_file=open(fname, 'r')
	sweep_ascii=sweep_file.readlines()
	sweep_file.close()
	info_line=sweep_ascii[1]
	info_dict=parse_info_line(info_line)
	blocks=split_cpol_rays(info_line, sweep_ascii)
	#assuming all rays have the same number of gates
	dict0=parse_block(blocks[0])
	gated_vars=[ 'UZ', 'ZD', 'SW', 'VR', 'CZ', 'range', 'PH', 'RH']
	sweep_dict=dict([(key,zeros([len(blocks),len(dict0[key])], dtype=float)) for key in set(dict0.keys())&set(gated_vars) ])
	sweep_dict.update(dict([(key,zeros([len(blocks)],dtype=float)) for key in set(dict0.keys())-set(gated_vars) ]))
	for block in blocks:
		bdict=parse_block(block)
		ray_num=bdict['ray_no']
		for key in set(bdict.keys())&set(gated_vars):
			sweep_dict[key][ray_num-1, :]=bdict[key]
		for key in set(bdict.keys())-set(gated_vars):
			sweep_dict[key][ray_num-1]=bdict[key]
	sweep_dict.update(info_dict)
	sweep_dict.update({'radar_name':radar_name})
	return sweep_dict


def append_xyz(sweep_dict):
	xar=zeros(sweep_dict['range'].shape, dtype=float)
	yar=zeros(sweep_dict['range'].shape, dtype=float)
	zar=zeros(sweep_dict['range'].shape, dtype=float)
	for i in range(sweep_dict['range'].shape[0]):
		for j in range(sweep_dict['range'].shape[1]):
			x,y,z=radar_coords_to_cart(sweep_dict['range'][i,j],sweep_dict['Azmth'][i] ,sweep_dict['Elev'][i] , debug=False)
			xar[i,j]=x
			yar[i,j]=y
			zar[i,j]=z
	sweep_dict.update({'xar':xar, 'yar':yar,'zar':zar})
	return sweep_dict

#Re=6371.0*1000.0
#h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
#s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
#	p_r=4.0*Re/3.0
#	rm=rng*1000.0
#	z=(rm**2 + p_r**2 + 2.0*rm*p_r*sin(ele*pi/180.0))**0.5 -p_r
#arc length
#	s=p_r*arcsin(rm*cos(ele*pi/180.)/(p_r+z))
#	if debug: print "Z=", z, "s=", s
#	y=s*cos(az*pi/180.0)
#	x=s*sin(az*pi/180.0)
	
#def unit_vector(x,y,h,debug=False, **kwargs):
def append_ijk(sweep_dict):
	#[i,j] i is azmth, j is gates
	iar=zeros(sweep_dict['range'].shape, dtype=float)
	jar=zeros(sweep_dict['range'].shape, dtype=float)
	kar=zeros(sweep_dict['range'].shape, dtype=float)
	Re=6371.0*1000.0
	p_r=4.0*Re/3.0
	for az_num in range(sweep_dict['range'].shape[0]):
		r=sweep_dict['range'][az_num,:]*1000.0
		ele=sweep_dict['Elev'][az_num]
		az=sweep_dict['Azmth'][az_num]
		z=(r**2 + p_r**2 + 2.0*r*p_r*sin(ele*pi/180.0))**0.5 -p_r
		s=p_r*arcsin(r*cos(ele*pi/180.)/(p_r+z))
		dh=mathematics.dy(z)
		ds=mathematics.dy(s)
		angle=arctan(dh/ds)
		k_comp=sin(angle)
		r_comp=cos(angle)
		i_comp=r_comp*sin(az*pi/180.0)
		j_comp=r_comp*cos(az*pi/180.0)
		for gate in range(sweep_dict['range'].shape[1]):
			if sweep_dict['VR'][az_num, gate]==-999: i_comp[gate]=j_comp[gate]=k_comp[gate]=-999.
		#print iar.shape
		#print i_comp.shape
		iar[az_num,:]=i_comp
		jar[az_num,:]=j_comp
		kar[az_num,:]=k_comp
	sweep_dict.update({'i_comp':iar, 'j_comp':jar,'k_comp':kar})
	return sweep_dict


def construct_lassen_scan(**kwargs):
	path=kwargs.get('path','/media/iriver/data/cpol_ppi/20060122_0357/')
	files=listdir(path)
	files.sort()
	scan_list=[]
	for file in files:
		print "Loading File ", files.index(file)+1, " of ", len(files), "files"  
		scan_list.append(append_ijk(append_xyz(read_cpol(fname=path+file))))
		
	return scan_list


def construct_uf_scan(**kwargs):
	psth={}
	if "radar_name" in kwargs.keys(): psth.update({"radar_name":kwargs["radar_name"]}) 
	path=kwargs.get('path','/bm/gdata/scollis/berrimah/20060122_035004/')
	files=listdir(path)
	files.sort()
	scan_list=[]
	for file in files:
		print "Loading File ", files.index(file)+1, " of ", len(files), "files"  
		scan_list.append(append_ijk(append_xyz(read_ufascii(fname=path+file, **psth))))
	return scan_list


#def plot_ppi(sweep_dict, **kwargs):
#	f=figure()
#	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
#	fig_name=kwargs.get('fig_name', 'cappi_output.png')
#	Re=6371.0*1000.0
#	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
#	lons=radar_loc[1]+360.0*sweep_dict['xar']/(rad_at_radar*2.0*pi)
#	lats=radar_loc[0] + 360.0*sweep_dict['yar']/(Re*2.0*pi)
#	map= Basemap(projection='merc',lat_0=lats.mean(), lon_0=lons.mean(),llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max() , urcrnrlon=lons.max(), resolution='l',area_thresh=1., lat_ts=lats.mean())
#	xx, yy = map(lons, lats)
	#map.drawcoastlines()
	#map.drawcountries()
#	map.drawmapboundary()
#	map.readshapefile('/flurry/home/scollis/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
#	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
#	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
#	map.contourf(xx,yy,sweep_dict['CZ'], levels=linspace(10,62,14))
#	title('dBz')
#	colorbar()
#	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
#	close(f)




#def plot_ppi_vel(sweep_dict, **kwargs):
#	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
#	fig_name=kwargs.get('fig_name', 'cappi_output.png')
#	f=figure()
#	#gp_loc=[-12.2492,  131.0444]
#	Re=6371.0*1000.0
#	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
#	lons=radar_loc[1]+360.0*sweep_dict['xar']/(rad_at_radar*2.0*pi)
#	lats=radar_loc[0] + 360.0*sweep_dict['yar']/(Re*2.0*pi)
#	map= Basemap(projection='merc',lat_0=lats.mean(), lon_0=lons.mean(),llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max() , urcrnrlon=lons.max(), resolution='l',area_thresh=1., lat_ts=lats.mean())
#	xx, yy = map(lons, lats)
	#map.drawcoastlines()
	#map.drawcountries()
#	map.drawmapboundary()
#	map.readshapefile('/flurry/home/scollis/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
#	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
#	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
#	map.contourf(xx,yy,sweep_dict['VR'], levels=linspace(-15,15,31))
#	title('dBz')
#	colorbar()
#	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
#	close(f)



def scan_list_to_ldata(slist, parm):
	xar=[]
	yar=[]
	zar=[]
	par=[]
	#todo: put in a dh/dr element for the k-velocity component
	print "hi"
	for scn in slist:
		for i in range(scn['xar'].shape[0]):
			for j in range(array(scn['xar']).shape[1]):
				xar.append(scn['xar'][i,j])
				yar.append(scn['yar'][i,j])
				zar.append(scn['zar'][i,j])
				par.append(scn[parm][i,j])
	return xar, yar, zar, par

def plot_cappi(xx,yy,data, **kwargs):
	print "making a Cafe-Cappi!"
	print xx.shape
	print yy.shape
	print data.shape
	f=figure()
	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
	fig_name=kwargs.get('fig_name', 'cappi_output.png')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*xx/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*yy/(Re*2.0*pi)
	print lats.shape
	print lons.shape
	map= Basemap(projection='merc',lat_0=lats.mean(), lon_0=lons.mean(),llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max() , urcrnrlon=lons.max(), resolution='l',area_thresh=1., lat_ts=lats.mean())
	longr, latgr=meshgrid(lons,lats)
	print longr.shape
	print latgr.shape
	xg, yg = map(longr, latgr)
	print xg.shape
	print yg.shape
	#map.drawcoastlines()
	#map.drawcountries()
	map.drawmapboundary()
	map.readshapefile('/flurry/home/scollis/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	map.contourf(xg,yg,data, levels=linspace(10,62,14))
	title('dBz')
	colorbar()
	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
	close(f)

def plot_cappi_vel(xx,yy,data, **kwargs):
	print "making a Cafe-Cappi!"
	print xx.shape
	print yy.shape
	print data.shape
	f=figure()
	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
	fig_name=kwargs.get('fig_name', 'cappi_output.png')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*xx/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*yy/(Re*2.0*pi)
	print lats.shape
	print lons.shape
	map= Basemap(projection='merc',lat_0=lats.mean(), lon_0=lons.mean(),llcrnrlat=lats.min(), llcrnrlon=lons.min(), urcrnrlat=lats.max() , urcrnrlon=lons.max(), resolution='l',area_thresh=1., lat_ts=lats.mean())
	longr, latgr=meshgrid(lons,lats)
	print longr.shape
	print latgr.shape
	xg, yg = map(longr, latgr)
	print xg.shape
	print yg.shape
	#map.drawcoastlines()
	#map.drawcountries()
	map.drawmapboundary()
	map.readshapefile('/flurry/home/scollis/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	map.contourf(xg,yg,data, levels=linspace(-15,15,31))
	title('dBz')
	colorbar()
	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
	close(f)



if __name__ == "__main__":
	print "hi"
	scanme=construct_lassen_scan(path='/bm/gdata/scollis/20060122_0357/')
	plot_ppi(scanme[5])
	xar, yar, zar, par=scan_list_to_ldata(scanme, 'CZ')
	xara=array(xar)
	yara=array(yar)
	zara=array(zar)
	para=array(par)
	para[where(para <0)]=0
	xxx=linspace(-150,150,75)*1000.0
	yyy=linspace(-150,150,75)*1000.0
	#def make_cappi(xa,ya,za,pa,xs, ys, z, npts=2):
	use=where(abs(zara-2500.0) < 1500.0)
	ca=transpose(make_cappi(xara[use],yara[use],zara[use],para[use],xxx,yyy,2.5*1000.0))
	f=figure()
	print xxx.shape
	print yyy.shape
	print ca.shape
	contourf(xxx,yyy,ca, levels=linspace(10,62,14))
	colorbar()
	savefig('/flurry/home/scollis/cap.png')
	close(f)
	plot_cappi(xxx,yyy,ca)










#array([float(blocks[0][i+2].split()[0]) for i in range(len(blocks[0])-2)]