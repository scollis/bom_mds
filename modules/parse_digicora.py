import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import *
from pylab import num2date, datestr2num, date2num
from numpy import interp, nan
from numpy import where as nwhere
import climt

def float_conv(strg):
	try:
		val=float(strg)
	except ValueError:
		val=nan
	return val

def read_digicora(fname):
#fname='/bm/gdata/scollis/digicora/Brisbane_Ap_200802081100.digicora'
	sonde_file=open(fname, 'r')
	sonde_ascii=sonde_file.readlines()
	sonde_file.close()
	buf=[]
	#easiest just to prescan
	for item in sonde_ascii:
		if item[0] !='[': 
			buf.append(item)
		elif item!='[$]\n':
			hflag=item=='[Header]\n' #we are reading in the header at the moment
			sonde_flag=item=='[Trace]\n' #we are reading the trace
			buf=[] #reset the buffer as we are at the start of a new data stream
		elif item=='[$]\n':
			if sonde_flag: sonde_items=copy(buf)
			if hflag: header_items=copy(buf)
	titles=(sonde_items[0].replace(',',' ')).split()
	data_dict=dict([(titles[i],array([float_conv((sonde_items[j+1].replace(',',' ')).split()[i]) for j in arange(len(sonde_items)-2)+1])) for i in range(len(titles))])
	for item in header_items:
		pair=((item.replace(' ','_')).replace('=',' ')).split()
		try:
			name=pair[0][0:(pair[0].index('['))]
		except ValueError:
			name=pair[0]
		try:
			data_dict.update({name:float(pair[1])})
		except ValueError:
			strg=pair[1].replace('"', '')
			if (('date' in name.lower()) or ('time' in name.lower())):
				try:
					data_dict.update({name:num2date(datestr2num(strg))})
				except ValueError:
					data_dict.update({name:strg})
			else:
				data_dict.update({name:strg})
	return data_dict



def append_z(sonde, key='alt(m)', pkey='Pres', tkey='Temp'):
	sonde.update({key:climt.thermodyn.z(sonde[pkey][::-1]*100.0, sonde[tkey][::-1]+273.0)[::-1]})
	return sonde
 

def digi_to_twpice(digi):
	new_sonde={}
	pairs={'aifstime':'date_list', 'AbsHum':'AbsHum', 'Temp':'tdry(degs)', 'Geop':'Geop', 'w_units':'w_units', 'Dewp':'dp(degs)', 'stationnumber':'stationnumber', 'levels':'levels', 'Pres':'press(hPa)', 'stationname':'stationname', 'fields':'fields', 'Spd':'wspd(m/s)', 'Dir':'wdir(degs)', 'alt(m)':'alt(m)'}
	for key in digi.keys():
		if key=='aifstime':
			new_sonde.update({pairs[key]:[digi[key]]*len(digi['Pres'])})
		else:
			new_sonde.update({pairs[key]:digi[key]})
	return new_sonde

#conv_sonde={'press(hPa)':data_dict['Pres'], 'dp(degs)':data_dict['Dewp'], 'tdry(degs)':data_dict['Temp']}
#pres.plot_sonde(conv_sonde, 'BrisbaneAP200802081100.png')
#kwargs={}
#sonde=data_dict

def fix_sounding(sonde, **kwargs):
	new_sonde={}
	level_index=kwargs.get('level_index','Pres')
	for item in sonde.keys():
		try:
			glev=[]
			gdata=[]
			for i in range(len(sonde[item])):
				if sonde[item][i]!=-9999.0:
					glev.append(sonde[level_index][i])
					gdata.append(sonde[item][i])
			if len(glev) > 3:
				new_sonde.update({item:interp(sonde[level_index][-1] - sonde[level_index], sonde[level_index][-1] - array(glev), array(gdata))})
		except TypeError: #no length to object
			new_sonde.update({item:sonde[item]})
	return new_sonde

#conv_sonde={'press(hPa)':new_sonde['Pres'], 'dp(degs)':new_sonde['Dewp'], 'tdry(degs)':new_sonde['Temp']}
