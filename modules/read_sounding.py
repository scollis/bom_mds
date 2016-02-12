###################################################################
#read_sounding.py:                                                #
###################################################################
#Reads TWP-ICE style soundings                                    #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, JUNE 2008       #
###################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import *
from pylab import num2date, datestr2num, date2num
from numpy import interp, nan
from numpy import where as nwhere

def float_conv(strg):
	try:
		val=float(strg)
	except ValueError:
		val=nan
	return val


def parse_sounding_block(sounding_block):
	headers=sounding_block[0].split()[0:17]
	start_date_str=sounding_block[1].split()[0]+" "+sounding_block[1].split()[1]
	data_dict=dict([(headers[i],array([float_conv(sounding_block[j+1].split()[i]) for j in range(len(sounding_block)-1)])) for i in range(len(headers))])
	date_list=[num2date(datestr2num(sounding_block[i+1].split()[0]+" "+sounding_block[i+1].split()[1])) for i in range(len(sounding_block)-1)]
	data_dict.update({'date_list':array(date_list)})
	return data_dict



def read_sounding(fname):
	sonde_file=open(fname, 'r')
	sonde_ascii=sonde_file.readlines()
	sonde_file.close()
	header_pos=[]
	for line_number in range(len(sonde_ascii)):
		if sonde_ascii[line_number][0:4] == 'yyyy': header_pos.append(line_number)#print sonde_ascii[line_number]
	end_pos=array(header_pos[1:(len(header_pos))], dtype=integer)
	end_pos=append(end_pos, len(sonde_ascii)-1)
	sonde_list=[]
	for i in range(len(header_pos)):
		sounding_block=sonde_ascii[header_pos[i]:end_pos[i]]
		sonde_list=append(sonde_list, parse_sounding_block(sounding_block))
	return sonde_list

def read_sounding_within_a_day(fname, dateobj):
	sonde_file=open(fname, 'r')
	sonde_ascii=sonde_file.readlines()
	sonde_file.close()
	header_pos=[]
	for line_number in range(len(sonde_ascii)):
		if sonde_ascii[line_number][0:4] == 'yyyy': header_pos.append(line_number)#print sonde_ascii[line_number]
	end_pos=array(header_pos[1:(len(header_pos))], dtype=integer)
	end_pos=append(end_pos, len(sonde_ascii)-1)
	sonde_list=[]
	for i in range(len(header_pos)):
		sounding_block=sonde_ascii[header_pos[i]:end_pos[i]]
		date_of_block=datestr2num(sounding_block[1].split()[0]+" "+sounding_block[1].split()[1])
		if abs(date_of_block-date2num(dateobj)) < 1.0 : sonde_list=append(sonde_list, parse_sounding_block(sounding_block))
	return sonde_list


def interp_sounding(sonde_dict, levs):
	interp_vars=['dpcorr(deg)', 'alt(m)', 'lon(degs)', 'press(hPa)', 'lat(degs)', 'wspd(m/s)', 'tdry(degs)', 'mixr(g/kg)', 'dp(degs)', 'rh_tlag(%)', 'rh(%)', 'rh_corr(%)', 'wdir(degs)', 'zenith(deg)', 'mixrc(g/kg)']
	interp_dict=dict([(var,interp(levs, sonde_dict['alt(m)'], sonde_dict[var])) for var in set(interp_vars)&set(sonde_dict.keys())])
	return interp_dict


def get_two_best_sondes(date_str, **kwargs):
	sonde_file=kwargs.get('sonde_file', '/data/twpice/darwin.txt')
	#outdir=kwargs.get('outdir', '/flurry/home/scollis/bom_mds/dealias/')
	sonde_file=kwargs.get('sonde_file', '/data/twpice/darwin.txt')
	outdir=kwargs.get('outdir', '/home/scollis/bom_mds/dealias/')
	tim_date=num2date(datestr2num(date_str))
	sonde_list=read_sounding_within_a_day(sonde_file, tim_date)
	launch_dates=[sonde['date_list'][0] for sonde in sonde_list]
	#print launch_dates
	launch_date_offset=[date2num(sonde['date_list'][0])- date2num(tim_date)  for sonde in sonde_list]
	sonde_made_it=False
	candidate=0
	while not(sonde_made_it):
		best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[candidate]]
		candidate=candidate+1
		sonde_made_it=best_sonde['alt(m)'][-1] > 18000.
		if not sonde_made_it: print "Sonde Burst at ", best_sonde['alt(m)'][-1], "m rejecting"
	print "Sonde Burst at ", best_sonde['alt(m)'][-1], "m Accepting"
	sonde_made_it=False
	while not(sonde_made_it):
		sec_best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[candidate]]
		candidate=candidate+1
		sonde_made_it=sec_best_sonde['alt(m)'][-1] > 18000.
		if not sonde_made_it: print "Sonde Burst at ", sec_best_sonde['alt(m)'][-1], "m rejecting"
	print "Sonde Burst at ", sec_best_sonde['alt(m)'][-1], "m Accepting"
	print 'Time of radar: ', tim_date, ' Time of  best sonde_launch: ', best_sonde['date_list'][0], ' Time of sonde_termination: ', best_sonde['date_list'][-1]
	print 'Time of radar: ', tim_date, ' Time of second sonde_launch: ', sec_best_sonde['date_list'][0], ' Time of sonde_termination: ', best_sonde['date_list'][-1]
	for i in range(len(sonde_list)):
		best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[i]]
		print 'Time of radar: ', tim_date, ' Time of  best sonde_launch: ', best_sonde['date_list'][0], ' Offset', abs(date2num(best_sonde['date_list'][0])-date2num(tim_date))*24.0
	return best_sonde, sec_best_sonde



def is_sonde_good(sonde, min_ht, keylist):
	sonde_good=sonde['alt(m)'][-1] > min_ht
	if not(sonde_good): print "Sonde Burst at ", sonde['alt(m)'][-1], "m rejecting"
	for i in range(len(keylist)):
		occur=len(nwhere(array(sonde[keylist[i]])==-9999)[0])
		sonde_good=sonde_good and not(-9999 in sonde[keylist[i]])
		if (-9999 in sonde[keylist[i]]): print "Dirty data in ",keylist[i]," rejecting ", occur, " occurances"
		if occur < 5 and  (-9999 in sonde[keylist[i]]): print  "at ", sonde['alt(m)'][nwhere(array(sonde[keylist[i]])==-9999)[0]]
	return sonde_good



def get_two_best_conc_sondes(date_str, **kwargs):
	req_vars=kwargs.get('req_vars', [ 'alt(m)', 'press(hPa)',  'wspd(m/s)', 'tdry(degs)',  'wdir(degs)'])
	min_ht=kwargs.get('min_ht', 17000.)
	sonde_file=kwargs.get('sonde_file', '/data/twpice/darwin.txt')
	#outdir=kwargs.get('outdir', '/flurry/home/scollis/bom_mds/dealias/')
	sonde_file=kwargs.get('sonde_file', '/data/twpice/darwin.txt')
	outdir=kwargs.get('outdir', '/home/scollis/bom_mds/dealias/')
	tim_date=num2date(datestr2num(date_str))
	sonde_list=read_sounding_within_a_day(sonde_file, tim_date)
	launch_dates=[sonde['date_list'][0] for sonde in sonde_list]
	#print launch_dates
	launch_date_offset=[date2num(sonde['date_list'][0])- date2num(tim_date)  for sonde in sonde_list]
	sonde_made_it=False
	candidate=0
	while not(sonde_made_it):
		sonde_num=argsort(abs(array(launch_date_offset)))[candidate]
		best_sonde=sonde_list[sonde_num]
		candidate=candidate+1
		sonde_made_it=is_sonde_good(best_sonde, min_ht, req_vars)
		#sonde_made_it=best_sonde['alt(m)'][-1] > min_ht
		#if not sonde_made_it: print "Sonde Burst at ", best_sonde['alt(m)'][-1], "m rejecting"
	#print "Sonde Burst at ", best_sonde['alt(m)'][-1], "m Accepting"
	if launch_date_offset[sonde_num] < 0.:
		print "Best Sonde flight is before the radar scan"
		first_sonde=best_sonde
		candidate=sonde_num +1
		sonde_made_it=False
		while not(sonde_made_it):
			second_sonde=sonde_list[candidate]
			candidate=candidate+1
			sonde_made_it=is_sonde_good(second_sonde, min_ht, req_vars)
			#sonde_made_it=second_sonde['alt(m)'][-1] > min_ht
			#if not sonde_made_it: print "Sonde Burst at ", second_sonde['alt(m)'][-1], "m rejecting"
			if (candidate==len(sonde_list)) and not(sonde_made_it):
				print "I'm running off the end of the list, Bailing"
				second_sonde=first_sonde
				sonde_made_it=True
	elif launch_date_offset[sonde_num] > 0.:
		print "Best Sonde flight is after the radar scan"
		second_sonde=best_sonde
		candidate=sonde_num -1
		sonde_made_it=False
		while not(sonde_made_it):
			first_sonde=sonde_list[candidate]
			candidate=candidate-1
			sonde_made_it=is_sonde_good(first_sonde, min_ht, req_vars)
			#sonde_made_it=first_sonde['alt(m)'][-1] > min_ht
			#if not sonde_made_it: print "Sonde Burst at ", first_sonde['alt(m)'][-1], "m rejecting"
			if candidate==-1 and not(sonde_made_it):
				print "I'm running off the start of the list, Bailing"
				first_sonde=second_sonde
				sonde_made_it=True
	elif launch_date_offset[sonde_num]==0.:
		print "Sonde launched at exactly the same time of the scan... Huey be praised"
		first_sonde=second_sonde=best_sonde
	print 'Time of radar: ', tim_date
	print  ' Time of  first sonde_launch: ', first_sonde['date_list'][0], ' Time of first sonde_termination: ', first_sonde['date_list'][-1]
	print "Sonde Burst at ", first_sonde['alt(m)'][-1], "m Accepting"
	print  ' Time of second sonde_launch: ', second_sonde['date_list'][0], ' Time of second sonde_termination: ', second_sonde['date_list'][-1]
	print "Sonde Burst at ", second_sonde['alt(m)'][-1], "m Accepting"
	
	for i in range(len(sonde_list)):
		best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[i]]
		print 'Time of radar: ', tim_date, ' Time of  best sonde_launch: ', best_sonde['date_list'][0], ' Offset', abs(date2num(best_sonde['date_list'][0])-date2num(tim_date))*24.0
	
	return first_sonde, second_sonde


def interp_sonde_time(first_sonde, second_sonde, date_obj, levs):
	if hasattr(levs, 'sum'):#Numpy array
		levels=levs
	elif hasattr(levs, 'append'):#list
		levels=array(levs)
	else:#a number of levels to construct
		top_level=min(first_sonde['alt(m)'][-1],second_sonde['alt(m)'][-1])
		bottom_level=max(first_sonde['alt(m)'][0],second_sonde['alt(m)'][0])
		levels=linspace(bottom_level, top_level, levs)
	first_sonde_int=interp_sounding(first_sonde, levels)
	second_sonde_int=interp_sounding(second_sonde, levels)
	time_between_sondes=date2num(second_sonde['date_list'][0])-date2num(first_sonde['date_list'][0])
	if time_between_sondes==0.0:
		time_interp_sonde=first_sonde_int
	else:
		
		first_wt=(date2num(date_obj)-date2num(first_sonde['date_list'][0]))/time_between_sondes
		second_wt=(date2num(second_sonde['date_list'][0])- date2num(date_obj))/time_between_sondes
		print "First sonde weighting ", first_wt, " Second sonde weighting ", second_wt, " Sum (should be 1) ", first_wt+ second_wt
		time_interp_sonde=dict([(key,first_sonde_int[key]*first_wt+second_sonde_int[key]*second_wt) for key in first_sonde_int.keys()])
	return time_interp_sonde


def make_deal_sonde(date_str, **kwargs):
	#sonde_file=kwargs.get('sonde_file', '/bm/gdata/scollis/twpice/darwin.txt')
	outdir=kwargs.get('outdir', '/home/scollis/bom_mds/dealias/')
	tim_date=num2date(datestr2num(date_str))
	nl=kwargs.get('nl', 990)
	#sonde_list=read_sounding.read_sounding_within_a_day(sonde_file, tim_date)
	#launch_dates=[sonde['date_list'][0] for sonde in sonde_list]
	#launch_date_offset=[date2num(sonde['date_list'][0])- date2num(tim_date)  for sonde in sonde_list]
	#best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[0]]
	req=[ 'alt(m)',  'wspd(m/s)',  'wdir(degs)']
	first_sonde,second_sonde = get_two_best_conc_sondes(date_str, req_vars=req)
	interp_sonde=interp_sonde_time(first_sonde, second_sonde, tim_date, nl)
	#print 'Time of radar: ', tim_date, ' Time of sonde_launch: ', best_sonde['date_list'][0], ' Time of sonde_termination: ', best_sonde['date_list'][-1]
	#levels_onto=linspace(best_sonde['alt(m)'][0],  best_sonde['alt(m)'][-1], 900)
	#interp_sonde=read_sounding.interp_sounding(best_sonde, levels_onto)
	days=(int(date2num(tim_date))- datestr2num('01/01/'+date_str[0:4])) +1.0
	yr4="%(y)4d" %{'y':tim_date.year}
	date_dict={'ddd':days, 'HH':tim_date.hour, 'MM':tim_date.minute}
	post_name="%(ddd)03d%(HH)02d%(MM)02d_ti_sounding.txt" %date_dict
	dfname=yr4[2:4]+post_name
	fname=kwargs.get('fname', dfname)
	#days=(int(date2num(best_sonde['date_list'][0]))- datestr2num('01/01/'+date_str[0:4])) +1.0
	#yr4="%(y)4d" %{'y':best_sonde['date_list'][0].year}
	#date_dict={'ddd':days, 'HH':best_sonde['date_list'][0].hour, 'MM':best_sonde['date_list'][0].minute}
	#post_name="%(ddd)03d%(HH)02d%(MM)02d_sounding.txt" %date_dict
	#dfname=yr4[2:4]+post_name
	#fname=kwargs.get('fname', dfname)
	print "opening ", outdir+fname
	f=open(outdir+fname, 'w')
	f.write('H\n')
	for i in range(interp_sonde['alt(m)'].shape[0]):
		items=['alt(m)', 'press(hPa)', 'press(hPa)', 'wdir(degs)', 'wspd(m/s)', 'press(hPa)', 'press(hPa)', 'press(hPa)', 'press(hPa)']
		printdict=dict([(str(n+1),interp_sonde[items[n]][i]) for n in range(len(items))])
		prstr="%(1)f %(2)f %(3)f %(4)f %(5)f %(6)f %(7)f %(8)f %(9)f \n" %printdict
		f.write(prstr)
	f.write('here endeth the file \n')
	f.close()
	return fname

