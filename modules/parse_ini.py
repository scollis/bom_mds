#Parse_ini: Text file parameter parsing program
#Scott Collis
#CAWCR June 2008
__author__ = "Scott Collis"
__version__ = "1.0"

def parse_ini(filename, **kwargs):
	com_str=kwargs.get('com_str','#')
	ini_file=open(filename, 'r')
	ini_ascii=ini_file.readlines()
	ini_file.close()
	ini_dict={}
	for line in ini_ascii:
		if not(line[0]==com_str):
			item_list=line.split()
			if len(item_list)==2:
				try:
					ini_dict.update({item_list[0]:float(item_list[1])})
				except ValueError:
					ini_dict.update({item_list[0]:item_list[1]})
			elif  len(item_list)>2:
				items=[]
				for item in item_list[1:len(item_list)]:
					try: 
						items.append(float(item))
					except ValueError:
						items.append(item)
				ini_dict.update({item_list[0]: items})
	return ini_dict
