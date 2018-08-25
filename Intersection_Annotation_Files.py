import itertools
file_list=["D4_0.clean.variant.txt", "D4_1.clean.variant.txt", "D5_0.clean.variant.txt", "D5_1.clean.variant.txt"]


def getHeader(sample_input_list):
	input_file_name=file_list[0]
	input_file = open(input_file_name, 'r')
	fileheader = input_file.readline().strip()
	return fileheader

def input_sample_data_to_key_value_list(file_list):
	print("Input sample data to value list \n")
	all_dict_key_value_array=[]
	all_key_set_array=[]
	file_series_number=0
	for input_file_name in file_list:
		input_file=open(input_file_name,'r')
		sample_dataset=[str(calling_line).strip().split("\t") for calling_line in input_file if len(calling_line.strip())!=0]
		sample_dict_key_value=dict()

		not_First_line=0
		for select_values in sample_dataset:
			if not_First_line!=0:
				sample_dict_key_value[select_values[0]+"\t"+select_values[1]+"\t"+select_values[2]+"\t"+select_values[4]]="\t".join(select_values)
			not_First_line=not_First_line+1

		all_dict_key_value_array.append(sample_dict_key_value)
		all_dict_key_array=[]
		all_dict_key_array.append(sample_dict_key_value.keys())
		########### Print how many calls in samples ##############
		print(input_file_name+": has "+str(len(sample_dict_key_value)))+" calls"
		##########################################################
		file_series_number=file_series_number+1
	input_file.close()
	#print(all_dict_key_value_array)
	return all_dict_key_value_array

def convert_list_to_set_list(sample_input_list):
	print("Convert list to set \n")
	all_set_list=[]
	for each_sample in sample_input_list:
		sample_set=set()
		for each_value in each_sample:
			sample_set.add(each_value)
		all_set_list.append(sample_set)
	return all_set_list

def double_layer_list_to_set(sample_input_list):
	list_to_set=set()
	for first_layer_set in sample_input_list:
		temp_set=set()
		for second_layer_value in first_layer_set:
			temp_set.add(second_layer_value)
		value_to_frozenset=frozenset(temp_set)
		list_to_set.add(value_to_frozenset)
	return list_to_set

def construct_combination(sample_input_list):
	print("Construct combinations \n")
	num_list=range(len(sample_input_list))
	combinations_list=[]
	for set_size in range(1, len(num_list)+1):
		combinations_single_layer_list=[list(item) for item in itertools.combinations(num_list, set_size)]
		combinations_list.append(combinations_single_layer_list)
	return num_list, combinations_list

def intersection(comb_list,all_set_list):
	print("Process intersection \n")
	process_set_result_dict={}
	for set_size in range(len(comb_list)):
		for each_set_list in comb_list[set_size]:
			process_list_set=[]
			process_complement_list_set=[]
			complement_set_list=list(set(num_list).difference(set(each_set_list)))
			############ Print each set of select region and complement region #########
			print("Set:"+str(each_set_list))
			for each_sample_in_set in each_set_list:
				process_list_set.append(all_set_list[each_sample_in_set])
			for each_sample_in_complement_set in complement_set_list:
				process_complement_list_set.append(all_set_list[each_sample_in_complement_set])
			intersection_result=set.intersection(*process_list_set)
			union_result=set.union(*process_list_set)
			if len(process_complement_list_set) >= 1:
				complement_union_result=set.union(*process_complement_list_set)
			else:
				complement_union_result=set()
			purely_specific_region=intersection_result-complement_union_result
			print("Union region calling numbers is :"+str(len(union_result)))
			print("Purely specific region calling numbers is :"+str(len(purely_specific_region))+"\n")
			process_set_result_dict[tuple(each_set_list)]=purely_specific_region
	return process_set_result_dict

def lookup_whole_value_from_key(set_key_dict,key_whole_value_dict_array):
	region_calling_value_dict=dict()
	print("Lookup whole values from key \n")
	for purely_specific_region, calling_key_set in set_key_dict.items():
		print(purely_specific_region)
		calling_value_set=set()
		for calling_key in calling_key_set:
			if_has_replace=0
			for sample_number in purely_specific_region:
				if if_has_replace==0:
					calling_value=key_whole_value_dict_array[sample_number][calling_key]
					calling_value_set.add(calling_value)
					if_has_replace=if_has_replace+1
		region_calling_value_dict[purely_specific_region]=calling_value_set
	print("\n\n")
	return region_calling_value_dict

def write_dict_to_files(input_dict):
	print("Write dict to files \n")
	for purely_specific_region, calling_value_set in input_dict.items():
		output_filename=str(purely_specific_region)+".txt"
		print(str(purely_specific_region)+" has "+str(len(calling_value_set))+" counts.")
		outputfile=open(output_filename,'w')
		outputfile.write(fileheader+"\n")
		for calling_value in calling_value_set:
			outputfile.write(str(calling_value)+"\n")

def sortedDictValues(input_dict): 
	keys = inpu_dict.keys()
	keys.sort()
	return [dict[key] for key in keys] 

fileheader=getHeader(file_list)
all_dict_key_value_array=input_sample_data_to_key_value_list(file_list)
all_set_list=convert_list_to_set_list(all_dict_key_value_array)
num_list, comb_list=construct_combination(file_list)
process_set_result_dict=intersection(comb_list,all_set_list)
region_calling_value_dict=lookup_whole_value_from_key(process_set_result_dict,all_dict_key_value_array)
write_dict_to_files(region_calling_value_dict)
