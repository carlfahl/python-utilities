#
# Utility functions
#
# Copyright 2011 Fahlstrom Research LLC
#
# Author : Carl A. Fahlstrom
#
# Date : 06-25-11
#

# List of all the command line flags that are supported for
# specifying printout of documentation.
#
help_li = ['--help', '-help', '-h', '--h', '-info', '--info']

def show_help(argv, doc_st):
	'''
	Function to print out the doc string passed to the function.
	After printing the help string the program exits.
	'''
	import sys

	help_li = ['--help', '-help', '-h', '--h', '-info', '--info']

	try:
		if sys.argv[1] in help_li:
			print doc_st
			sys.exit()
	except:
		print doc_st
		sys.exit()

def print_doc_string(meth):
	'''
	After printing the help string the program exits.
	'''
	import sys

	try:
		print meth.__doc__
		sys.exit()
	except:
		print "Failed to print doc string."
		sys.exit()

def find_closest_index(ar, val):
	'''
	'''

	accending = True

	if ar[0] > ar[1]:
		accending = False

	ret_index = 0	

	for elem in ar:
		if accending:
			if elem > val:
				ret_index = list(ar).index(elem)
				break
		else:
			if elem < val:
				ret_index = list(ar).index(elem)
				break

	return ret_index


def deriv(data):
	'''
	Returns the derivitive of input data based on midpoint
	method.
	'''

	li = []

	li.append(data[1] - data[0])

	for i in xrange(1,len(data)-1):
		li.append(data[i+1] - data[i-1])

	li.append(data[-1] - data[-2])

	return li


def select_cols(data, index_list, sep=None):
	'''
	select_cols(data, index_list, sep) - returns specified columns form a file.

	Input data is the output of the readlines function.

	index_list is a list of the column numbers desired (starts at 0).

	sep (optional) is the text character that seperates data; if not specified 
	seperation is at any whitespace.

	Copyright 2011 Fahlstrom Research LLC	

	Author : Carl A. Fahlstrom

	Date : 06-25-2011

	Version : 0.1
	'''

#	cdef int i

	li = []

	for elem in data:
		li2 = []
		tmp_li = elem.split(sep)
		for i in index_list:
			try:
				li2.append(tmp_li[i])
			except:
				print "Index not valid: ", i
				continue
		li.append(li2)

	return li

def convert_li_typ(in_list, out_type):
	'''
	'''

	for i in xrange(len(in_list)):
		in_list[i] = __builtins__[out_type](in_list[i])
	
	return in_list

def read_file(filename, split_lines=False, split_char=None):
	'''
	General purpose function for reading data from a text file.
	Includes and option to return a list of split data from a 
	text file of column data.
	'''
	
	try:
		batch_file = open(filename, 'r+')
		bfile_lines = batch_file.readlines()
		batch_file.close()
	except:
		print "Not a valid file"

	if split_lines:
		li = []
		for elem in bfile_lines:
			li.append(elem.split(split_char))
		return li
	else:
		return bfile_lines

def write_file(filename, lines):
	'''
	'''

	try:
		outfile = file(filename, 'w')
		outfile.write('\n'.join(lines))
		outfile.close()
	except:
		print 'Writing to File ', filename, ' Failed'

def write_file_2(filename, vec1, vec2):
        '''
        '''

        li = []

        for i in range(len(vec1)):
                li.append(str(vec1[i])+"\t"+str(vec2[i]))

        write_file(filename, li)

def isCommon(nest_list):
	'''
	'''
	list_length = len(nest_list)
	result = []
	for elem in nest_list[0]:
		for i in range(1,list_length):
			if elem not in nest_list[i]:
				break
			elif i == list_length-1:
				result.append(elem)
			else:
				continue
	return result

