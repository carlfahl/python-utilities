#!/usr/bin/python

#
# File : v2_geo.py
#
# Produces a geoemtry relavent to the calculation of electronic coupling constants of indole to backbone
# amide (produces n or c terminal amide) from a protein. 
#
# 04/20/2013 - added the ability to process the case of NATA
#
# 04/22/2013 - added special case of PRO
#
# Copyright 2011-2013 Fahlstrom Research LLC
#
# Author : Carl A. Fahlstrom
#
#

import sys
from utils_fr import read_file, write_xyz, process_opts
from utils_fr_vector import new_vec

def v2_geo(infile, outfile, n_ter, resnum):
	'''
	v2_geo -

	Carl A. Fahlstrom - 09Aug06
	'''

	# read the pdb file
	#
	pro_coord =  read_file(infile, True)

	# Ring atoms and other atoms of Trp that will always be included.
	#
	ring_ats = ['CA', 'HA', 'CB', 'HB1', 'HB2', '1HB', '2HB', 'CG', 'CD1', 'HD1', 'CD2', 'NE1', 'HE1', 'CE2', \
		'CE3', 'HE3', 'CZ2', 'HZ2', 'CZ3', 'HZ3', 'CH2', 'HH2']

	# Other amino acids; ACE and NH2 added for the case of NATA and similar
	#
	other_aa = ['GLY', 'ALA', 'SER', 'CYS', 'THR', 'ILE', 'VAL', 'MET', 'ASP', 'ASN', 'LEU', 'LYS', \
		'GLU', 'GLN', 'PRO', 'ARG', 'HIS', 'PHE', 'TYR', 'ACE', 'NH2']

	x_trp = []
	y_trp = []
	z_trp = []
	atoms = []

	# If the first residue is the TRP to be extracted
	# it must be n-terminal
	#
	if (pro_coord[5][3] == 'TRP'):
		n_ter = True

	x_sa = 0.0	
	y_sa = 0.0	
	z_sa = 0.0	

	x_sc = 0.0
	y_sc = 0.0
	z_sc = 0.0

	no_rep_o = False

	for elem in pro_coord:
		try:
			if elem[0] not in ['ATOM', 'HETATM']:
				continue
			if ((elem[3] == 'TRP') and (int(elem[4]) == resnum)):
				if elem[2] in ring_ats:
					x_trp.append(float(elem[5]))
					y_trp.append(float(elem[6]))
					z_trp.append(float(elem[7]))
					atoms.append(elem[2][0])
					if elem[2] == 'CA':
						x_sa = float(elem[5])
						y_sa = float(elem[6])
						z_sa = float(elem[7])
				if n_ter:
					if (elem[2] in ['C', 'O']):
						x_trp.append(float(elem[5]))
						y_trp.append(float(elem[6]))
						z_trp.append(float(elem[7]))
						atoms.append(elem[2][0])
					if elem[2] == 'N':
						x_rep_t = float(elem[5])
						y_rep_t = float(elem[6])
						z_rep_t = float(elem[7])
				else:
					if (elem[2] in ['N', 'H']):
						x_trp.append(float(elem[5]))
						y_trp.append(float(elem[6]))
						z_trp.append(float(elem[7]))
						atoms.append(elem[2][0])
					if elem[2] == 'C':
						x_rep_t = float(elem[5])
						y_rep_t = float(elem[6])
						z_rep_t = float(elem[7])
	
			elif (elem[3] in other_aa):
				# Special conditions for NATA 
				#
				if ((elem[3] == 'ACE') and (n_ter == True)):
					continue
				if ((elem[3] == 'NH2') and (n_ter == True)):
					x_trp.append(float(elem[5]))
					y_trp.append(float(elem[6]))
					z_trp.append(float(elem[7]))
					atoms.append(elem[2][0])
					no_rep_o = True
					continue
				# CH3 is for NATA and similar, special case for NATA ACE, TRP, NH2 all have
				# resnum of 1
				#
				if ((elem[2] in ['CH3']) and (elem[3] == 'ACE') and (n_ter == False)):
					x_rep_o = float(elem[5])
					y_rep_o = float(elem[6])
					z_rep_o = float(elem[7])
				if (n_ter and (int(elem[4]) == resnum+1)) or ((not n_ter) and (int(elem[4]) == resnum-1)):
					if elem[2] in ['CA']:
						x_rep_o = float(elem[5])
						y_rep_o = float(elem[6])
						z_rep_o = float(elem[7])
				# H1 is for NATA and similar
				#
				if (n_ter and (int(elem[4]) == resnum+1)):
					if (elem[2] in ['N', 'H', 'H1']):
						x_trp.append(float(elem[5]))
						y_trp.append(float(elem[6]))
						z_trp.append(float(elem[7]))
						atoms.append(elem[2][0])
						if elem[2] == 'N':
							x_sc = float(elem[5])
							y_sc = float(elem[6])
							z_sc = float(elem[7])
					# Special Case for Pro
					#
					if (elem[3] == 'PRO' and elem[2] == 'CD'):
						x_rep_p = float(elem[5])
						y_rep_p = float(elem[6])
						z_rep_p = float(elem[7])
				# Needed to add conditional for ACE for the case of NATA
				# 
				elif (((not n_ter) and (int(elem[4]) == resnum-1)) or ((not n_ter) and elem[3] == 'ACE')):
					if (elem[2] in ['C', 'O']):
						x_trp.append(float(elem[5]))
						y_trp.append(float(elem[6]))
						z_trp.append(float(elem[7]))
						atoms.append(elem[2][0])
						if elem[2] == 'C':
							x_sc = float(elem[5])
							y_sc = float(elem[6])
							z_sc = float(elem[7])
			else:
				break

		except:
			continue

	# Replace an atom with H
	#
	atoms.append('H')
	new_c = new_vec([x_rep_t,x_sa], [y_rep_t,y_sa], [z_rep_t,z_sa])
	x_trp.append(new_c[1])
	y_trp.append(new_c[2])
	z_trp.append(new_c[3])

	# Replace an atom with H if not NATA
	#
	if not no_rep_o:
		atoms.append('H')
		new_c = new_vec([x_rep_o,x_sc], [y_rep_o,y_sc], [z_rep_o,z_sc])
		x_trp.append(new_c[1])
		y_trp.append(new_c[2])
		z_trp.append(new_c[3])

	# Special case for Pro
	#
	try:
		new_c = new_vec([x_rep_p,x_sc], [y_rep_p,y_sc], [z_rep_p,z_sc])
		atoms.append('H')
		x_trp.append(new_c[1])
		y_trp.append(new_c[2])
		z_trp.append(new_c[3])
	except:
		pass

	li = []
	for i in xrange(len(atoms)):
		li.append(atoms[i]+'\t'+ str(x_trp[i])+'\t'+str(y_trp[i])+'\t'+str(z_trp[i]))

	write_xyz(outfile, li)

if __name__ == "__main__":
	if sys.argv[1] == "--help":
		print v2_geo.__doc__
	else:
		st2bool = {'True':True, 'False':False}
		opts = {'-if':None, '-of':None, '-nter':'False', '-resnum':1}
		opts = process_opts(opts, sys.argv[1:])
		v2_geo(opts['-if'], opts['-of'], st2bool[opts['-nter']], int(opts['-resnum']))

