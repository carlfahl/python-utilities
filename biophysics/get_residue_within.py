#!/usr/bin/python

#
# File :
#
# Copyright 2013 Fahlstrom Research LLC
#
# Author : Carl A. Fahlstrom
#

import sys
import string
from utils_fr import read_file, write_file
from utils_fr_vector import new_vec, unit_vec
from GuiBase import GuiBase

class GetResWithin:
	'''
	'''
	
	def __init__(self):
		'''
		'''

		self.res_to_get = []
		self.hetatm_res = []
		self.dist_stand =  0.0
		self.write_list = []
		self.atoms = {}
		self.resid = {}
		self.resname = {}
		self.chainid = {}
		self.coords_x = {}
		self.coords_y = {}
		self.coords_z = {}
		self.cs = {}
		self.cas = {}
		self.ns = {}

	def set_dist(self, dist):
		'''
		'''

		self.dist_stand = dist

	def setup_inputs(self, infile):
		'''
		'''

		pro_pdb_file = read_file(infile, True)

		for line in pro_pdb_file:
			if line[0] in ['ATOM','HETATM']:
				self.resid[int(line[1])] = int(line[5])
				self.resname[int(line[5])] = line[3]
				self.chainid[int(line[1])] = line[4]
				self.atoms[int(line[1])] = line[2]
				self.coords_x[int(line[1])] = float(line[6])
				self.coords_y[int(line[1])] = float(line[7])
				self.coords_z[int(line[1])] = float(line[8])
				if line[0] == 'ATOM':
					if line[2] == 'C':
						self.cs[int(line[5])] = int(line[1])
					if line[2] == 'N':
						self.ns[int(line[5])] = int(line[1])
					if line[2] == 'CA':
						self.cas[int(line[5])] = int(line[1])

	def get_residues_within(self, infile, center_at, dist_stand, st_line, end_line, out_file, debug):
		'''
		get_residues_within -

		Carl A. Fahlstrom - 09Aug06
		'''

		pro_pdb_file = read_file(infile, True)

		for line in pro_pdb_file:
			if line[0] in ['ATOM','HETATM']:
				self.resid[int(line[1])] = int(line[5])
				self.resname[int(line[5])] = line[3]
				self.chainid[int(line[1])] = line[4]
				self.atoms[int(line[1])] = line[2]
				self.coords_x[int(line[1])] = float(line[6])
				self.coords_y[int(line[1])] = float(line[7])
				self.coords_z[int(line[1])] = float(line[8])
				if line[0] == 'ATOM':
					if line[2] == 'C':
						self.cs[int(line[5])] = int(line[1])
					if line[2] == 'N':
						self.ns[int(line[5])] = int(line[1])
					if line[2] == 'CA':
						self.cas[int(line[5])] = int(line[1])

		if type(center_at) is type([]):
			x0 = []
			y0 = []
			z0 = []
			for elem in center_at:
				x0.append(self.coords_x[elem])
				y0.append(self.coords_y[elem])
				z0.append(self.coords_z[elem])
		else:
			x0 = self.coords_x[center_at]
			y0 = self.coords_y[center_at]
			z0 = self.coords_z[center_at]

		for i,j in self.coords_x.items():
			if type(center_at) is type([]):
				for k in range(len(x0)):
					dist = (((x0[k] - j)**2 + (y0[k] - self.coords_y[i])**2 + (z0[k] - self.coords_z[i])**2)**0.5)
					if dist < dist_stand:
						if self.resid[i] not in self.res_to_get:
							self.res_to_get.append(self.resid[i])
			else:
				dist = (((x0 - j)**2 + (y0 - self.coords_y[i])**2 + (z0 - self.coords_z[i])**2)**0.5)
				if dist < dist_stand:
					if self.resid[i] not in self.res_to_get:
						self.res_to_get.append(self.resid[i])

		res_to_cap_2 = [] # list of residues that require caping of both n and c terminus
		res_to_cap_n = [] # list of residues that require caping of both n and c terminus
		res_to_cap_c = [] # list of residues that require caping of both n and c terminus
		for elem in self.res_to_get:
			print self.resname[elem], "  ", elem
			if elem not in [494,496]:
				if (elem-1 not in self.res_to_get) and (elem+1 not in self.res_to_get):
#					res_to_cap_2.append(elem)
					self.cap_res_2(elem)
				elif elem-1 not in self.res_to_get:
#					res_to_cap_n.append(elem)
					self.cap_res_2(elem, True, False)
				elif elem+1 not in self.res_to_get:
#					res_to_cap_c.append(elem)
					self.cap_res_2(elem, False, True)

		self.get_residues(self.res_to_get)

		self.write_output(out_file)

	def write_output(self, out_file):
		'''
		'''

		num_ats = len(self.write_list)
		self.write_list.insert(0," ")
		self.write_list.insert(0," "+str(num_ats))

		write_file(out_file, self.write_list)

	def cap_res_2(self, resnum, cap_n=True, cap_c=True):
		'''
		'''

		if cap_c:
			# Add an Oxygen to the C terminal end of a lone amino acid
			vec1 = unit_vec([self.coords_x[self.ns[resnum+1]], self.coords_x[self.cs[resnum]]], \
				[self.coords_y[self.ns[resnum+1]], self.coords_y[self.cs[resnum]]], \
				[self.coords_z[self.ns[resnum+1]], self.coords_z[self.cs[resnum]]])

			O_coord_x = vec1[1] * 1.38 + self.coords_x[self.cs[resnum]]
			O_coord_y = vec1[2] * 1.38 + self.coords_y[self.cs[resnum]]
			O_coord_z = vec1[3] * 1.38 + self.coords_z[self.cs[resnum]]

			self.write_list.append("O\t" + str(O_coord_x) + "\t" + str(O_coord_y) + "\t" + str(O_coord_z))

			vec2 = unit_vec([self.coords_x[self.cas[resnum+1]], self.coords_x[self.ns[resnum]]], \
				[self.coords_y[self.cas[resnum+1]], self.coords_y[self.ns[resnum]]], \
				[self.coords_z[self.cas[resnum+1]], self.coords_z[self.ns[resnum]]])

			H_coord_x = vec1[1] * 1.38 + self.coords_x[self.cs[resnum]] + vec2[1] * 0.98
			H_coord_y = vec1[2] * 1.38 + self.coords_y[self.cs[resnum]] + vec2[2] * 0.98
			H_coord_z = vec1[3] * 1.38 + self.coords_z[self.cs[resnum]] + vec2[3] * 0.98

			self.write_list.append("H\t" + str(H_coord_x) + "\t" + str(H_coord_y) + "\t" + str(H_coord_z))

		if cap_n:
			vec3 = new_vec([self.coords_x[self.cs[resnum-1]], self.coords_x[self.ns[resnum]]], \
        	                [self.coords_y[self.cs[resnum-1]], self.coords_y[self.ns[resnum]]], \
        	                [self.coords_z[self.cs[resnum-1]], self.coords_z[self.ns[resnum]]], 1.01)

			self.write_list.append("H\t" + str(vec3[1]) + "\t" + str(vec3[2]) + "\t" + str(vec3[3]))

	def get_residues(self, resnum):
		'''
		'''

		for i,j in self.coords_x.items():
			if self.resid[i] in resnum:
				k = 0
				while self.atoms[i][k] not in string.uppercase:
					k += 1
				self.write_list.append(self.atoms[i][k] + "\t" + str(j) + "\t" + str(self.coords_y[i]) + "\t" + str(self.coords_z[i]))

	def add_atom_center(self, atomlist):
		'''
		'''

		vecs = []
		for elem in atomlist[1:]:
			vec.append([self.coords_x[atomlist[0]] - self.coords_x[elem], \
				self.coords_y[atomlist[0]] - self.coords_y[elem], \
				self.coords_z[atomlist[0]] - self.coords_z[elem]])

class GetResGui(GuiBase):
	'''
	GUI under development
	'''

	def __init(self):
		'''
		'''

		self.w = gtk.window(GTK_TOPLEVEL)

if __name__ == "__main__":
	if sys.argv[1] == "--help":
		print add_hydrogens.__doc__
        else:
		x = GetResWithin()
		if sys.argv[1] == '-gui':
			y = GetResGui()
		# opts
		# ================
		# -of : output file
		# -if : Input file
		# -ca : center atom(s)
		# -car : center atoms range
		# -sl : 
		# -el : 
		# -dist : distance parameter
                opts = {'-of':None, '-if':None, '-res':0, '-ca':0, '-car':None, '-sl':0, '-el':0, '-dist':0, '-debug':False}
		for i in range(len(sys.argv[1:])):
			if sys.argv[1+i] in opts.keys():
				elem = sys.argv[1+i]
                                idex = i+2
				if sys.argv[idex][-1] == ',':
					opts[elem] = []
				else:
                                	opts[elem] = sys.argv[idex]
				while sys.argv[idex][-1] == ',':
					if int(sys.argv[idex][:-1]) not in opts[elem]:
						opts[elem].append(int(sys.argv[idex][:-1]))
					idex += 1
					if sys.argv[idex][-1] == ',':
						opts[elem].append(int(sys.argv[idex][:-1]))
					else:
						opts[elem].append(int(sys.argv[idex]))
		if type(opts['-ca']) is type(''):
			opts['-ca'] = int(opts['-ca'])
		if type(opts['-res']) is type(''):
			opts['-res'] = [int(opts['-res'])]
#		x.get_residues_within(opts['-if'], opts['-ca'], float(opts['-dist']), int(opts['-sl']), int(opts['-el']), opts['-of'], opts['-debug'])
		x.setup_inputs(opts['-if'])
		x.get_residues(opts['-res'])
#		for elem in opts['-res']:
#			x.cap_res_2(elem)
		x.write_output(opts['-of'])

