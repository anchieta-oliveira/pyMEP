#!/usr/bin/env python3

# Description
###############################################################################

""" 
Module to store and manipulate data from PDB files.

This module contains the PDB class, which is used to store and manipulate
data from PDB files. The class has methods to read, write, and show data 
from PDB files.

Example:
--------

    >>> pdb = PDB()
    >>> pdb.read('example.pdb')
    >>> pdb.show()

"""

'''
PSF Files:
---------------------------------------
https://www.ks.uiuc.edu/Training/TutorialsOverview/namd/namd-tutorial-unix-html/node23.html
'''

# Imports
###############################################################################
import logging
import numpy as np



# License
###############################################################################
'''
pymep
Authors: -
[The Federal University of Rio de Janeiro]
Contact info:
Carlos Chagas Filho Institute of Biophysics
Laboratory for Molecular Modeling and Dynamics
Av. Carlos Chagas Filho 373 - CCS - bloco G1-19,
Cidade Universitária - Rio de Janeiro, RJ, CEP: 21941-902
E-mail address: -
This project is licensed under Creative Commons license (CC-BY-4.0) (Ver qual)
'''

# Classes
###############################################################################
class PSF:
	""" Class to store and manipulate data from PDB files.

	Attributes
	----------
	text : str
		PDB file content.

	"""

	def __init__(self, path: str = "") -> None:
		''' Constructor for PDB class.
		
		Parameters
		----------
		path : str, optional
			Path to the PDB file. Default is an empty string. If path is not provided the PDB object is created empty.
		'''

		self.text:str = ""
		self.data:np.array = np.array([])
		self.atoms:list = []
		self.resids:list = []
		self.name = "None"
		self.inid = np.array([], dtype=np.int32)
		self.resinid = np.array([], dtype=np.int32)
		self.serials = np.array([], dtype=np.int32)
		self.names = np.array([], dtype=np.str_)
		self.resnames = np.array([], dtype=np.str_)
		self.resseq = np.array([], dtype=np.int32)
		self.charges = np.array([], dtype=np.float32)
		if path != "":
			self.read(path)