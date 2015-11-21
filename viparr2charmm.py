########################################################################
# This script converts viparr3 parameters for small molecule ligands
# into a CHARMM-streamfile. Created by sam.hertig@gmail.com.
#
# USAGE:
#
# python thisscriptname.py input/viparrdir output/charmmdir
#
# Optionally, placing a stream file obtained from paramchem in the 
# viparr3 directory will change atom names and types.
#
# The input directory has to be in Viparr3 format. Use 
# "ConvertViparr1.py" to convert a directory from Viparr1 to 
# Viparr3 format.
#
########################################################################

from os.path import join
import os, optparse
from itertools import product
from decimal import Decimal as dec
import networkx as nx
from networkx.algorithms import isomorphism

# Create datastructures found in "templates"-file:
class Mass:
	""" These go into the .rtf section with lines starting with 'MASS'. """
	def __init__(self, atomtype, atomnumber, atommass, memo=None):
		self.type = atomtype
		self.nr = atomnumber
		self.mass = atommass
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '	
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				atomtype=typedict[self.type]
			except KeyError:
				#atomtype = self.type
				print 'An extra mass parameter involving atom type %s was found, but the corresponding atom seem to be missing. Skipping.' %(self.type)
				return
		else:
			atomtype = self.type	
		filehandle.write( 'MASS %d %-7s %s %s \n' %(self.nr, atomtype, self.mass, self.memo) )
			
class Atom:
	""" These go into the .rtf section with lines starting with 'ATOM'. """
	def __init__(self, atomname, atomtype, atomcharge, element, memo=None):
		self.name = atomname 
		self.type = atomtype
		self.charge = atomcharge
		self.element =  element # only used if atomtypes need to be converted using  paramchem file
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '
	def convertelement(self):
		""" Converts atom number (integer) to element (string). """
		elements = {
			'1': 'H',
			'6': 'C',
			'7': 'N',
			'8': 'O',
			'9': 'F',
			'15': 'P',
			'13':'Al', 
			'16': 'S',
			'17': 'Cl',
			'35': 'Br',
			'53': 'I',
		}
		self.element = elements[self.element]
	def writeto(self, filehandle, namedict=None, typedict=None):
		if typedict:
			atomtype = typedict[self.type]
			atomname = namedict[self.name]
		else:
			atomtype = self.type
			atomname = self.name
		filehandle.write( 'ATOM %-7s%-7s %s %s \n' %(atomname, atomtype, self.charge, self.memo) )

class Bond:
	""" These go into the .rtf section with lines starting with 'BOND'.
	Assumes that these have no memo/comment. """
	def __init__(self, bondedatom1, bondedatom2):
		self.atom1 = bondedatom1 
		self.atom2 = bondedatom2
	def writeto(self, filehandle, namedict=None):
		if namedict:
			atom1 = namedict[self.atom1]
			atom2 = namedict[self.atom2]
		else:
			atom1 = self.atom1
			atom2 = self.atom2	
		filehandle.write( 'BOND %7s %7s \n' %(atom1, atom2) )			

class Imp:
	""" These go into the .rtf section with lines starting with 'IMPR'. 
	Assumes that these have no memo/comment. """
	def __init__(self, impatom1, impatom2, impatom3, impatom4):
		self.atom1 = impatom1 
		self.atom2 = impatom2
		self.atom3 = impatom3 
		self.atom4 = impatom4	
	def writeto(self, filehandle, namedict=None):
		if namedict:
			atom1 = namedict[self.atom1]
			atom2 = namedict[self.atom2]
			atom3 = namedict[self.atom3]
			atom4 = namedict[self.atom4]			
		else:
			atom1 = self.atom1
			atom2 = self.atom2
			atom3 = self.atom3
			atom4 = self.atom4				
		filehandle.write( 'IMPR %7s %7s %7s %7s \n' %(atom1, atom2, atom3, atom4) )		


# Create datastructures found in other files:

class BondParams:
	""" Bonds are specified in "stretch_harm". """
	def __init__(self, atomtype1, atomtype2, bondparam1, bondparam2, memo=None):
		self.type1 = atomtype1
		self.type2 = atomtype2
		self.param1 = bondparam1
		self.param2 = bondparam2
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '			
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				type1 = typedict[self.type1]
				type2 = typedict[self.type2]
			except KeyError:
				#type1 = self.type1
				#type2 = self.type2
				print 'An extra bond parameter involving atom types %s %s was found, but the corresponding atoms seem to be missing. Skipping.' %(self.type1, self.type2)
				return
		else:
			type1 = self.type1
			type2 = self.type2	
		filehandle.write( '%-7s%-7s %s  %s  %s \n' %(type1, type2, self.param2, self.param1, self.memo) )	

class AngleParams:
	""" Angles are specified in "angle_harm", and UB terms are found in "ureybradley_harm". """
	def __init__(self, atomtype1, atomtype2, atomtype3, angleparam1, angleparam2, memo=None):
		self.type1 = atomtype1
		self.type2 = atomtype2
		self.type3 = atomtype3
		self.param1 = angleparam1
		self.param2 = angleparam2
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '			
		self.extraparam1 = None # Urey Bradley term
		self.extraparam2 = None # Urey Bradley term
		self.extramemo = None # Urey Bradley memo string
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				type1 = typedict[self.type1]
				type2 = typedict[self.type2]
				type3 = typedict[self.type3]
			except KeyError:
				#type1 = self.type1
				#type2 = self.type2	
				#type3 = self.type3
				print 'An extra angle parameter involving atom types %s %s %s was found, but the corresponding atoms seem to be missing. Skipping.' %(self.type1, self.type2, self.type3)
				return
		else:
			type1 = self.type1
			type2 = self.type2	
			type3 = self.type3	
		if self.extraparam1:
			#if self.extramemo:
			#	self.memo = self.memo + '; UB:' + self.extramemo
			filehandle.write( '%-7s%-7s%-7s  %s  %s  %s  %s  %s \n' %(type1, type2, type3, self.param2, self.param1, self.extraparam2, self.extraparam1, self.memo) )
		else:
			filehandle.write( '%-7s%-7s%-7s  %s  %s  %s \n' %(type1, type2, type3, self.param2, self.param1, self.memo) )		

class DihedralParams:
	""" Dihedrals are listed in "dihedral_trig". """
	def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, dihedparam1, dihedparam2, dihedparam3, dihedparam4, dihedparam5, dihedparam6, dihedparam7, dihedparam8, memo=None):
		self.type1 = atomtype1
		self.type2 = atomtype2
		self.type3 = atomtype3
		self.type4 = atomtype4
		self.phi0 = dihedparam1
		self.fcarray = [dihedparam2,dihedparam3,dihedparam4,dihedparam5,dihedparam6,dihedparam7,dihedparam8]
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '			
	def convertparams(self):
		""" Tries to figure out which linear combination of fc1 to fc6 yields the correct fc0 term, 
		using coefficients -1 or 1. This also determines the phase (0 or 180deg). """
		tol = 0.001
		self.a = []
		self.b = []
		self.c = []
		phi0 = dec(self.phi0)
		fc0 = dec(self.fcarray[0])
		nonzerofcs = [dec(fc) for fc in self.fcarray[1:] if dec(fc) != 0]
		nonzerinds = [i for i,fc in enumerate(self.fcarray[1:]) if dec(fc) != 0]
		nonzeronr = len(nonzerofcs)
		combos = list(product(*[[-1,1]]*nonzeronr))
		assert len(combos) == 2**nonzeronr
		combofound = False
		for c in combos:
			sum = dec(0.0)
			for sign, fc in zip(c, nonzerofcs):
				sum+=sign*fc
			if abs(sum - fc0) < tol:
				winningcombo = c
				combofound = True
				break
		if not combofound:
			print 'ERROR: No linear combo found when converting fc terms. Try setting a lower tolerance.'
		for i, sign, fc in zip(nonzerinds, winningcombo, nonzerofcs):
			self.b.append(i+1)
			if phi0 == 0.0 or phi0 == 180.00:
				self.a.append(sign*fc)
				if sign == -1:
					self.c.append(str(180.00))
				elif sign == 1:	
					self.c.append(str(0.00))
			else:
				self.a.append(fc)
				self.c.append(phi0)	
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				type1 = typedict[self.type1]
				type2 = typedict[self.type2]
				type3 = typedict[self.type3]
				type4 = typedict[self.type4]
			except KeyError:
				#type1 = self.type1
				#type2 = self.type2	
				#type3 = self.type3	
				#type4 = self.type4
				print 'An extra dihedral parameter involving atom types %s %s %s %s was found, but the corresponding atoms seem to be missing. Skipping.' %(self.type1, self.type2, self.type3, self.type4)
				return
		else:
			type1 = self.type1
			type2 = self.type2	
			type3 = self.type3	
			type4 = self.type4	
		if len(self.a)>0:
			for i in range(len(self.a)):
				filehandle.write( '%-7s%-7s%-7s%-7s  %s  %5d  %s  %s \n' %(type1, type2, type3, type4, self.a[i], self.b[i], self.c[i], self.memo) )	
		else:
			filehandle.write( '%-7s%-7s%-7s%-7s  %s  %5d  %s  %s \n' %(type1, type2, type3, type4, 0.0, 0, 0.0, self.memo) )		
	
class ImproperParams:
	""" These are listed in "improper_harm". """
	def __init__(self, atomtype1, atomtype2, atomtype3, atomtype4, impparam1, impparam2, memo=None):
		self.type1 = atomtype1
		self.type2 = atomtype2
		self.type3 = atomtype3
		self.type4 = atomtype4
		self.param1 = impparam1
		self.param2 = impparam2
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '			
	def convertparams(self):
		self.a = self.param2
		self.b = '0'
		self.c = '0.00' #self.param1
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				type1 = typedict[self.type1]
				type2 = typedict[self.type2]
				type3 = typedict[self.type3]
				type4 = typedict[self.type4]
			except KeyError:
				#type1 = self.type1
				#type2 = self.type2	
				#type3 = self.type3	
				#type4 = self.type4
				print 'An extra improper parameter involving atom types %s %s %s %s was found, but the corresponding atoms seem to be missing. Skipping.' %(self.type1, self.type2, self.type3, self.type4)
				return		
		else:
			type1 = self.type1
			type2 = self.type2	
			type3 = self.type3	
			type4 = self.type4
		filehandle.write( '%-7s%-7s%-7s%-7s %s  %s %s %s \n' %(type1, type2, type3, type4, self.a, self.b, self.c, self.memo) )

class NonbondedParams:
	""" These terms are specified in the files "vdw1" and "vdw1_14". """
	def __init__(self, atomtype, nbparam1, nbparam2, memo=None):
		self.type = atomtype
		self.param1 = nbparam1
		self.param2 = nbparam2
		if memo:
			self.memo = '! ' + memo
		else:
			self.memo = ' '			
		self.extraparam1 = None # from vdw1_14 file
		self.extraparam2 = None # from vdw1_14 file		
		self.extramemo = None
	def convertparams(self):
		""" Scales nonbonded terms by 2^(-5/6) """
		factor = dec(1.78179743628)# = 2^(5/6) = (2.0 / 2.0**(1./6.))
		self.diafcopa1 = (-1 * (dec(self.param1).as_tuple().exponent)) - 16 # chosen such that we get the same nr of sigfigs after ff_charm_to_viparr
		self.diafcopa2 = -1 * (dec(self.param2).as_tuple().exponent)
		self.a = dec(self.param1) / factor
		self.b = - dec(self.param2)
		if self.extraparam1:
			self.c = dec(self.extraparam1) / factor
			self.d = - dec(self.extraparam2)
	def writeto(self, filehandle, typedict=None):
		if typedict:
			try:
				atype = typedict[self.type]
			except KeyError:
				#atype = self.type	
				print 'An extra nonbonded parameter involving atom type %s was found, but the corresponding atom seems to be missing. Skipping.' %(self.type)
				return
		else:
			atype = self.type
		if self.extraparam1:
			filehandle.write( '%-7s 0.000000 %s  %s 0.000000 %s  %s  %s \n' %(atype, round(self.b, self.diafcopa2), round(self.a, self.diafcopa1), round(self.d, self.diafcopa2), round(self.c, self.diafcopa1), self.memo) )						
		else:	
			filehandle.write( '%-7s 0.000000 %s  %s  %s \n' %(atype, round(self.b, self.diafcopa2), round(self.a, self.diafcopa1), self.memo) )


def readviparr3files(viparrdir):
	""" Reads in all viparr3 files and populates lists with instances of corresponding objects. """
	# Read "mass"-file:
	try:
		mfname = join(viparrdir, 'mass')
		massfile = open(mfname, 'r')
		gotatomtypes = False
		gotparams = False
		j = 0
		for line in massfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype = splittedline[i+1]
					gotatomtypes = True
				elif w == 'amu':
					amu = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')
					except IndexError:
						memo = None	
			if gotatomtypes and gotparams:
				j += 1			
				m = Mass(atomtype,j,amu,memo) 
				masses.append(m)
				gotatomtypes = False
				gotparams = False
		massfile.close()		
	except IOError:
		print 'WARNING: No Viparr "mass" file found in directory %s' % viparrdir
	# Read "templates" file:
	tfname = join(viparrdir, 'templates')
	impropers = False
	bonds = False
	atoms = False	
	try:
		templatesfile = open(tfname, 'r')
		for line in templatesfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			if len(splittedline) == 0:
				continue
			elif 'cmap' in splittedline:
				print 'WARNING: CMAP terms found. This is not supported, and these terms will be skipped'
				#raise IOError	
				continue
			elif len(splittedline) == 1:
				if 'atoms' in splittedline[0]:
					atoms = True
					continue
				elif 'bonds' in splittedline[0]:
					bonds = True
					atoms = False
					continue
				elif 'impropers' in splittedline[0]:
					impropers = True
					bonds = False
					atoms = False
					continue			
				else:
					ligandname = splittedline[0]
					continue
			else: # Something is a bit weird; no line breaks between words and params... !!!
				if 'atoms' in splittedline[1]:
					print 'WARNING: "templates"-file is probably missing line breaks between words and params. Doublecheck your charges, bond and improper lists.'
					print 'WARNING: making a guess how to properly parse atoms in "templates"-file.'
					ligandname = splittedline[0]
					splittedline = list(splittedline[2:])
					atoms = True
				if 'bonds' in splittedline[0]:
					print 'WARNING: making a guess how to properly parse bonds in "templates"-file.'
					splittedline = list(splittedline[1:])
					atoms = False
					bonds = True	
				if 'impropers' in splittedline[0]:
					print 'WARNING: making a guess how to properly parse impropers in "templates"-file.'
					splittedline = list(splittedline[1:])
					atoms = False
					bonds = False
					impropers = True									
				elif not atoms and not bonds and not impropers:
					print "ERROR: Couldn't parse 'templates'-file."
					raise IOError
			if atoms:
				try:
					memo = ' '.join(splittedline[(4):])
				except IndexError:
					memo = None
				a = Atom(splittedline[0], splittedline[3], splittedline[2], splittedline[1], memo)
				a.convertelement()
				atomlist.append(a)
			elif bonds:
				bondlist.append( Bond(splittedline[0],splittedline[1]) )
			elif impropers:
				imps.append( Imp(splittedline[0],splittedline[1],splittedline[2],splittedline[3]) )
		templatesfile.close()		
	except IOError:
		raise IOError('ERROR: No Viparr "templates" file found in directory %s, or error while parsing that file. Aborting.' % viparrdir)
	# Read "stretch_harm"-file:
	try:
		bfname = join(viparrdir, 'stretch_harm')
		bondsfile = open(bfname, 'r')
		gotatomtypes = False
		gotparams = False
		for line in bondsfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype1 = splittedline[i+1]
					atomtype2 = splittedline[i+2]
					gotatomtypes = True
				elif w == 'r0':
					r0 = splittedline[i+1]
				elif w == 'fc':
					fc = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')
					except IndexError:
						memo = None		
			if gotatomtypes and gotparams:	
				b = BondParams(atomtype1,atomtype2,r0,fc,memo) 
				bonparams.append(b)
				gotatomtypes = False
				gotparams = False
		bondsfile.close()		
	except IOError:
		print 'ERROR: No Viparr "stretch_harm" file found.'
		raise IOError('Conversion needs "stretch_harm"-file. Make sure that your folder is in Viparr3 and not Viparr1 format.')
	# Read "angle_harm"-file:
	try:
		afname = join(viparrdir, 'angle_harm')
		anglesfile = open(afname, 'r')
		gotatomtypes = False
		gotparams = False
		for line in anglesfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':
					atomtype1 = splittedline[i+1]
					atomtype2 = splittedline[i+2]
					atomtype3 = splittedline[i+3]
					gotatomtypes = True
				elif w == 'theta0':
					theta0 = splittedline[i+1]
				elif w == 'fc':
					fc = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')
					except IndexError:
						memo = None
					break							
			if gotatomtypes and gotparams:	
				a = AngleParams(atomtype1, atomtype2, atomtype3, theta0, fc, memo)
				anglparams.append(a)
				gotatomtypes = False
				gotparams = False
		anglesfile.close()		
	except IOError:
		print 'ERROR: No Viparr "angle_harm" file found.'
		raise IOError('Conversion needs "angle_harm"-file.')
	# Read "ureybradley_harm"-file:
	try:
		ubfname = join(viparrdir, 'ureybradley_harm')
		ubfile = open(ubfname, 'r')
		gotatomtypes = False
		gotparams = False
		for line in ubfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'fc':
					fc = splittedline[i+1]
				elif w == 'r0':	
					r0 = splittedline[i+1]
					gotparams =True
				if w == 'type':
					atomtype1 = splittedline[i+1]
					atomtype2 = splittedline[i+2]
					atomtype3 = splittedline[i+3]
					gotatomtypes = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')						
					except IndexError:
						memo = None						
			if gotatomtypes and gotparams:	
				for a in anglparams:
					if a.type1 == atomtype1 and a.type2 == atomtype2 and a.type3 == atomtype3:
						a.extraparam1 = r0
						a.extraparam2 = fc
						a.extramemo = memo
						gotatomtypes = False
						gotparams = False
						break
		ubfile.close()				
	except IOError:
		print 'WARNING: No Viparr "ureybradley_harm" file found in directory %s' % viparrdir
	# Read "dihedral_trig"-file:
	try:
		dhfname = join(viparrdir, 'dihedral_trig')
		dihedfile = open(dhfname, 'r')
		gotatomtypes = False
		gotparams = False		
		for line in dihedfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype1 = splittedline[i+1]
					atomtype2 = splittedline[i+2]
					atomtype3 = splittedline[i+3]
					atomtype4 = splittedline[i+4]
					gotatomtypes = True
				elif w == 'phi0':
					phi0 = splittedline[i+1]
				elif w == 'fc0':
					fc0 = splittedline[i+1]
				elif w == 'fc1':
					fc1 = splittedline[i+1]
				elif w == 'fc2':
					fc2 = splittedline[i+1]
				elif w == 'fc3':
					fc3 = splittedline[i+1]
				elif w == 'fc4':
					fc4 = splittedline[i+1]
				elif w == 'fc5':
					fc5 = splittedline[i+1]
				elif w == 'fc6':
					fc6 = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')
					except IndexError:
						memo = None						
			if gotatomtypes and gotparams:
				d = DihedralParams(atomtype1, atomtype2, atomtype3, atomtype4, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6, memo)
				diheparams.append(d)
				gotatomtypes = False
				gotparams = False
		dihedfile.close()
	except IOError:
		print 'ERROR: No Viparr "dihedral_trig" file found.'
		raise IOError('Conversion needs "dihedral_trig"-file.')
	# Read "improper_harm"-file:
	try:
		impfname = join(viparrdir, 'improper_harm')
		impfile = open(impfname, 'r')
		gotatomtypes = False
		gotparams = False
		for line in impfile:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype1 = splittedline[i+1]
					atomtype2 = splittedline[i+2]
					atomtype3 = splittedline[i+3]
					atomtype4 = splittedline[i+4]
					gotatomtypes = True
				elif w == 'phi0':
					theta0 = splittedline[i+1]
				elif w == 'fc':
					fc = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')					
					except IndexError:
						memo = None						
			if gotatomtypes and gotparams:	
				imp = ImproperParams(atomtype1, atomtype2, atomtype3, atomtype4, phi0, fc, memo)
				impparams.append(imp)
				gotatomtypes = False
				gotparams = False
		impfile.close()
	except IOError:
		print 'WARNING: No Viparr "improper_harm" file found in directory %s' % viparrdir
	# Read "vdw1"-file:
	try:
		vdw1fname = join(viparrdir, 'vdw1')
		vdw1file = open(vdw1fname, 'r')
		gotatomtypes = False
		gotparams = False
		for line in vdw1file:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype1 = splittedline[i+1]
					gotatomtypes = True
				elif w == 'sigma':
					sigma = splittedline[i+1]
				elif w == 'epsilon':
					epsilon = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')					
					except IndexError:
						memo = None						
			if gotatomtypes and gotparams:	
				nb = NonbondedParams(atomtype1, sigma, epsilon, memo)
				nonbondparams.append(nb)
				gotatomtypes = False
				gotparams = False
		vdw1file.close()		
	except IOError:
		print 'WARNING: No Viparr "vdw1" file found in directory %s' % viparrdir
	# Read "vdw1_14"-file:
	try:
		vdw1_14fname = join(viparrdir, 'vdw1_14')
		vdw1_14file = open(vdw1_14fname, 'r')
		gotatomtypes = False
		gotparams = False		
		for line in vdw1_14file:
			strippedline = ''.join(c for c in line if c not in '"[]{},:')
			splittedline = strippedline.split()
			for i,w in enumerate(splittedline):
				if w == 'type':	
					atomtype1 = splittedline[i+1]
					gotatomtypes = True
				elif w == 'sigma':
					sigma = splittedline[i+1]
				elif w == 'epsilon':
					epsilon = splittedline[i+1]
					gotparams = True
				elif w == 'memo':
					try:
						ind1 = line.index('"memo": "')
						memo = ((line[ind1+9:].rstrip('\n')).rstrip(',')).rstrip('"}')					
					except IndexError:
						memo = None						
			if gotatomtypes and gotparams:	
				for nb in nonbondparams:
					if nb.type == atomtype1:
						nb.extraparam1 = sigma
						nb.extraparam2 = epsilon
						nb.extramemo = memo
						gotatomtypes = False
						gotparams = False
						break				
		vdw1_14file.close()			
	except IOError:
		print 'WARNING No Viparr "vdw1_14" file found in directory %s' % viparrdir
	# We need the name of the molecule later:	
	return ligandname


def readparamchem(viparrdir, viparratoms, viparrbonds):
	""" Maps new atomnames from paramchem file to old ones found in viparr file.
	Uses the VF2 algorithm from the networkx package. """

	class MoleculeMatcher(isomorphism.GraphMatcher):
		def semantic_feasibility(self, G1_node, G2_node):
			""" Overrides the default method that always returns True. Nodes have to be atom instances."""
			return G1_node.element == G2_node.element

	# Create Graph with new atomnames and types from paramchem:
	found = False
	for file in os.listdir(viparrdir):
		if file.endswith(".str"):
			pcfname = join(viparrdir, file)
			found = True
			print 'Paramchem stream-file "%s" found, thus converting atomtypes and names.' % pcfname
			break
	if not found:
		print 'WARNING: No paramchem stream-file found in directory %s, thus keeping original atomtypes and names.' % viparrdir
		return None			
	pcfile = open(pcfname, 'r')
	charmmatoms = []
	elements = {
		'HG': 'H', 
		'CG': 'C',
		'NG': 'N',
		'OG': 'O',
		'FG': 'F',
		'PG': 'P',
		'AL':'Al', 
		'SG': 'S', 
		'CL': 'Cl',
		'BR': 'Br',
		'IC': 'I', 
	}
	charmmGraph = nx.Graph()
	for l in pcfile:
		splittedline = l.split()
		if len(splittedline) == 0:
			continue
		elif splittedline[0] == 'ATOM':
			try:
				element = elements[ (splittedline[2][0:2]) ] 
			except KeyError: 	
				print 'ERROR: Unknown element for atom type %s' % splittedline[2]
				raise
			charmmatoms.append( Atom(splittedline[1], splittedline[2], splittedline[3], element) )
		elif splittedline[0] == 'BOND':
			for a in charmmatoms:
				if a.name == splittedline[1]:
					atom1 = a
				elif a.name == splittedline[2]:
					atom2 = a
			charmmGraph.add_edge(atom1,atom2)		
	pcfile.close()		
	# Create Graph with old atomnames and atomtypes from viparr:
	viparrGraph = nx.Graph()
	for b in viparrbonds:
		for a in viparratoms:
			if a.name == b.atom1:
				atom1 = a
			elif a.name == b.atom2:
				atom2 = a
		viparrGraph.add_edge(atom1,atom2)
	# Create mapping:
	Iso = MoleculeMatcher(viparrGraph, charmmGraph)
	isomorphy = Iso.is_isomorphic()
	if not isomorphy:
		print "ERROR: the molecule specified in the .str-file %s is not isomorphic with the one specified by viparr." % pcfname
		raise AssertionError
	return Iso.mapping


def createdicts(atomlist, mapping):
	""" Creates dictionaries to look up new paramchem/charmm atom names and types.
	'mapping' is the mapping dictionary from the networkx isomorphism algorithm. """
	namedict, typedict = {}, {}
	for a in atomlist:
		typedict[a.type] = mapping[a].type
		namedict[a.name] = mapping[a].name
	return namedict, typedict	


def writestream(filehandle, namedict, typedict):
	""" Writes .str-file by combining .par and .rtf parts. """
	filehandle.write( '* Parameters automatically converted from Viparr directory %s \n' %(dirold) )
	filehandle.write( '*\nread rtf card append\n\n' )
	writertf(filehandle, namedict, typedict)
	filehandle.write( '\nread param card flex append\n' )
	writepar(filehandle, namedict, typedict)
	filehandle.write( '\nEND\nRETURN\n' )


def writertf(filehandle, namedict, typedict):
	""" Writes .rtf-file; loops over lists of mass parameters, atom charges, bonds and impropers. """
	filehandle.write( '* Parameters automatically converted from Viparr directory %s \n' %(dirold) )
	if namedict:
		filehandle.write( '* Atomnames and -types from streamfile in %s \n*\n\n' %(dirold) )
	filehandle.write( '36 1\n\n' ) # Not sure what to put here (CHARMM version number)
	if len(masses) > 0:
		for m in masses:
			m.writeto(filehandle, typedict)
		filehandle.write( '\n' )
	filehandle.write( 'RESI %s %4.3f \n' %(ligandname, float(totalcharge)) )
	filehandle.write( 'GROUP\n' )	
	for ac in atomlist:
		ac.writeto(filehandle, namedict, typedict)
	for bc in bondlist:
		bc.writeto(filehandle, namedict)	
	for ic in imps:
		ic.writeto(filehandle, namedict)
	filehandle.write( 'END \n' )


def writepar(filehandle, namedict, typedict):
	""" Writes .par-file; loops over lists of bond parameters, angle parameters, dihedral paramters, improper parameters and nonbonded parameters. """
	filehandle.write( '* Parameters automatically converted from Viparr directory %s \n*\n\n' %(dirold) )
	filehandle.write( 'BONDS\n' ) 
	for b in bonparams:
		b.writeto(filehandle, typedict)
	filehandle.write( '\nANGLES\n' ) 	
	for a in anglparams:
		a.writeto(filehandle, typedict)
	filehandle.write( '\nDIHEDRALS\n' ) 	
	for d in diheparams:
		d.convertparams()
		d.writeto(filehandle, typedict)
	filehandle.write( '\n' ) 
	if len(impparams) > 0:
		filehandle.write( 'IMPROPER\n' ) 	
		for i in impparams:
			i.convertparams()
			i.writeto(filehandle, typedict)
		filehandle.write( '\n' )
	if len(nonbondparams) > 0:	
		filehandle.write( 'NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch - \n' ) # Unsure
		filehandle.write( 'cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5 \n' )
		for nb in nonbondparams:
			nb.convertparams()
			nb.writeto(filehandle, typedict)	


if __name__ == '__main__':
	# Handle arguments and paths:
	opt = optparse.OptionParser()
	opts, args = opt.parse_args()
	if len(args) != 2:
		opt.error('Please provide input (viparr3) and output directories (charmm) as arguments, seperated by a space.\n\
Optionally, place a stream file obtained from paramchem in the viparr3 directory to change atom names and types.')
	dirold=args[0]
	dirnew=args[1]
	assert dirold != dirnew
	if not os.path.exists(dirold):
		raise AssertionError('Input directory does not exist.')
	if os.path.exists(dirnew):
		raise AssertionError('Output directory already exists.')
	else:	
		os.makedirs(dirnew)
	print 'Input directory (Viparr3): %s' % dirold
	print 'Output directory (CHARMM): %s' % dirnew
	# Lists that will be populated with instances of the objects above
	masses = []
	atomlist = []
	bondlist = []
	imps = []
	bonparams = []
	anglparams = []
	diheparams = []
	impparams = []
	nonbondparams = []
	# Read in Viparr3 files
	try:
		ligandname = readviparr3files(dirold)
	except:
		os.rmdir(dirnew)
		raise		
	# Calculate total charge
	totalcharge = 0.0
	for ac in atomlist:
		totalcharge += float(ac.charge)
	# Convert atom types found in Viparr to CHARMM-compatible atom 
	# types by reading in .str frile generated by paramchem, if present
	mapping = readparamchem(dirold, atomlist, bondlist)
	if mapping:
		namedict, typedict = createdicts(atomlist, mapping)
	else:
		namedict, typedict = None, None
	# Write to .str-file
	outputfilesprefix = join(dirnew, ligandname)
	streamfile = outputfilesprefix+'.str'
	stream = open(streamfile, 'w')
	writestream(stream, namedict, typedict)
	# Write to .par-file
	parfile = outputfilesprefix+'.par'
	par = open(parfile, 'w')
	writepar(par, namedict, typedict)
	# Write to .rtf-file
	rtffile = outputfilesprefix+'.rtf'
	rtf = open(rtffile, 'w')
	writertf(rtf, namedict, typedict)
	print 'Conversion successful.'
	