#!/usr/bin/env python
### Written by Andrey I. Frolov, Institute of Solution Chemistry RAS, Russia. March 2012.
# Andrey I. Frolov, Jan 2014, ISC RAS, Ivanovo, Russia
#
#   Copyright 2014 Andrey I. Frolov
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

# This is the main program. It parses the input options and runs the corresponding command file. 

# Last modified 28 oct 2014

# v 1.3: 28 oct 2014
# - Bug fixed: Different LJ parameters can correspond to the same atom type in ffld_server output (see "types" collumn). Now the program define the atomtype name as a combination of "type" and "vdw" collumns of ffld_server output. Therefore, atomtypes are unique. Only conversion from ffld_server format was affected. Many solutes with many different hydrogen atomtypes were affected.  
#   Function ReadFFLDServerOutput fixed.
# - Changes: now longer atomtype names are defined for the conversion from ffld_server output.

# v 1.2: 26 jul 2014
# - Bug fixed: wrong automatic assignment of pairs. Only special solutes were affected: 
#   e.g. with five-member rings. Also most of compounds with nitrogen atoms. Function GenPairsFromBondConnectivity fixed.

# v 1.1: 
# - added convertion Fourier type dihedral potentials to Ryckaert-Bellemans (RB)
# - set the convertion to RB type by default for the output of Gromacs topology
# - started tracking versions

import sys, math, commands, os, shutil, string
from molecule_class import *

from func_namd_io import *
from func_gromacs_io import *
from func_ffldserver_io import *
from func_amber_io import *
from func_charmm_io import *
from func_rismmol_io import *

###------------------------------------------------------
### Options parser
from optparse import OptionParser

usage = "ffconv.py performs conversion of molecular topology and force field files between different formats. \n Example usage: $FFCONVPATH/ffconv.py --fromprog FFLDSERVER --toprog GROMACS --ffname oplsaa --moltopfn paracetamol.ffld --ffconvpath $FFCONVPATH "

parser = OptionParser(usage=usage)
parser.add_option("--ffconvpath", 
                  default='',
                  dest="ffconvpath", type="string", 
                  help="path to ffconv")
parser.add_option("--fromprog", 
                  default='',
                  dest="fromprog", type="string", 
                  help="program to convert from")
parser.add_option("--toprog", 
                  default='',
                  dest="toprog", type="string", 
                  help="program to convert to")
parser.add_option("--ffname", 
                  dest="ffname", type="string", 
                  default='',
                  help="name of the force field")
parser.add_option("--moltopfn", 
                  dest="moltopfn", type="string", 
                  default='',
                  help="file with molecular topology in FROMPROG format")
parser.add_option("--molcoorfn", 
                  dest="molcoorfn", type="string", 
                  default='',
                  help="file with coordinates")
parser.add_option("--FFfn", 
                  dest="FFfn", type="string", 
                  default='',
                  help="file with force field in in FROMPROG format")
parser.add_option("--execfile", 
                  dest="customfn", type="string", 
                  default='',
                  help="custom file with ffconv commands")
parser.add_option("--molname", 
                  dest="molname", type="string", 
                  default='',
                  help="name of the molecule. (to generate output file names)")


(opts, args) = parser.parse_args()


### Set options
ffconvpath = opts.ffconvpath
fromprog = opts.fromprog
toprog = opts.toprog
ffname = opts.ffname

moltopfn = opts.moltopfn
molcoorfn = opts.molcoorfn
FFfn = opts.FFfn
molname = opts.molname

# Get molname from the filename
if len(molname)==0:
    tmp = moltopfn.split('/'); tmp=tmp[-1]; tmp=tmp.split('.'); molname=tmp[0]
    

# VARS
customfn = opts.customfn
if len(customfn)==0: 
    lCustomCmdFile = False
else:
    lCustomCmdFile = True
    

execfn = ffconvpath + '/ffconv_'+fromprog+'_to_'+toprog+'.runcmd'


### Run commands in the cmd file
if not lCustomCmdFile:
    # Run predefined conversion cmds
    execfile(execfn) 

else:
    # Run cmd from custom file
    execfile(customfn) 
    
sys.exit()

