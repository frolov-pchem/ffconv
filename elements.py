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
from molecule_class import *

Elements = {
 'H':
     {'Mass':1.0079, 'PeriodicTableNum':1 },
 'C':
     {'Mass':12.011, 'PeriodicTableNum':6 },
 'O':
     {'Mass':15.999, 'PeriodicTableNum':8 },
 'N':
     {'Mass':14.007, 'PeriodicTableNum':7 },
 'P':
     {'Mass':30.973762, 'PeriodicTableNum':15 },
 'S':
     {'Mass':32.07, 'PeriodicTableNum':16 },
 'F':
     {'Mass':18.9984032, 'PeriodicTableNum':9 },
 'CL':
     {'Mass':35.453, 'PeriodicTableNum':17 },
 'BR':
     {'Mass':79.904, 'PeriodicTableNum':35 },
 'I':
     {'Mass':126.90447, 'PeriodicTableNum':53 },
 'AL':
     {'Mass':26.9815386, 'PeriodicTableNum':13 }
}

def TrimElementName(s):
    El = s[:2].upper()
    if El in Elements.keys():
        return El
    else:
        if El[:1] in Elements.keys():
            return El[:1]
        else:
            sys.stderr.write("--- !Warning in "+inspect.stack()[0][3]+": could not get an element name from the string ["+str(s)+"]. Maybe extend the elements list? Known elements: ["+str(Elements)+"]. Returning none. \n")
            return 'none'
        
        

