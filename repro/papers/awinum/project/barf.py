import os
HOME = os.getenv('HOME')
print(HOME)
APPS = os.path.join(HOME,'Library')
print(APPS)
CWPROOT = os.path.join(APPS,'CWP/44r16')
print(CWPROOT)
RSFROOT = os.path.join(APPS,'madagascar/rsfroot')
print(RSFROOT)
TRIP = os.path.join(HOME,'work/trip')
print(TRIP)
PYRVL = os.path.join(TRIP,'pyrvl')
print(PYRVL)
PYTHONPATH = os.getenv('PYTHONPATH')
print(PYTHONPATH)

import data
data.hello()



