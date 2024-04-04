import os
HOME = os.getenv('HOME')
APPS = os.path.join(HOME,'Library')
CWPROOT = os.path.join(APPS,'CWP/44r16')
RSFROOT = os.path.join(APPS,'madagascar/rsfroot')
TRIP = os.path.join(HOME,'work/trip')
MPIROOT = '/usr/local'
os.environ['CWPROOT'] = CWPROOT
os.environ['RSFROOT'] = RSFROOT
os.environ['TRIP'] = TRIP
os.environ['MPIROOT'] = MPIROOT

