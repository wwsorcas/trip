import scenv
import data
import linalg

import sys

NX  = int(sys.argv[1])
NZ  = int(sys.argv[2])
DX  = float(sys.argv[3])
FAC = float(sys.argv[4])
RAD = float(sys.argv[5])

m = sys.argv[6]

data.model(bulkfile=m, bulk=4.0, nx=NX, nz=NZ,
   dx=DX, dz=DX, lensfac=FAC, lensradd=RAD)

