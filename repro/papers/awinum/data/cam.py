import numpy as np
import math
import data
import linalg
import sys

def cam(rad=800, inner=2.0, outer=4.0, n1=201, n2=401, d1=20, d2=20):
    try:

        # create rsf files
        data.rsffile('cambulk.rsf', 'Bulk_modulus', 'GPa', n1, n2, d1, d2, val=1.0)
        data.rsffile('cambuoy.rsf', 'Buoyancy', 'cc/g', n1, n2, d1, d2, val=1.0)

        # create array
        z = np.ones(n1*n2)
        cam = z.reshape(n2,n1)

        # loop through array scaling as appropriate
        for i in range(n1):
            z = i*20.0 - 2000.0
            for j in range(n2):
                x = j*20 - 4000.0
                d = math.sqrt(x*x + z*z)
                if d<rad:
                    cam[j,i] *= inner
                else:
                    cam[j,i] *= outer

        # write bulkmod to rsf file
        linalg.ndarraytorsfdata(cam,'cambulk.rsf')

    except Exception as ex:
        print(ex)
        raise Exception('called from cam')

if __name__ == '__main__':
    args = ", ".join(sys.argv[1:])
    cmd = 'cam(' + args + ')'
    exec(cmd)

