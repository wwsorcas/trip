import numpy as np
import linalg
import data
import os

def runtests():
    try: 
        # setup
        RSFROOT = os.getenv('RSFROOT')
        if not os.path.exists(RSFROOT):
            raise Exception('root M8R directory RSFROOT = ' + RSFROOT + ' not found')
        cmdin = os.path.join(RSFROOT,'bin/sfin')
        cmdat = os.path.join(RSFROOT,'bin/sfattr')
        cmdrm = os.path.join(RSFROOT,'bin/sfrm')
        if data.rsffile('nprin.rsf', 'zebras', 'animal', 101, 201, 20, 10, val=2.0) != 0:
            raise Exception('rsffile failed')

        print('\ninput file:')
        ret=os.system(cmdin + ' < nprin.rsf')                  
        ret=os.system(cmdat + ' < nprin.rsf')                  
        
        # extract ndarray
        x = linalg.ndarrayfromrsfdata('nprin.rsf')

        print('\nextracted ndarray:')
        print('size = ' + str(x.size) + ' max = ' + str(x.max()) + ' min = ' + str(x.min()))
        if data.rsffile('nprout.rsf', 'roses', 'flower', 101, 201, 20, 10, val=3.0) != 0:
            raise Exception('rsffile failed')

        print('\noutput before overwrite')
        ret=os.system(cmdin + ' < nprout.rsf')                  
        ret=os.system(cmdat + ' < nprout.rsf')

        linalg.ndarraytorsfdata(x,'nprout.rsf')
        
        print('\noutput after overwrite')
        ret=os.system(cmdat + ' < nprout.rsf')

        # clean up
        os.system(cmdrm + ' nprin.rsf')
        os.system(cmdrm + ' nprout.rsf')
        
    except Exception as ex:
        print(ex)

    # marmousi read - write, change unit from m/s (wrong) to m/ms (right)
    try:
        cmddd=os.path.join(RSFROOT,'bin/sfdd')
        os.system(cmddd + ' form=native < marm_vel.HH > vel.rsf')
        linalg.simplot('vel.rsf',addcb=True)
        if data.rsffile('velcopy.rsf', 'Velocity', 'm/ms', 751, 2301, 4, 4, val=1.0) != 0:
            raise Exception('rsffile failed')
        # extract ndarray
        v = linalg.ndarrayfromrsfdata('vel.rsf')
        linalg.ndarraytorsfdata(v,'velcopy.rsf')
        linalg.simplot('velcopy.rsf',addcb=True)
        print('\nMarmousi vel after correction of unit m/s to m/ms:')
        os.system(cmdin + ' < velcopy.rsf')
        os.system(cmdat + ' < velcopy.rsf')

        os.system(cmdrm + ' vel.rsf')
        os.system(cmdrm + ' velcopy.rsf')
        
    except Exception as ex:
        print(ex)        

runtests()




