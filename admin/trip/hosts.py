import os

def getParHosts():

    penv = {'stampede2' : { 'batch'    : 'slurm',
                                'queue'    : 'normal',
                                'acct'     : 'FDTD3D-Cont',
                                'mail'     : 'symes@caam.rice.edu',
                                'launcher' : 'export OMP_NUM_THREADS=1;' +
                                'export I_MPI_ROOT=' +
                                str(os.getenv('I_MPI_ROOT')) +
                                '; ibrun'
                         },
            'davinci' :  { 'batch'    : 'slurm',
                           'queue'    : 'trip',
                           'acct'     : 'trip',
                           'mail'     : 'symes@caam.rice.edu',
                           'launcher' : 'srun'
                         },
            'macbook' :  {},
            'obelix'  :  {},
            'getafix' :  {},
            }
    return penv
