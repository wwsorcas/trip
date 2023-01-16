####### IWAVE/IWAVE++ numerical experiment SConstruct template ########
#
# Stores results (Flow targets) in working directory, from two 
# kinds of flows:
# 
# I. flows invoking non-IWAVE commands: these should be written 
# in normal Madagascar fashion, and placed in the 'local commands'
# section of this SConstruct.
# 

# II. IWAVE flows: defined in the flow loop at the bottom of this
# file, which the user should not alter. This utility assumes that
# IWAVE/IWAVE++ commands will read their input from parameter files,
# which should reside in this working directory. A parameter file
# is identified by the suffix ".par". 

# Note that other commands which read all input data from a single
# parameter file with suffix ".par" should work the same way - will
# refer to all such commands as "IWAVE commands".

# For parfile 'foo.par', 'foo' is a key in the jobs dictionary. an
# entry in jobs looks like

# 'foo' : { 'cmd' : <command name>, 
#           'src' : <dependency list>, 
#           'tgt' : <target dictionary>,
#           'exe' : (optional) <execution environment dictionary> }
# 
# That is, the root parfile name serves as the keyword, and the
# associated value is another dictionary. The four keywords of 
# the value dictionary and their associated values are
#
# 'cmd' : path to the IWAVE/IWAVE++ command 
# 'src' : results on which this job depends - each result is
#         a file, which is a target of a Flow defined in this 
#         SConstruct, and resides in this directory
# 'tgt' : another dictionary, in which each keyword is a 
#         keyword in the parfile specifying an output file, and 
#         the corresponding value is the filename under 
#         which the result should be archived in this 
#         working directory. Each such file is a Flow
#         target. So for example if foo[tgt] contains the pair 
#         'datafile : mydata.su', then the flow for target 'mydata.su'
#         moves the file associated to keyword 'datafile' in the 
#         job parfile to the file mydata.su from whereever it currently
#         resides (usually a scratch directory for storing various process
#         outputs from another flow) to "mydata.su" in this working 
#         directory.

# 'exe' : (optional) yet another dictionary, specifying execution
#         environment.  If given, the keyword 'platf', specifies
#         execution platform, which falls into one of three major
#         categories:

#         SERIAL - If either 'exe' is not a key, or if the value of
#         'platf' is '', execution takes place in the foreground, by
#         serial (non-mpi) command line invocation of the command.

#         COMMAND-LINE MPI: A second predefined value of 'platf' is
#         'mpi', for command-line mpi execution: for this parallel
#         mode, you must specify integer values for the keyword 'ppn'
#         (processes per node - by convention, this form of
#         parallelism takes place on a single "node").

#         Thus foo:exe defined as
#                     'exe'     : { 'platf' : 'mpi',
#                                   'ppn'   : '16',
#                                 }
#         with foo:cmd = <cmd> results in the foreground execution of 
#             mpirun -np 16 <cmd> par=foo.par

#         BATCH: Other 'platf' values should be descriptive names for
#         environments in which the user intends to run iwave jobs
#         under batch control. Additional parameters for batch
#         submissions are:

#             wall  = Wallclock time limits (hr:mn:se)
#             nodes = number of nodes
#             ppn   = processes/cores/threads per node

#         Each batch environment name must be a keyword in the 'penv'
#         dictionary appearing just before the jobs dictionary. 'penv'
#         defines several standard inputs needed by all batch
#         systems. Here is an example:

#             penv = { (...)
#                     'stampede' : { 'batch' : 'slurm',
#                                    'queue' : 'normal',
#                                    'acct'  : 'FDTD3D-Cont',
#        	                     'mail'  : 'symes@caam.rice.edu',
#                                  }
#                      (...)
#                    }
#         The value of penv:platf:batch must be a keyword in the
#         bsys dictionary, defined in the common code section below
#         bsys defines various implementation-indedpendent attributes
#         of batch systems - see the code below for the currently 
#         accommodated systems, and feel free to add another and 
#         share!
#
#         So foo:exe defined as follows:
#                     'exe'     : { 'platf' : 'stampede',
#                                   'wall'  : '00:59:00',
#                                   'nodes' : '4',
#                                   'ppn'   : '16',
#                                  }

#         would cause the a batch file foo.bat to be written, and
#         submitted for execution on 4 nodes, 16 cores per node, of
#         environment 'stampede' via invocation of 'sbatch'

#         Thus 'bsys:batch' specifies attributes common to all
#         implementations of the batch system 'batch', 'penv:platform'
#         specifies batch parameters common to all uses of 'platform'
#         (including a choice of batch system), and 'jobs:jobname:exe'
#         specifies parameters particular to the job 'jobname',
#         including the platform on which it is to be executed.
 
# IWAVE commands execute in scratch (sub)directories generated by the flows that 
# invoke them; these flows generate the work directory names from the parfiles. 
# For example, the parfile foo.par generates a commond which executes in 
# foo.work. The subdirectories appear in WORKPATH, specified in the 
# common definitions section immediately following

import os
from trip.hosts import getParHosts

THISPATH        = os.getcwd()

def getCommand(jobdict, locdict):

    envdict = getParHosts()
    
#    print jobdict
#    print envdict

    if 'exe' not in jobdict.keys():
        jobdict['exe'] = getXenv(envdict,locdict)

    OMP_NUM_THREADS = os.getenv('OMP_NUM_THREADS')
    if not OMP_NUM_THREADS:
        OMP_NUM_THREADS = 1
    
    # batch execution - input to sfbatch
    if isBatch(jobdict,envdict):

        workcmd = jobdict['pre'] + '; ' + \
            'export OMP_NUM_THREADS=' + str(OMP_NUM_THREADS) + '; ' + \
            os.path.join(os.getenv('RSFROOT'),'bin/sfbatch') + \
            ' exe=' + '"' + \
            envdict[jobdict['exe']['platf']]['launcher'] + ' ' + \
            jobdict['cmd'] + '"' +\
            ' job=' + jobdict['job'] + \
            ' wall=' + jobdict['exe']['wall'] + \
            ' nodes=' + str(jobdict['exe']['nodes']) + \
            ' ppn=' + str(jobdict['exe']['ppn']) + \
            ' np=' +  str(getThreads(locdict)) + \
            ' batch=' + envdict[jobdict['exe']['platf']]['batch'] + \
            ' queue=' + envdict[jobdict['exe']['platf']]['queue'] + \
            ' mail=' + envdict[jobdict['exe']['platf']]['mail'] + \
            ' acct=' + envdict[jobdict['exe']['platf']]['acct'] + \
            ' path=' + os.getcwd() + '/' + jobdict['job'] + '.work' 

    elif isMPI(jobdict):

        # check to see if MPIRUN is defined in environment
        MPIRUN = os.getenv('MPIRUN')
        if not MPIRUN:
        # else try standard choice            
            MPIROOT = os.getenv('MPIROOT')
            MPIRUN = os.path.join(MPIROOT,'bin/mpirun')
        # print 'MPIRUN=' + MPIRUN
        # this should be checked for existence
        workcmd = jobdict['pre'] + '; ' + \
            'export OMP_NUM_THREADS=' + str(OMP_NUM_THREADS) + '; ' + \
            MPIRUN + ' -np ' + str(jobdict['exe']['ppn']) + \
            ' ' + jobdict['cmd']
        # print workcmd
    else:      

        workcmd = jobdict['pre'] + '; ' + \
            'export OMP_NUM_THREADS=' + str(OMP_NUM_THREADS) + '; ' + \
            jobdict['cmd']

    return workcmd

def isBatch(jobdict,envdict):
    if ('exe' in jobdict.keys()):
        # check spec of platform
        if ('platf' in jobdict['exe'].keys() and \
            'wall' in jobdict['exe'].keys() and \
            'nodes' in jobdict['exe'].keys() and \
            'ppn' in jobdict['exe'].keys()):
            # batch case - spec'd in envdict 
            if (jobdict['exe']['platf'] in envdict.keys()) and ('batch' in envdict[jobdict['exe']['platf']]):
                return True
            else:
                return False
        else:
            return False
    else:
        return False

def hasMPIRUN(jobdict):
    MPIRUN = os.getenv('MPIRUN')
    if MPIRUN:
        return MPIRUN
    MPIROOT = os.getenv('MPIROOT')
    if MPIROOT and os.path.exists(os.path.join(MPIROOT,'bin/mpirun')):
        return os.path.join(MPIROOT,'bin/mpirun')
    else:
        print('command line MPI via mpirun not available')
        return ''

def isMPI(jobdict):
# test for presence of MPI info, for command line MPI
    if ('exe' in jobdict.keys()):
        # check spec of platform
        if ('platf' in jobdict['exe'].keys() and \
            'ppn' in jobdict['exe'].keys()):
            if (jobdict['exe']['platf'] == 'mpi'):
                if hasMPIRUN(jobdict):
                    return True
                return False
        else:
            return False
    else:
        return False
    
def getThreads(lenv):
    penv = getParHosts()
    exe  = getXenv(penv,lenv)
    nthread=1
    if ('ppn' in exe.keys()):
        nthread = nthread * int(exe['ppn'])
    if ('nodes' in exe.keys()):
        nthread = nthread * int(exe['nodes'])
    return nthread

def getPlatform(penv):
    HOSTNAME = str(os.getenv('HOSTNAME'))
#    print HOSTNAME
    hna = HOSTNAME.split(".")
#    hna = Split(HOSTNAME,".")
    for host in penv.keys():
        if (host in hna):
            return host
    return ''

def printPlatform(penv):
    print(getPlatform(penv))

def getXenv(penv,lenv):
    platf = getPlatform(penv);
    if platf in lenv.keys():
        if 'batch' in penv[platf].keys():
            if ('nodes' in lenv[platf].keys() and
                'ppn' in lenv[platf].keys() and
                'wall' in lenv[platf].keys()):
                exe = {'platf': platf, 'nodes': lenv[platf]['nodes'], 'ppn': lenv[platf]['ppn'], 'wall': lenv[platf]['wall'] }
                return exe
            else:
                exe = {}
                return exe
        else:
            if 'ppn' in lenv[platf].keys():
                if lenv[platf]['ppn'] > 1:
                    exe = {'platf': 'mpi', 'ppn': lenv[platf]['ppn'] }
                    return exe
                else:
                    exe = {}
                    return exe
    else:
        exe = {}
        return exe
                               
def tripExec(jobs, lenv):
    for i in range(len(jobs)):
        cmd = getCommand(jobs[i], lenv)
        if cmd == None:
            print('Error return from newbatch.tripExec - cannot set up jobs['+str(i)+']')
        else:
            Flow(jobs[i]['tgt'], jobs[i]['src'], cmd,
                 stdin=0, stdout=-1, workdir=jobs[i]['job']+'.work')

def printJobs(jobs, lenv):
    for i in range(len(jobs)):
        cmd = getCommand(jobs[i], lenv)
        if cmd == None:
            print('Note: jobs['+str(i)+' not defined]')
        else:
            print('cmd = ' + cmd + '\n')
            print('src = ' + ' '.join(jobs[i]['src']) + '\n')
            print('tgt = ' + ' '.join(jobs[i]['tgt']) + '\n')

def getnum(filename, key):
#    print('enter getnum')
#    print(filename)
#    print(key)
    val=0.0
    pathfile=os.path.join(os.getcwd(),filename)
    if os.path.exists(pathfile):
        f = open(pathfile,'r')
        for line in f:
#            print(line)
            alist = (line.strip('\n')).split('=')
#            print(alist)
            if alist[0] == key:
                val = float(alist[1])
        f.close()
#        print('val='+str(val))
#        print('exit getnum')
    return val            
