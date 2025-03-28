from rsf.proj import *
from trip.proj import getCommand
from trip.hosts import getParHosts
import os

def exec(jobs, lenv):
    penv = getParHosts()
    for i in range(len(jobs)):
        cmd = getCommand(jobs[i], penv, lenv)
        if cmd == None:
            print 'Error return from trip.exec.tripExec - cannot set up jobs['+str(i)+']'
        else:
            Flow(jobs[i]['tgt'], jobs[i]['src'], cmd,
                 stdin=0, stdout=-1, workdir=jobs[i]['job']+'.work')

