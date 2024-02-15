#!/usr/bin/env python3

import data
import sys
import os
import data

fkeys = {    'data.bpfilt': ['file','nt','dt','s','f1','f2','f3','f4','sx','sz'],
             'data.bpfiltgather': ['file','nt','dt','s','f1','f2','f3','f4','ntr','sxstart','szstart','dsx','dsz'],
             'data.delta': ['file','nt','dt','sx','sz'],
             'data.deltagather': ['file','nt','dt','ntr','sxstart','szstart','dsx','dsz'],
             'data.rechdr': ['file','nt','dt','ntr','rx','rz','sx','sz','drx','delrt','nshot','dsx','fixed'],
             'data.rsffile': ['file', 'datatype', 'unit', 'n1', 'n2', 'd1', 'd2', 'val'],
             'data.model': ['bulkfile', 'bulk', 'nx', 'nz', 'dx', 'dz', 'lensfac', 'buoy', 'lensradd', 'lensradt'],
        }

if __name__ == '__main__':
    try:

        if len(sys.argv) < 2:
            raise Exception('must be called with at least one argument ' + \
                                '= name of function')
        fcn = sys.argv[1]
        if not fcn in fkeys:
            raise Exception('function ' + fcn + ' not accessible')
        cmd = fcn + '('

        arglist = fkeys[fcn]

        args = dict(arg.split('=') for arg in sys.argv[2:])
        for k in range(len(fkeys[fcn])):
            if not fkeys[fcn][k] in args:
                raise Exception('argument ' + fkeys[fcn][k] + ' for function ' + fcn + ' not provided')

            if len(fkeys[fcn][k]) == 1:
                cmd + fkeys[fcn][k][0]
        if fcn=='data.bpfilt':
            data.bpfilt(args['file'],
                       int(args['nt']),
                       float(args['dt']),
                       float(args['s']),
                       float(args['f1']),
                       float(args['f2']),
                       float(args['f3']),
                       float(args['f4']),
                       float(args['sx']),
                       float(args['sz']))
        if fcm=='data.bpfiltgather':
            
        
    except Exception as ex:
        print(ex)
