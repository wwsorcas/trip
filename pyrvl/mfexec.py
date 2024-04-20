#!/usr/bin/env python3

import os
import sys
import scenv

if __name__ == '__main__':   
    try:

#        print('check enough args')
        if len(sys.argv) < 2:
            raise Exception('must be called with at least one argument ' + \
                                '= name of function')
                                
#        print('sys.argv[1]')
#        print(sys.argv[1])

        arga = sys.argv[1].split('.')
#        print(arga)
        if len(arga) < 2:
            raise Exception('first arg must have form mod.fcn')

#        print('import ' + arga[0])
        exec('import ' + arga[0])

        argl = ', '.join(sys.argv[2:])

        cmd = sys.argv[1] + '(' + argl + ')'
            
#        print('function call:')
#        print(cmd)

        exec(cmd)
        
    except Exception as ex:
        print(ex)
        sys.exit(1)
        
        

        
