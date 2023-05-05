import linalg
import os
import m8r

def su_to_array(f):



    '''

    input .su file and output ndarray

    '''



    RSFROOT = os.getenv('RSFROOT')

    if not os.path.exists(RSFROOT):

        print('su_to_array')

        print('root M8R directory RSFROOT = ' + RSFROOT + ' not found')

        print('M8R not correcly installed')

        return False

    if linalg.sanity(f,'su'):

        ff = f[0:len(f)-3]+'.rsf'

        cmd = os.path.join(RSFROOT,'bin/sfsuread')

        ret=os.system(cmd + ' read=data endian=0 < ' + f + ' > ' + ff)

        if ret != 0:

            print('su_to_array')

            print('command:')

            print(cmd + ' read=data endian=0 < ' + f + ' > ' + ff)

            print('failed with return value ' + str(ret))

            return False

        asprat=1

    elif linalg.sanity(f,'rsf'):

            ff=f

    else:

        print('su_to_array: arg does not name file will either .su or .rsf suffix')

        return False



    # output as array





    inp=m8r.Input(ff)

    datafile=str(inp.get('in'))

    data = inp.read()

    if linalg.sanity(f,'su'):

        os.unlink(ff)

        os.unlink(datafile[2:len(datafile)-3])

    return data

x = su_to_array('/var/tmp/baru.su')
