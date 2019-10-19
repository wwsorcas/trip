import os

def getnum(filename, key):
    val=0.0
    pathfile=os.path.join((os.getcwd(),filename)
    if os.path.exists(pathfile):
        f = open(pathfile,'r')
        for line in f:
            alist = (f.read().strip('\n')).split('=')
            if alist[0] == key:
                val = float(alist[1])
            f.close()
    return val
