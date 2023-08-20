import linalg
import numpy as np

try:
    # 1D basic case
    x = np.zeros([21,31])
#    linalg.ndarraytorsf(x,'nprsf1.rsf')

    # add unit
#    linalg.ndarraytorsf(x,'nprsf2.rsf', unit='GPa')

    # add axis unit
#    linalg.ndarraytorsf(x,'nprsf3.rsf', unit='GPa', units=['m'])

    # 2D basic case
#    linalg.ndarraytorsf(x,'nprsf4.rsf', n=[21,31])

    # 2D decorated
    linalg.ndarraytorsf(x,'nprsf5.rsf', n=[21,31], d=[10, 10], o=[0.0, 1500], unit='GPa', units=['m', 'm'])
    
except Exception as ex:
    print(ex)


