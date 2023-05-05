import numpy as np
import array

# setup
model=np.zeros((4))
model[1]=1.0
model[2]=2.0
model[3]=3.1415927

# use ndarray.tofile
model.tofile('doubles')

# use array.tofile
modelarray=array.array('f')
modellist=model.tolist()
modelarray.fromlist(modellist)
modelfile = open('singles','wb')
modelarray.tofile(modelfile)
modelfile.close()
