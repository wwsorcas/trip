import matplotlib.pyplot as plt
import numpy as np

fileName = 'coronalSlice750.npy'
inSlice = np.load(fileName)
plt.figure()
plt.imshow(np.squeeze(inSlice[0]), cmap='jet', origin = 'lower')
plt.title(fileName + ':  rho')
plt.colorbar()
plt.figure()
plt.imshow(np.squeeze(inSlice[1]), cmap='jet', origin = 'lower')
plt.title(fileName +':  velP')
plt.colorbar()
plt.figure()
plt.imshow(np.squeeze(inSlice[2]), cmap='jet', origin = 'lower')
plt.title(fileName + ':  tauP')
plt.colorbar()

plt.show()

