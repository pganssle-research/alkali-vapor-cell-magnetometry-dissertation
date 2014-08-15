from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

def contour3d(x, y, z, D):
	# Plot a stack of contours.
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	for ii in range(0,len(z)):
		ax.contour(x,y,D[:,:,ii],offset=z[ii])
		
	ax.set_zlim(min(z),max(z))
	plt.draw()