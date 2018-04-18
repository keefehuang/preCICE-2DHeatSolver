import matplotlib.pyplot as plt
import numpy as np
import configuration_file as conf
from mpl_toolkits.mplot3d import Axes3D
import configuration_file as config


def plotVar(ax, sharedBound, maxT):
	plt.plot(sharedBound)


def plotTempLeft(ax, temp, maxT):
	N = conf.n_elem
	x = np.linspace(-1,0, num=N+2)
	y = np.linspace(-1,0, num=N+1)
	X,Y = np.meshgrid(x, y, indexing='xy')
	ax.plot_surface(X,Y,temp)
	ax.set_zlim([0,maxT+5])

def plotTempRight(ax, temp, maxT):
	N = conf.n_elem
	x = np.linspace(0,1, num=N+2)
	y = np.linspace(0,1, num=N)
	X,Y = np.meshgrid(x, y, indexing='xy')
	ax.plot_surface(X,Y,temp)
	ax.set_zlim([0,maxT+5])

def plotTempWhole(ax, temp, maxT):
	N = conf.n_elem
	x = np.linspace(0,1, num=N+2)
	y = np.linspace(0,1, num=2*N+3)
	X,Y = np.meshgrid(x, y, indexing='xy')
	ax.plot_surface(X,Y,temp)
	ax.set_zlim([0,maxT+5])

def doPlottingLeft(ax, sharedBound, temp, t, maxT):
	# plotVar(ax, sharedBound, maxT)
	plotTempLeft(ax, temp, maxT)
	plt.title(t)
	plt.pause(0.1)

def doPlottingRight(ax, sharedBound, temp, t, maxT):
	# plotVar(ax, sharedBound, maxT)
	plotTempRight(ax, temp, maxT)
	plt.title(t)
	plt.pause(0.1)

def doPlottingWhole(ax, temp, t, maxT):
	plotTempWhole(ax, temp, maxT)
	plt.title(t)
	plt.pause(0.1)

if __name__ == '__main__':
	N = config.n_elem
	init_temp = config.init_temp
	bound_temp = config.bound_temp
	t = 0
	maxT = 40

	fig = plt.figure(1)
	ax = fig.add_subplot(111, projection='3d')

	fig2 = plt.figure(2)
	ax2 = fig.add_subplot(111, projection='3d')

	botNodes = (init_temp * np.ones((N-1)*N)).reshape((N-1),N)
	botNodes = np.pad(botNodes, (1,1), 'constant', constant_values=(bound_temp,bound_temp))[0:N][:]

	midNodes = np.copy(botNodes[N-1][:])

	topNodes = (init_temp * np.ones(N*N)).reshape(N,N)
	topNodes = np.pad(topNodes, (1,1), 'constant', constant_values=(bound_temp,bound_temp))[0:N+1][:]

	print(topNodes)
	doPlottingLeft(ax, midNodes, topNodes, t, maxT)
	doPlottingRight(ax2, midNodes, topNodes, t, maxT)