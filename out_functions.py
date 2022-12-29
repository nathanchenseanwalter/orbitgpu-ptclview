import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import readEqdsk
import KS2D
from matplotlib import animation
from scipy.stats.kde import gaussian_kde
import seaborn as sns
from shapely.geometry import Polygon, Point
import geopandas as gpd
from scipy.spatial import cKDTree as KDTree

from glob_vars import *

def KLdivergence(dat):
    """Compute the Kullback-Leibler divergence between two multivariate samples. 
    Parameters
    ----------
    x : 2D array (n,d)
    Samples from distribution P, which typically represents the true
    distribution.
    y : 2D array (m,d)
    Samples from distribution Q, which typically represents the approximate
    distribution.
    Returns
    -------
    out : float
    The estimated Kullback-Leibler divergence D(P||Q).
    References
    ----------
    Pérez-Cruz, F. Kullback-Leibler divergence estimation of
    continuous distributions IEEE International Symposium on Information
    Theory, 2008.
    """
    x = dat[["xgc","zgc"]][:numpoints:skip].to_numpy()
    y = dat[["xflr","zflr"]][:numpoints:skip].to_numpy()
    
    
    # Check the dimensions are consistent
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    
    n,d = x.shape
    m,dy = y.shape
    
    assert(d == dy)
    
    
    # Build a KD tree representation of the samples and find the nearest neighbour
    # of each point in x.
    xtree = KDTree(x)
    ytree = KDTree(y)
    
    # Get the first two nearest neighbours for x, since the closest one is the
    # sample itself.
    r = xtree.query(x, k=2, eps=.01, p=2)[0][:,1]
    s = ytree.query(x, k=1, eps=.01, p=2)[0]
    
    # There is a mistake in the paper. In Eq. 14, the right side misses a negative sign
    # on the first term of the right hand side.
    return -np.log(r/s).sum() * d / n + np.log(m / (n - 1.))

def KS_test(dat):
    arrgc = dat[["xgc","zgc"]][:numpoints:skip].to_numpy()
    arrflr = dat[["xflr","zflr"]][:numpoints:skip].to_numpy()
    print(KS2D.ks2d2s(arrgc,arrflr))
    
def points_out(dat,xval,zval):
    #TODO: implement gc edge check for points clearly stuck on edge
    array_length = len(dat[xval])
    contained_points = np.empty(array_length,dtype=bool)
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    polygon = Polygon(list(zip(x_in,z_in)))
    # plt.plot(*polygon.exterior.xy) # plots polygon
    dat['coords'] = list(zip(dat[xval],dat[zval]))
    dat['coords'] = dat['coords'].apply(Point)
    points = gpd.GeoDataFrame(dat, geometry='coords')
    for i in range(array_length):
        contained_points[i] = polygon.contains(points['coords'][i])
    count = np.count_nonzero(contained_points)
    print(xval + ": " + str((array_length - count)/array_length))
    
def tok_scatter(dat,xval,zval,cm):
    fs=(8,8)
    tok_x_lim = (0,2)
    tok_z_lim = (-2,2)
    
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    x_out = eq_dat_dict['rlim']
    z_out = eq_dat_dict['zlim']
    
    mylabels = ['co-pass c','ctr-pass c','trap c','stag','conf potato','trap potato']
    
    fig, ax = plt.subplots()
    fig.set_size_inches(fs)
    line_in = Line2D(x_in,z_in)
    line_out = Line2D(x_out,z_out)
    
    scatter = plt.scatter(dat[xval],dat[zval],s=5,c=dat.otp.astype('category').cat.codes,cmap=cm)
    plt.scatter(x=-eq_dat_dict['bcentr'],y=0,c='r')
    
    ax.add_line(line_in)
    ax.add_line(line_out)
    
    ax.set_xlabel(xval)
    ax.set_xlim(tok_x_lim)
    ax.set_ylabel(zval)
    ax.set_ylim(tok_z_lim)
    
    plt.legend(handles=scatter.legend_elements()[0],title="species",labels=mylabels)
    plt.show()
    plt.close()
    
def tok_scatter_otp(dat,xval,zval):
    
    plt.rcParams['figure.constrained_layout.use'] = False
    
    fs=(12,8)
    s = 1
    a = 0.1
    tok_x_lim = (0,2)
    tok_z_lim = (-2,2)
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    x_out = eq_dat_dict['rlim']
    z_out = eq_dat_dict['zlim']
    
    fig, ax = plt.subplots(nrows=2,ncols=3,sharex=True, sharey=True)
    
    fig.set_size_inches(fs)
    
    for i in range(2):
        for j in range(3):
            ax[i,j].plot(x_in,z_in,c='C1',alpha=0.3)
            ax[i,j].plot(x_out,z_out,c='C1',alpha=0.3)
            ax[i,j].set_xlim(tok_x_lim)
            ax[i,j].set_ylim(tok_z_lim)
    
    ax[0,0].scatter(dat[dat.otp == 1][xval],dat[dat.otp == 1][zval],s=s,alpha=a)
    ax[0,1].scatter(dat[dat.otp == 3][xval],dat[dat.otp == 3][zval],s=s,alpha=a)
    ax[0,2].scatter(dat[dat.otp == 5][xval],dat[dat.otp == 5][zval],s=s,alpha=a)
    ax[1,0].scatter(dat[dat.otp == 7][xval],dat[dat.otp == 7][zval],s=s,alpha=a)
    ax[1,1].scatter(dat[dat.otp == 8][xval],dat[dat.otp == 8][zval],s=s,alpha=a)
    ax[1,2].scatter(dat[dat.otp == 9][xval],dat[dat.otp == 9][zval],s=s,alpha=a)
    
    ax[0,0].title.set_text('Co-Passing')
    ax[0,1].title.set_text('Counter-Passing')
    ax[0,2].title.set_text('Trapped')
    ax[1,0].title.set_text('Stagnation')
    ax[1,1].title.set_text('Confined Potato')
    ax[1,2].title.set_text('Trapped Lost')
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel("x (m)")
    plt.ylabel("z (m)")
    
    plt.show()
    
def update(num, xs, ys, zs, line):
    line.set_data(xs[0:num],ys[0:num])
    line.set_3d_properties(zs[0:num])
    return line,
    
def traj3Dani(dat,x1,x2,x3):
    xs = dat['x_flr'][:numpoints:skip]
    ys = dat['y_flr'][:numpoints:skip]
    zs = dat['zflr'][:numpoints:skip]
    
    # plt.rcParams['animation.ffmpeg_path'] ='C:\\ffmpeg\\bin\\ffmpeg.exe'
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #plotting red dotted lines with tiny markers
    line, = ax.plot(xs, ys, zs)
    #and on top of it goes a scatter plot with different markersizes
    ax.set_xlabel('x (m)')
    ax.set_xlim3d([-1.5, 1.5])
    ax.set_ylabel('y (m)')
    ax.set_ylim3d([-1.5, 1.5])
    ax.set_zlabel('z (m)')
    ax.set_zlim3d([-1.50, 1.50])
    
    ani = animation.FuncAnimation(fig, update, fargs=(xs, ys, zs, line), frames= 2000, interval=60, blit=False)
    FFwriter = animation.FFMpegWriter(fps=30)
    ani.save('torus.mp4', writer=FFwriter)
    
    plt.show()
    
def traj3D(dat):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #plotting red dotted lines with tiny markers
    # ax.plot(dat['x_flr'][0:numpoints], dat['y_flr'][0:numpoints], dat['zflr'][0:numpoints], linewidth=0.5)
    ax.plot(dat['x_gc'][0:numpoints], dat['y_gc'][0:numpoints], dat['zgc'][0:numpoints], linewidth=1)
    ax.scatter(dat['x_gc'][numpoints], dat['y_gc'][numpoints], dat['zgc'][numpoints], s=2,c='C1')
    #and on top of it goes a scatter plot with different markersizes
    ax.set_xlabel('x (m)')
    ax.set_xlim3d([-1.50, 1.50])
    ax.set_ylabel('y (m)')
    ax.set_ylim3d([-1.50, 1.50])
    ax.set_zlabel('z (m)')
    ax.set_zlim3d([-1.50, 1.50])
    plt.show()
    
def gc_flrdiff(dat):
    delx = dat['xflr'] - dat['xgc']
    delz = dat['zflr'] - dat['zgc']
    fig, ax = plt.subplots()
    plt.plot(dat['time'][0:numpoints],delx[0:numpoints],label='delta x')
    plt.plot(dat['time'][0:numpoints],delz[0:numpoints],label='delta z')
    ax.set_xlabel('time')
    ax.set_ylabel('position diff')
    ax.legend()
    plt.show()
    
def trajtop(dat):
    fig, ax = plt.subplots()
    fig.set_size_inches((5,5))
    plt.plot(dat['x_flr'][0:numpoints],dat['y_flr'][0:numpoints],linewidth=.5,label='flr')
    plt.plot(dat['x_gc'][0:numpoints],dat['y_gc'][0:numpoints],linewidth=.5,label='gc')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    plt.gca().set_aspect('equal')
    ax.legend()
    plt.show()
    
def trajside(dat):
    fig, ax = plt.subplots()
    fig.set_size_inches((5,5))
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    line_in = Line2D(x_in,z_in)
    ax.plot(dat['xflr'][0:numpoints],dat['zflr'][0:numpoints],linewidth=.5,label='flr')
    plt.plot(dat['xgc'][0:numpoints],dat['zgc'][0:numpoints],linewidth=.5,label='gc')
    ax.add_line(line_in)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('z (m)')
    ax.set_xlim([0, 1.75])
    ax.set_ylim([-1.5, 1.5])
    ax.legend()
    plt.show()
    
def heatmap(dat,xtype,ztype,plottype):
    x = dat[xtype]
    y = dat[ztype]
    
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    x_out = eq_dat_dict['rlim']
    z_out = eq_dat_dict['zlim']
    line_in = Line2D(x_in,z_in,color='white')
    line_out = Line2D(x_out,z_out,color='white')
    
    xlim = [0,2]
    ylim = [-2,2]
    
    k = gaussian_kde(np.vstack([x, y]))
    xi, yi = np.mgrid[xlim[0]:xlim[1]:x.size**0.5*1j,ylim[0]:ylim[1]:y.size**0.5*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot()
    
    # alpha=0.5 will make the plots semitransparent
    if(plottype):
        pc = ax.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap='hot')
    else:
        pc = ax.contourf(xi, yi, zi.reshape(xi.shape),cmap='hot')
    ax.add_line(line_in)
    ax.add_line(line_out)
    
    fig.colorbar(pc)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xtype)
    ax.set_ylabel(ztype)
    
    plt.show()
    
def enptchheatmap(dat,xtype,ztype,plottype):
    x = dat[xtype][:numpoints:skip]
    y = dat[ztype][:numpoints:skip]
    
    
    xlim = [-1,1]
    ylim = [y.min(),y.max()]
    
    k = gaussian_kde(np.vstack([x, y]))
    xi, yi = np.mgrid[xlim[0]:xlim[1]:x.size**0.5*1j,ylim[0]:ylim[1]:y.size**0.5*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot()
    
    # alpha=0.5 will make the plots semitransparent
    if(plottype):
        pc = ax.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap='hot')
    else:
        pc = ax.contourf(xi, yi, zi.reshape(xi.shape),cmap='hot')
    
    cbar = fig.colorbar(pc)
    cbar.set_label("probability density function", rotation=270, labelpad=15)
    ax.set_xlim(xlim)
    ax.set_xlabel('pitch')
    ax.set_ylabel('energy (keV)')
    
    plt.show()
    
def snsheatmap(dat,xtype,ztype):
    
    x_in = eq_dat_dict['rbbbs']
    z_in = eq_dat_dict['zbbbs']
    x_out = eq_dat_dict['rlim']
    z_out = eq_dat_dict['zlim']
    line_in = Line2D(x_in,z_in,color='white')
    line_out = Line2D(x_out,z_out,color='white')
    
    xlim = [0,2]
    ylim = [-2,2]
    
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot()
    
    sns.kdeplot(x=dat[xtype][::10],y=dat[ztype][::10],fill=True,thresh=0,levels=100,cmap="mako",cbar=True,ax=ax)
    ax.add_line(line_in)
    ax.add_line(line_out)
    
    ax.set_facecolor('black')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    plt.show()
    
def trajTok(dat,animate):
    dat['x_flr'] = dat['xflr'] * np.cos(dat['zet'])
    dat['y_flr'] = dat['xflr'] * np.sin(dat['zet'])
    dat['x_gc'] = dat['xgc'] * np.cos(dat['zet'])
    dat['y_gc'] = dat['xgc'] * np.sin(dat['zet'])
    if animate == True:
        traj3Dani(dat,'x','y','zflr')
    else:
        traj3D(dat)
        trajtop(dat)

def fixvals(dat):
    norm = ekev / engn
    dat.en = dat.en * norm
    dat.thet = dat.thet % (2 * math.pi) - math.pi
    try:
        dat.rmu = dat.rmu * norm
        dat.zet = dat.zet % (2 * math.pi) - math.pi
    except:
        pass
    dat.xgc = dat.xgc / 100
    dat.zgc = dat.zgc / 100
    dat.xflr = dat.xflr / 100
    dat.zflr = dat.zflr / 100

def checkvals(dat):
    if dat.xgc.max() > 100:
        fixvals(dat)