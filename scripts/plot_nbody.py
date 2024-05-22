import matplotlib.pyplot as plt
import numpy as np 
import sys, os 
import imageio

def read_ascii_file(filename):
    data = {'timestep': [], 'time': [], 'id': [], 
            'm': [], 'r': [],
            'x': [], 'y': [], 'z': [], 
            'vx': [], 'vy': [], 'vz': [],
            'ax': [], 'ay': [], 'az': [],
            }
    # load data 
    data['timestep'], data['id'] = np.loadtxt(filename, dtype = int, usecols = [0,2], unpack=True)
    data['time'], data['m'], data['r'], data['x'], data['y'], data['z'], data['vx'], data['vy'], data['vz'], data['ax'], data['ay'], data['az'] = np.loadtxt(filename, dtype = np.float64, usecols = [1,3,4,5,6,7,8,9,10,11,12,13], unpack=True)
    nparts = np.max(data['id'])
    nsteps = np.max(data['timestep'])+1
    return nparts,nsteps,data

def plot_scatter(nsteps, nparts, data, 
                 cleanplots = True,
                 basefilename = 'plots/nbody-scatter-plot', 
                 ):
    cmap = plt.get_cmap('viridis')
    sample_values = np.linspace(0, 1, nsteps)  # Example: 10 sample values ranging from 0 to 1
    sample_colors = cmap(sample_values)
    fig = plt.figure(figsize=(10, 10))
    alphaval = 1.0/(nparts**(1.0/6.0))
    print('Generating frames...')
    xlim = [np.min(data['x']),np.max(data['x'])]
    ylim = [np.min(data['y']),np.max(data['y'])]
    zlim = [np.min(data['z']),np.max(data['z'])]
    slim = [np.min(data['m']),np.max(data['m'])]
    for i in range(nsteps):
        c = sample_colors[i]
        c = 'Blue'
        gs = plt.GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        indices = np.arange(i*nparts,(i+1)*nparts,1)
        
        #indices = np.array(np.where(data['timestep']==i+1)[0], dtype=np.int32)
        xdata = data['x'][indices]
        ydata = data['y'][indices]
        sdata = (data['m'][indices]-slim[0])/(slim[1]-slim[0])*10.0
        ax1.scatter(xdata, ydata, s=sdata, marker='o', facecolor=c, ec=c, alpha=alphaval)
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.annotate('X-Y', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='black')
        xdata = data['x'][indices]
        ydata = data['z'][indices]
        ax3.scatter(xdata, ydata, s=sdata, marker='o',facecolor=c, ec=c, alpha=alphaval)
        ax3.set_xlim(xlim)
        ax3.set_ylim(zlim)
        ax3.annotate('X-Z', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='black')
        xdata = data['z'][indices]
        ydata = data['y'][indices]
        ax2.scatter(xdata, ydata, s=sdata, marker='o',facecolor=c, ec=c, alpha=alphaval)
        ax2.set_xlim(ylim)
        ax2.set_ylim(zlim)
        ax2.annotate('Y-Z', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='black')
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(basefilename + f'.t-{i}.png')
        fig.clear()
    print('Produced frames, generating movies')
    # Create a GIF and movie from the images
    images = []
    for i in range(nsteps):
        images.append(imageio.v3.imread(basefilename + f'.t-{i}.png'))
    imageio.v3.imwrite(basefilename + '.gif', images)
    imageio.v3.imwrite(basefilename + '.mp4', images, fps=5)
    if cleanplots:
        for i in range(nsteps):
            os.remove(basefilename + f'.t-{i}.png')
    print('Done')

def plot_density(nsteps, nparts, data, 
                 nbins = 25,
                 cleanplots = True, 
                 basefilename = 'plots/nbody-density-plot',
                 ):
    fig = plt.figure(figsize=(10, 10))
    print('Generating frames...')
    xlim = [np.min(data['x']),np.max(data['x'])]
    ylim = [np.min(data['y']),np.max(data['y'])]
    zlim = [np.min(data['z']),np.max(data['z'])]
    for i in range(nsteps):
        gs = plt.GridSpec(2, 2, figure=fig)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, 0])
        indices = np.arange(i*nparts,(i+1)*nparts,1)
        
        #indices = np.array(np.where(data['timestep']==i+1)[0], dtype=np.int32)
        xdata = data['x'][indices]
        ydata = data['y'][indices]
        ax1.hexbin(xdata, ydata, gridsize=nbins, cmap='viridis')
        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        ax1.annotate('X-Y', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='white')
        xdata = data['x'][indices]
        ydata = data['z'][indices]
        ax3.hexbin(xdata, ydata, gridsize=nbins, cmap='viridis')
        ax3.set_xlim(xlim)
        ax3.set_ylim(zlim)
        ax3.annotate('X-Z', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='white')
        xdata = data['z'][indices]
        ydata = data['y'][indices]
        ax2.hexbin(xdata, ydata, gridsize=nbins, cmap='viridis')
        ax2.set_xlim(ylim)
        ax2.set_ylim(zlim)
        ax2.annotate('Y-Z', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='white')
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.savefig(basefilename + f'.t-{i}.png')
        fig.clear()
    print('Produced frames, generating movies')
    # Create a GIF and movie from the images
    images = []
    for i in range(nsteps):
        images.append(imageio.v3.imread(basefilename + f'.t-{i}.png'))
    imageio.v3.imwrite(basefilename + '.gif', images)
    imageio.v3.imwrite(basefilename + '.mp4', images, fps=1)
    if cleanplots:
        for i in range(nsteps):
            os.remove(basefilename + f'.t-{i}.png')
    print('Done')

def plot_orbits(nparts, data, 
                filename = 'plots/nbody-orbit-plot.png',
                subsample = 0.01):
    xlim = [np.min(data['x']),np.max(data['x'])]
    ylim = [np.min(data['y']),np.max(data['y'])]
    zlim = [np.min(data['z']),np.max(data['z'])]
    fig = plt.figure(figsize=(10, 10))
    gs = plt.GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    #ax4 = fig.add_subplot(gs[1, 1])
    colors = ['Red', 'Blue', 'Green', 'Orange', 'Purple', 'LightGreen', 'Cyan', 'Magenta', 'Grey']

    if (nparts > 10):
        nsample = np.random.choice(np.arange(nparts), np.int64(np.ceil(subsample * nparts)))
    else:
        nsample = np.arange(nparts)
    for i in nsample:
        c = colors[i%len(colors)]
        indices = np.array(np.where(data['id']==i+1)[0], dtype=np.int32)
        xdata = data['x'][indices]
        ydata = data['y'][indices]
        ax1.plot(xdata, ydata, linewidth=1, color=c, alpha=0.5, zorder=1)
        ax1.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5)
        ax1.scatter(xdata[0], ydata[0], linewidth=2, marker='o', facecolor='None', ec=c)
        ax1.scatter(xdata[-1], ydata[-1], linewidth=2, marker='*',facecolor='None', ec=c)
        xdata = data['x'][indices]
        ydata = data['z'][indices]
        ax3.plot(xdata, ydata, linewidth=1, color=c, alpha=0.5, zorder=1)
        ax3.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5)
        ax3.scatter(xdata[0], ydata[0], linewidth=2, marker='o', facecolor='None', ec=c)
        ax3.scatter(xdata[-1], ydata[-1], linewidth=2, marker='*',facecolor='None', ec=c)
        xdata = data['z'][indices]
        ydata = data['y'][indices]
        ax2.plot(xdata, ydata, linewidth=1, color=c, alpha=0.5, zorder=1)
        ax2.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5)
        ax2.scatter(xdata[0], ydata[0], linewidth=2, marker='o', facecolor='None', ec=c)
        ax2.scatter(xdata[-1], ydata[-1], linewidth=2, marker='*',facecolor='None', ec=c)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(filename)

def plot_phase(nsteps, nparts, data, 
                filename = 'plots/nbody-phase-plot',
                subsample = 0.01):
    xlim = [np.min(data['x']),np.max(data['x'])]
    ylim = [np.min(data['y']),np.max(data['y'])]
    zlim = [np.min(data['z']),np.max(data['z'])]
    fig = plt.figure(figsize=(10, 15))
    gs = plt.GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    ax5 = fig.add_subplot(gs[2, 0])
    ax6 = fig.add_subplot(gs[2, 1])
    colors = ['Red', 'Blue', 'Green', 'Orange', 'Purple', 'LightGreen', 'Cyan', 'Magenta', 'Grey']

    cm = np.zeros([nsteps,3])
    cmvel = np.zeros([nsteps,3])
    for i in np.arange(nsteps):
        indices = np.arange(i*nparts,(i+1)*nparts,1)
        mtot = np.sum(data['m'][indices])
        cm[i] = np.array([np.sum(data['m'][indices]*data['x'][indices]), 
                          np.sum(data['m'][indices]*data['y'][indices]), 
                          np.sum(data['m'][indices]*data['z'][indices])])/mtot
        cmvel[i] = np.array([np.sum(data['m'][indices]*data['vx'][indices]), 
                          np.sum(data['m'][indices]*data['vy'][indices]), 
                          np.sum(data['m'][indices]*data['vz'][indices])])/mtot
    cm = np.transpose(cm) 
    cmvel = np.transpose(cmvel)
    if (nparts > 10):
        nsample = np.random.choice(np.arange(nparts), np.int64(np.ceil(subsample * nparts)))
    else:
        nsample = np.arange(nparts)
    for i in nsample:
        c = colors[i%len(colors)]
        indices = np.array(np.where(data['id']==i+1)[0], dtype=np.int32)
        xdata = data['x'][indices]-cm[0]
        ydata = data['y'][indices]-cm[1]
        ax1.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax1.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax1.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax1.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)
        xdata = data['vx'][indices]-cmvel[0]
        ydata = data['vy'][indices]-cmvel[1]
        ax2.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax2.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax2.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax2.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)

        xdata = data['x'][indices]-cm[0]
        ydata = data['z'][indices]-cm[2]
        ax3.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax3.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax3.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax3.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)
        xdata = data['vx'][indices]-cmvel[0]
        ydata = data['vz'][indices]-cmvel[2]
        ax4.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax4.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax4.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax4.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)


        xdata = data['y'][indices]-cm[1]
        ydata = data['z'][indices]-cm[2]
        ax5.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax5.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax5.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax5.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)
        xdata = data['vy'][indices]-cmvel[1]
        ydata = data['vz'][indices]-cmvel[2]
        ax6.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax6.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax6.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax6.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)
        fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig(filename+'.png')

    fig = plt.figure(figsize=(10, 10))
    gs = plt.GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    for i in nsample:
        c = colors[i%len(colors)]
        indices = np.array(np.where(data['id']==i+1)[0], dtype=np.int32)
        posdata = np.array([data['x'][indices]-cm[0],data['y'][indices]-cm[1],data['z'][indices]-cm[2]])
        veldata = np.array([data['vx'][indices]-cmvel[0],data['vy'][indices]-cmvel[1],data['vz'][indices]-cmvel[2]])
        xdata = np.sqrt(posdata[0]*posdata[0]+posdata[1]*posdata[1]+posdata[2]*posdata[2])
        ydata = (posdata[0]*veldata[0]+posdata[1]*veldata[1]+posdata[2]*veldata[2])/xdata
        ttot = data['time'][indices[-1]]
        tdata = data['time'][indices]/ttot  # Myr
        ax2.plot(xdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax2.plot(xdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax2.scatter(xdata[0], ydata[0], marker='o', facecolor='None', ec=c)
        ax2.scatter(xdata[-1], ydata[-1], marker='*',facecolor='None', ec=c)
        ax2.annotate(r'$R-V_{r}$', xy=(0.5, 0.9), xycoords='axes fraction', fontsize=20, color='black')
        ax2.yaxis.tick_right()  # Move ticks to the right side
        ax2.yaxis.set_label_position('right')

        ax1.plot(tdata, xdata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax1.plot(tdata, xdata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax1.set_ylabel('Radial Distance from CM [pc]')

        ax3.plot(tdata, ydata, linewidth=2, color=c, alpha=0.25, zorder=1)
        ax3.plot(tdata, ydata, linewidth=0, color=c, marker='o', markersize=2, alpha=0.5, zorder=2)
        ax3.set_ylabel('Radial Velocity from CMVel [km/s]')
        ax3.set_xlabel(r'Normalized Time [$t_f={:.2f}~$Myr]'.format(ttot/3.154e13))
        fig.subplots_adjust(wspace=0, hspace=0)

    fig.savefig(filename+'-rad.png')

filename = sys.argv[1]  # Replace 'data.txt' with the actual filename
plottype = 'phase' 
if (len(sys.argv)>=3): 
    plottype=sys.argv[2]
nparts, nsteps, data = read_ascii_file(filename)
print('Plotting for ', nparts, 'particles over ', nsteps, 'producing',plottype)
if plottype == 'orbit':
    plot_orbits(nparts, data)
elif plottype == 'scatter':
    plot_scatter(nsteps, nparts, data, True)
elif plottype == 'density':
    plot_density(nsteps, nparts, data, 25, True)
elif plottype == 'phase':
    plot_phase(nsteps, nparts, data,)
