import matplotlib.pyplot as plt
import numpy as np 
import sys, os 

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

def check_comparison(nparts1:int, nparts2: int, 
                     nsteps1: int, nsteps2: int, 
                     data1: dict, data2: dict):
    if nparts1 != nparts2:
        print ('Comparing data from two systems containing different number of particles', nparts1, nparts2)
        print('Exiting')
        exit()
    if nsteps1 != nsteps2:
        print ('Comparing data from two systems containing over different number of steps', nsteps1, nsteps2)
        print('Exiting')
        exit()
    # check that ics are the same
    ic1 = np.array([data1['x'][:nparts1],data1['y'][:nparts1],data1['z'][:nparts1],data1['vx'][:nparts1],data1['vy'][:nparts1],data1['vz'][:nparts1]])
    ic2 = np.array([data2['x'][:nparts1],data2['y'][:nparts1],data2['z'][:nparts1],data2['vx'][:nparts1],data2['vy'][:nparts1],data2['vz'][:nparts1]])
    if not np.array_equal(ic1,ic2):
        print ('Comparing data from different ICs')
        print('Exiting')
        exit()


def plot_diff(data1 : dict, 
              data2 : dict, 
              basefilename = 'plots/comparison',
            ):
    fig = plt.figure(figsize=(15, 15))
    gs = plt.GridSpec(3, 3, figure=fig)
    axs = dict()
    axs['x'] = fig.add_subplot(gs[0, 0])
    axs['y'] = fig.add_subplot(gs[1, 0])
    axs['z'] = fig.add_subplot(gs[2, 0])
    axs['vx'] = fig.add_subplot(gs[0, 1])
    axs['vy'] = fig.add_subplot(gs[1, 1])
    axs['vz'] = fig.add_subplot(gs[2, 1])
    axs['ax'] = fig.add_subplot(gs[0, 2])
    axs['ay'] = fig.add_subplot(gs[1, 2])
    axs['az'] = fig.add_subplot(gs[2, 2])
    colors = ['Red', 'Blue', 'Green', 'Orange', 'Purple', 'LightGreen', 'Cyan', 'Magenta', 'Grey']
    nparts = np.max(data1['id'])
    indices = np.array(np.where(data1['id']==1)[0], dtype=np.int32)
    ttot = data1['time'][indices[-1]]
    tdata = data1['time'][indices]/ttot  # Myr
    for v in ['x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az']:
        for i in range(nparts):
            c = colors[i%len(colors)]
            ydata = 0.5*np.sqrt((data1[v][indices]-data2[v][indices])**2.0/(data1[v][indices]+data2[v][indices])**2.0)
            axs[v].plot(tdata, ydata, linewidth=1, color=c, alpha=0.5, zorder=1)
            axs[v].set_ylabel(r'Normalised $|\Delta$'+v+'|')
            axs[v].set_xlabel(r'Normalized Time [$t_f={:.2f}~$Myr]'.format(ttot/3.154e13))
    fig.subplots_adjust(hspace=0)
    fig.savefig(basefilename+'-time.png')

# read file names for comparison
filename1 = sys.argv[1]  
filename2 = sys.argv[2]  
nparts1, nsteps1, data1 = read_ascii_file(filename1)
nparts2, nsteps2, data2 = read_ascii_file(filename2)
check_comparison(nparts1,nparts2,nsteps1,nsteps2,data1,data2)
print('Processing and comparing inputs for ', filename1, ' and ', filename2)
plot_diff(data1, data2)

