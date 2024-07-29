import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, FancyArrowPatch

scaling = 0.5
offset = (1-scaling)/(2*scaling)
bbox_props = {
    'boxstyle': 'square, pad=0.1',
    'color': 'white'
}

def annotate_connection(fig, axs, indexA, indexB, tstart, tend, 
                        connectStart=True, arrowB=True, label=None):
    arrowIndex = indexB if arrowB else indexA
    if indexB < indexA:
        yA, yB = (0.1, 0.9)
    else:
        yA, yB = (0.9, 0.1)
    if connectStart:
        tconnector, t2 = (tstart, tend)
    else:
        tconnector, t2 = (tend, tstart)
    fig.add_artist(
        ConnectionPatch(xyA=(tconnector,yA), xyB=(tconnector,yB), 
                        coordsA='data', coordsB='data',
                        axesA=axs[indexA], axesB=axs[indexB],
                        color='gray')
    )
    axs[arrowIndex].add_artist(
        FancyArrowPatch(posA=(tstart,0.5), posB=(tend, 0.5), 
                        arrowstyle='<|-|>', facecolor='black',
                        mutation_scale=15, mutation_aspect=0.5)
    )
    axs[arrowIndex].add_artist(
        FancyArrowPatch(posA=(t2,0.9), posB=(t2, 0.1), 
                        arrowstyle='-', color='gray',
                        mutation_scale=15, mutation_aspect=0.5)
    )
    if label:
        axs[arrowIndex].text((tstart+tend)/2, 0.5, label, bbox=bbox_props,
                             fontfamily='monospace', fontsize='small',
                             ha='center', va='center', backgroundcolor='white',
                             )

time = np.arange(25.1, step=0.1, )
trig = np.zeros(time.shape)
trig_s = trig.copy()
trig[time < 0.1] = 1
trig_s[time < 9] = 1
trig_l = np.zeros(time.shape)
trig_l[(time >= 6) * (time < 106)] = 1
laser = np.ones(time.shape)
laser[(time < 8) + (time >= 18)] = 0
spark = np.zeros(time.shape)
spark[(time >= 9) * (time < 9.5)] = 1
spectro = np.zeros(time.shape)
spectro[(time >= 11) * (time < 15)] = 1
camera = np.zeros(time.shape)
camera[(time > 0)] = 1

signals = {
    'trig_main': trig,
    'trig_spark': trig_s,
    'trig_laser': trig_l,
    'laser': laser,
    'spark': spark,
    # 'camera': camera,
    'spectro': spectro
}

fig, axs = plt.subplots(len(signals), 1, sharex=True, figsize=(5.8, 3.5))
fig.subplots_adjust(hspace=0)
for index, key in enumerate(signals):
    axs[index].plot(time, signals[key], 'k-')
    axs[index].set_ylim(0-offset, 1+offset)
    axs[index].set_ylabel(key, rotation='horizontal', fontfamily='monospace', 
                          fontsize='small', ha='right')
    axs[index].set_yticks([])
    axs[index].set_xticks([])
    axs[index].spines[:].set_visible(False)

annotate_connection(fig, axs, 0, 2, 0, 6, label='6 ms')
annotate_connection(fig, axs, 2, 5, 6, 11, label='5 ms')
annotate_connection(fig, axs, 2, 3, 6, 8)
fig.add_artist(ConnectionPatch(xyA=(9,0.9), xyB=(9,0.1), 
                        coordsA='data', coordsB='data',
                        axesA=axs[1], axesB=axs[4],
                        alpha=0.3))

# fig.add_artist(
#     ConnectionPatch(xyA=(0,0.9), xyB=(0,0.1), 
#                     coordsA='data', coordsB='data',
#                     axesA=axs[0], axesB=axs[2],
#                     alpha=0.3)
# )
# axs[2].add_artist(
#     FancyArrowPatch(posA=(0,0.5), posB=(79, 0.5), 
#                     arrowstyle='<|-|>', facecolor='black',
#                     mutation_scale=15, mutation_aspect=0.5)
# )

axs[-1].spines.bottom.set_visible(True)
axs[-1].spines.bottom.set_bounds(time.min(), time.max())
ticklabels = np.linspace(time[0], time[-1], 6)
axs[-1].set_xticks(ticklabels, labels=ticklabels, fontfamily='monospace',
                   fontsize='small')

plt.xlabel('Time [ms]', fontfamily='monospace', fontsize='small')
plt.show()