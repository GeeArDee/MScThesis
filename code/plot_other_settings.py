def plot_other_settings(ax,legend_on,box_on,grid_on):
    #set ticks:
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1,direction='in')
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=2)

    if box_on:
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
    else:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    #set grid: 
    if grid_on:
        ax.grid(color='k', linestyle='--', linewidth=0.1,dashes=(60,20))

    #set legend
    if legend_on:
        legend=ax.legend(fontsize='x-small',frameon=False,facecolor='white', framealpha=1,loc='best',ncol=3)
        frame = legend.get_frame()
        frame.set_linewidth(0)