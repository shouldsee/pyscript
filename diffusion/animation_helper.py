import matplotlib.animation as animation

def plot_partial(frame_num, line,  y_data):
    # NOTE: there is no .set_data() for 3 dim data
    # need to set x y and z data separately
    #n_obs=frame_num*increment+1
    # n_obs=frame_num
    #if x_data is not None: line.set_xdata(x_data[frame_num])
    if y_data is not None: line.set_ydata(y_data[frame_num])
    #print frame_num,y_data[frame_num]
    fig=line.get_figure()
    fig.suptitle('{}'.format(frame_num))
    return line

def plot_time_series(fig1, line, results,frame_interval=5):
    anim=animation.FuncAnimation(fig1, plot_partial, len(results), fargs=(line,results),
                                  interval=frame_interval, blit=False,repeat=True)
    return anim
