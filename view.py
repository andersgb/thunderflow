from thunderflow import *
import pickle
from matplotlib import pyplot as plt
from matplotlib import animation

resfile = open('results.pkl', 'r')
steps = pickle.load(resfile)

def eos_var(var, q, ind = 0):
    eos_update(q)
    return p[var][ind]

fig = plt.figure()
ax = plt.axes(xlim=(0,N), ylim=(0, 700))
line, = ax.plot([], [], lw=2)

def init():
    line.set_data([], [])
    return line,

def animate(i):
    xs = range(0,N)
    y = [eos_var('state_t', steps[i][x]) for x in xrange(N)]
    line.set_data(xs, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)
anim.save('temperature.mp4', fps=30) #, extra_args=['-vcodec', 'libx264'])

# # velo0 = [steps[x][0][1]/steps[x][0][0] for x in xrange(num)]
# # velo1 = [steps[x][1][1]/steps[x][1][0] for x in xrange(num)]
# # velo9 = [steps[x][9][1]/steps[x][9][0] for x in xrange(num)]
# # velo8 = [steps[x][8][1]/steps[x][8][0] for x in xrange(num)]
# velo0 = [steps[0][x][1]/steps[0][x][0] for x in xrange(N)]
# velo20 = [steps[20][x][1]/steps[20][x][0] for x in xrange(N)]
# ts = map(lambda x: x*dt, range(num))
# # plt.plot(ts, velo0, label='cv0')
# # plt.plot(ts, velo1, label='cv1')
# # plt.plot(ts, velo8, label='cv8')
# # plt.plot(ts, velo9, label='cv9')
# # plt.xlim(0, dt*100)
# p0 = [eos_var('state_p', steps[0][x]) for x in xrange(N)]
# p20 = [eos_var('state_p', steps[20][x]) for x in xrange(N)]
# p180 = [eos_var('state_p', steps[180][x]) for x in xrange(N)]
# plt.plot(p0, label='time 0')
# plt.plot(p20, label='time 20dt')
# plt.plot(p180, label='time 180dt')
# plt.legend(loc='center')
# plt.show()
