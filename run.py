from thunderflow import *

s = sim(10)
num = 200
steps = [s.next() for _ in xrange(num)]
np.save('results.npy', steps)
