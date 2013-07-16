from thunderflow import *
import pickle

s = sim(10)
num = 200
steps = [s.next() for _ in xrange(num)]
outfile = open('results.pkl', 'w')
pickle.dump(steps, outfile)
