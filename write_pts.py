f = open('XY.pts', 'w')
for i in range(100):
    f.write('\n %i %i %06f %06f'%(i+1, 0, 15.5, 15.5 + 31.5 * i))
f.close()
