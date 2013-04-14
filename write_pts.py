f = open('XY.pts', 'w')
for i in range(100):
    f.write('\n %i %i %06f %06f'%(i+1, 0, 12.5 + 31 * i, -12.5))
f.close()
