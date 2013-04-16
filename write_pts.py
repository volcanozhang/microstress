"""
f = open('XY.pts', 'w')
for i in range(100):
    f.write('\n %i %i %06f %06f'%(i+1, 0, 12 + 31 * i, -12))
f.close()
"""
f = open('XY100x100.pts', 'w')
for i in range(101):
    for j in range(101):
        f.write('\n %i %i %06f %06f'%(i*100+j+1, 0, 12.5 + 25 * j, -12.5 - 25 * i))
f.close()
