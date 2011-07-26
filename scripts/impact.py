#!/usr/bin/python

import argparse
import sys
import os
import time
import LiHVIdata 
from convertlammps import *

if __name__ == '__main__':

    description = 'Li Hypervelocity Impact Script version July 12, 2011'
    epilog = 'Report bugs to feng@ugcs.caltech.edu'
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-v', type=int, help='impactvel in km/s', required=True)
    parser.add_argument('-np', type=int, default=1, help='Number of processors')
    parser.add_argument('-n', default='none', help='Notes for this run')
    parser.add_argument('-cube', type=int, default=4, help='Size of Li cubes')
    parser.add_argument('-slab', type=int, default=12, help='Size of Li slab')
    parser.add_argument('-thick', type=int, default=4, help='Thickness of Li slab')
    parser.add_argument('-sep', type=int, default=3, help='Sep. between cube faces (ang)')
    parser.add_argument('-seed', type=int, default=0, help='Random seed')
    args = parser.parse_args()

    start = time.time()

    # construct lammps input script
    
    basename = 'Li-HVI'
    impactvelkms = args.v
    impactvel = -float(impactvelkms) / 100.0 # km/s -> ang/fs
    notes = args.n
    numprocs = args.np
    if numprocs > 1 and numprocs % 2 != 0:
        print 'want even number of processors'
        sys.exit(0)
    processormapping = ['1 1 1','2 1 1','2 2 1','3 2 1','4 2 1']
    processors = processormapping[int(numprocs / 2)]
    
    # get next index from file 
    
    indexf = open('index', 'r')
    index = int(indexf.readlines()[-1].split(',')[0]) + 1
    indexf.close()
    indexs = '%04d' % index
    impactvels = '%0.2d' % impactvelkms
    name = '%s_%s_%s' % (basename, indexs, impactvels)    
    
    # write data file
    
    seed = args.seed if args.seed != 0 else None
    indices = LiHVIdata.write_data(args.cube, args.slab, args.thick, args.sep, 'data.%s' % name, seed)
    cubeidx = indices[0]
    blayeridx = indices[1] + 1 # + 1 since we have < blayeridx in input script
    print 'lammps data file written to data.%s' % name

    # write input script from template
    
    template = open('templates/in.Li-HVI', 'r')
    output = open('in.%s' % name, 'w')
    output.write('processors          %s\n' % processors)
    output.write('variable            sname index %s\n' % name)
    output.write('variable            impactvel equal %f\n' % impactvel)
    output.write('variable            cubeidx equal %s\n' % cubeidx)
    output.write('variable            blayeridx equal %s\n' % blayeridx)
    
    for line in template:
        output.write(line)
    template.close()
    output.close()
    print 'lammps input script written to in.%s' % name

    # write script to run lammps to sh
    
    output = open('%s.sh' % name, 'w')
    output.write('#!/bin/bash\n')
    if numprocs > 1:
        output.write("nice -10 mpirun -np %d lmp_openmpi -in in.%s\n" % (numprocs, name))
    elif numprocs == 1:
        output.write("nice -10 ./lmp_serial -in in.%s\n" % (name))    

    output.close()

    '''
    # run lammps using appropriate command
    
    if numprocs > 1:
        os.system("nice -10 mpirun -np %d lmp_openmpi -in in.%s" % (numprocs, name))
    elif numprocs == 1:
        os.system("nice -10 ./lmp_serial -in in.%s" % (name))    
    else:
        print 'Invalid numprocs (value %d)' % numprocs
        sys.exit(0)
    '''
    '''
    # convert lammpstrj to xyz viewable in vmd
    
    convert2xyz('%s.impact' % name, '%s.xyz' % name)
    '''

    elapsed = (time.time() - start) / 60 # seconds -> minutes
    print 'time elapsed: %4.1f minutes' % elapsed
    
    # increment index

    indexf = open('index', 'a')
    indexf.write('%s, %s, %4.1f, "%s"\n' % (indexs, impactvelkms, elapsed, notes))


