#!/usr/bin/python

import argparse
import random

def printf(format, *args): print format % args,

def LiSolid(nx, ny, nz):
    '''
    Generates a lithium solid with FCC nuclear positions 
    and interstitial electron positions.
    '''
    L = 4.419688383648  # eFF optimized (8.32 expt)
    # This part changes for different lattices
    xunit = (0, 0.5, 0, 0.5, 0.5, 0, 0, 0.5)
    yunit = (0, 0.5, 0.5, 0, 0, 0.5, 0, 0.5)
    zunit = (0, 0, 0.5, 0.5, 0, 0, 0.5, 0.5)

    r_elec = 0.3704240743 # 1s
    r_elec2 = 1.587531747 # 2s

    Lx, Ly, Lz = L, L, L

    xnuc, ynuc, znuc = [], [], []
    
    for z in xrange(nz):
        for y in xrange(ny):
            for x in xrange(nx):
                for i in xrange(len(xunit)):
                    xnuc.append(x * Lx + xunit[i] * L + 0.2645886245 * L)
                    ynuc.append(y * Ly + yunit[i] * L + 0.2645886245 * L)                    
                    znuc.append(z * Lz + zunit[i] * L + 0.2645886245 * L)             
                           
    numnuc = len(xnuc)

    # length of supercell
    xlo, xhi = 0, Lx * nx
    ylo, yhi = 0, Ly * ny
    zlo, zhi = 0, Lz * nz

    atoms = [] # [type q spin eradius x y z]
    for i in range(0, numnuc, 8):
        for j in range(4):
            # nuclei
            atoms.append([1, 3., 0, 0., xnuc[i+j], ynuc[i+j], znuc[i+j]])
            # core electrons
            atoms.append([2, 0., 1, r_elec, xnuc[i+j], ynuc[i+j], znuc[i+j]])
            atoms.append([2, 0., -1, r_elec, xnuc[i+j], ynuc[i+j], znuc[i+j]])
        # valence electrons
        spin = 1 if random.random() < 0.5 else -1
        atoms.append([2, 0.,  spin, r_elec2, xnuc[i+4], ynuc[i+4], znuc[i+4]])
        atoms.append([2, 0., -spin, r_elec2, xnuc[i+5], ynuc[i+5], znuc[i+5]])
        spin = 1 if random.random() < 0.5 else -1
        atoms.append([2, 0.,  spin, r_elec2, xnuc[i+6], ynuc[i+6], znuc[i+6]])
        atoms.append([2, 0., -spin, r_elec2, xnuc[i+7], ynuc[i+7], znuc[i+7]])

    return atoms, (xlo, xhi), (ylo, yhi), (zlo, zhi)

def write_data(cube_size, slab_size, slab_thickness, separation, filename, seed=None):
    '''
    Writes the lammps data file with a Li FCC cube of size $cube_size
    and a slab of thickness $slab_thickness separated in the z-direction
    by $separation given in angstroms to file $filename
    '''
    # set seed, seed=None is default
    random.seed(seed)
    
    # generate geometry
    cube = LiSolid(cube_size, cube_size, cube_size)
    slab = LiSolid(slab_size, slab_size, slab_thickness)
    natom = len(cube[0]) + len(slab[0])

    # xhi - xlo, yhi - ylo
    cube_dimensions = (cube[1][1] - cube[1][0], cube[2][1] - cube[2][0], cube[3][1] - cube[3][0])
    slab_dimensions = (slab[1][1] - slab[1][0], slab[2][1] - slab[2][0], slab[3][1] - slab[3][0])

    # move cube (x,y) to center of slab and create separation
    xshift = slab_dimensions[0] / 2.0 - cube_dimensions[0] / 2.0
    yshift = slab_dimensions[1] / 2.0 - cube_dimensions[1] / 2.0
    zshift = slab_dimensions[2] + separation # + cube_dimensions[2] / 2.0
    for i in range(len(cube[0])):
        cube[0][i][4] += xshift
        cube[0][i][5] += yshift
        cube[0][i][6] += zshift
    
    # sort slab atoms by z values to fix blayer
    slab[0].sort(lambda a, b : cmp(a[6], b[6]))
    
    # find natom in blayer
    zcutoff = slab[0][0][6]
    blayer = 0
    for atom in slab[0]:
        if atom[6] <= zcutoff:
            blayer += 1
    
    # lithium nucleus and electron mass in amu, respectively
    masses = (6.941000, 0.000548579867)
    # define simulation boundaries
    xlo, xhi = 0, slab[2][1]
    ylo, yhi = 0, slab[3][1]
    zlo, zhi = -100.0, 200.0 # note: z should be shrink-wrapped
    
    out = open('%s' % filename, 'w')
    xyz = open('%s.cube.xyz' % filename, 'w')

    out.write('Created by AJB + CF\n\n')
    out.write('%d atoms\n' % natom)
    out.write('2 atom types\n\n')
    out.write('%f %f xlo xhi\n' % (xlo, xhi))
    out.write('%f %f ylo yhi\n' % (ylo, yhi))
    out.write('%f %f zlo zhi\n\n' % (zlo, zhi))
    out.write('Masses\n\n')
    out.write('1 %f\n2 %f\n\n' % (masses))
    out.write('Atoms\n\n')

    xyz.write('%d\n' % (len(cube[0]) + len(slab[0])))
    xyz.write('Li_%d_%d_%d_%d\n' % (cube_size, slab_size, slab_thickness, separation))
    

    idx = 0 
    for atom in slab[0]:
        idx += 1
        out.write('%d %d %f %d %f %s %s %s \n' % 
            (idx, atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6]))
        if atom[0] == 2:
            name = 'N' if atom[2] == 1 else 'O'
            xyz.write('%s %f %f %f\n' % (name, atom[4], atom[5], atom[6]))
        else:
            xyz.write('%s %f %f %f\n' % ('Li', atom[4], atom[5], atom[6]))
    for atom in cube[0]:
        idx += 1
        out.write('%d %d %f %d %f %s %s %s \n' % 
            (idx, atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6]))
        if atom[0] == 2:
            name = 'N' if atom[2] == 1 else 'O'
            xyz.write('%s %f %f %f\n' % (name, atom[4], atom[5], atom[6]))
        else:
            xyz.write('%s %f %f %f\n' % ('Li', atom[4], atom[5], atom[6]))
    
    xyz.close()
    out.close()

    # cubeidx used in lammps input script to set lammps group
    cubeidx = len(slab[0])
    # blayeridx used in lammps input script to fix bottom layer
    # blayeridx = (len(slab[0]) + len(slab[1])) / slab_thickness
    blayeridx = blayer

    return cubeidx, blayeridx


if __name__ == '__main__':

    description = 'Li impact LAMMPS/eFF data preparation'
    epilog = 'Report bugs to feng@ugcs.caltech.edu'
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-c', type=int, help='cube size (unit cells)', required=True)
    parser.add_argument('-s', type=int, help='slab size (unit cells)', required=True)
    parser.add_argument('-t', type=int, help='slab thickness (unit cells)', required=True)
    parser.add_argument('-sep', type=int, help='separation (anstroms)', required=True)
    parser.add_argument('-o', help='output file', required=True)
    parser.add_argument('-seed', type=int, help='random seed', required=False, default=None)

    args = parser.parse_args()

    result = write_data(args.c, args.s, args.t, args.sep, args.o, seed=args.seed)

    print 'cubeidx: %d\nblayeridx: %d' % result




    
