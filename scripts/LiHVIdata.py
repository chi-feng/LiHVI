#!/usr/bin/python

import argparse
import random

# easier to translate printf from perl to python
def printf(format, *args): print format % args,

def appendAtom(atoms, atom):
    ''' appends a list [type q spin eradius x y z]
        to a list of dicionaries '''
    atoms.append({'type':atom[0], 'q':atom[1], 'spin':atom[2], 'eradius':atom[3], 
        'x':atom[4], 'y':atom[5], 'z':atom[6]})

def writeAtom(ofstream, atom, idx):
    ''' writes an atom to a lammps data file '''
    out.write('%d %d %f %d %f %s %s %s \n' % 
        (idx, atom['type'], atom['q'], atom['spin'], atom['eradius'], 
        atom['x'], atom['y'], atom['z'])) 

def LiSolid(nx, ny, nz):
    ''' Generates a lithium solid with FCC nuclear positions 
    and interstitial electron positions. Adapted from Li-solid.pl '''
    # eFF optimized lattice constant (8.32 expt)
    L = 4.419688383648 
    Lx, Ly, Lz = L, L, L
    # This part changes for different lattices
    xunit = (0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5)
    yunit = (0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5)
    zunit = (0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5)
    r_elec = 0.3704240743 # 1s radius
    r_elec2 = 1.587531747 # 2s radius
    
    xnuc, ynuc, znuc = [], [], []
    for z in xrange(nz):
        for y in xrange(ny):
            for x in xrange(nx):
                for i in xrange(len(xunit)):
                    xnuc.append(x * Lx + xunit[i] * L + 0.2645886245 * L)
                    ynuc.append(y * Ly + yunit[i] * L + 0.2645886245 * L)                    
                    znuc.append(z * Lz + zunit[i] * L + 0.2645886245 * L)  
    numnuc = len(xnuc)
    
    # bounds of supercell
    xlo, xhi = 0, Lx * nx
    ylo, yhi = 0, Ly * ny
    zlo, zhi = 0, Lz * nz
    atoms = [] # [type q spin eradius x y z]
    for i in range(0, numnuc, 8):
        for j in range(4):
            # nuclei
            appendAtom(atoms, [1, 3.,  0,     0., xnuc[i+j], ynuc[i+j], znuc[i+j]])
            # core electrons
            appendAtom(atoms, [2, 0.,  1, r_elec, xnuc[i+j], ynuc[i+j], znuc[i+j]])
            appendAtom(atoms, [2, 0., -1, r_elec, xnuc[i+j], ynuc[i+j], znuc[i+j]])
        # valence electrons
        spin = 1 if random.random() < 0.5 else -1
        appendAtom(atoms, [2, 0.,  spin, r_elec2, xnuc[i+4], ynuc[i+4], znuc[i+4]])
        appendAtom(atoms, [2, 0., -spin, r_elec2, xnuc[i+5], ynuc[i+5], znuc[i+5]])
        spin = 1 if random.random() < 0.5 else -1
        appendAtom(atoms, [2, 0.,  spin, r_elec2, xnuc[i+6], ynuc[i+6], znuc[i+6]])
        appendAtom(atoms, [2, 0., -spin, r_elec2, xnuc[i+7], ynuc[i+7], znuc[i+7]])

    return {'atoms':atoms, 'bounds':{'x':xhi, 'y':yhi, 'z':zhi}}

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
    natom = len(cube['atoms']) + len(slab['atoms'])

    # move cube (x,y) to center of slab and create separation in z
    for i in range(len(cube[0])):
        cube['atoms'][i]['x'] += slab['bounds']['x'] / 2.0 - cube['bounds']['x'] / 2.0
        cube['atoms'][i]['y'] += slab['bounds']['y'] / 2.0 - cube['bounds']['y'] / 2.0
        cube['atoms'][i]['z'] += slab['bounds']['z'] + separation 
    
    # sort slab atoms by z values to fix blayer
    slab['atoms'].sort(lambda atom1, atom2 : cmp(atom1['z'], atom2['z']))
    cube['atoms'].sort(lambda atom1, atom2 : cmp(atom1['z'], atom2['z']))
    
    # find natom in blayer
    zcutoff = slab['atoms'][0]['z']
    blayeridx = 0
    for atom in slab[0]:
        if atom[6] <= zcutoff:
            blayeridx += 1
    
    # lithium nucleus and electron mass in amu, respectively
    # masses = (6.941000, 0.000548579867)
    masses = (6.941000, 1.0)
    # define simulation boundaries
    xlo, xhi = 0, slab['bounds']['x'] 
    ylo, yhi = 0, slab['bounds']['y'] 
    zlo, zhi = -100.0, 200.0 # note: z should be shrink-wrapped
    
    out = open('%s' % filename, 'w')

    out.write('Created by LiHVIdata.py\n\n')
    out.write('%d atoms\n' % natom)
    out.write('2 atom types\n\n')
    out.write('%f %f xlo xhi\n' % (xlo, xhi))
    out.write('%f %f ylo yhi\n' % (ylo, yhi))
    out.write('%f %f zlo zhi\n\n' % (zlo, zhi))
    out.write('Masses\n\n1 %f\n2 %f\n\n' % (masses))
    out.write('Atoms\n\n')

    for i in xrange(len(slab['atoms']))
        writeAtom(out, slab['atoms'], i+1)
    for i in xrange(len(cube['atoms']))
        writeAtom(out, cube['atoms'], i+1+len(slab['atoms']))
    
    out.close()

    # cubeidx used in lammps input script to set lammps group
    cubeidx = len(slab['atoms'])
    
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




    
