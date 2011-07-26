#!/usr/bin/python

import sys

def convert2xyz(ifname, ofname):

    electron = 'H'
    nuclei = 'Li'

    infile = open(ifname, 'r')          
    outfile = open(ofname, 'w')
    
    timestep = 0
    natom = 0
    frames = 0
    
    for line in infile:
        if line.find("ITEM: TIMESTEP") != -1:
            timestep = int(infile.next())
            frames += 1
        if line.find("ITEM: NUMBER OF ATOMS") != -1:
            natom = int(infile.next())
        if line.find("ITEM: ATOMS") != -1:
            outfile.write('%d\n%s\n' % (natom, ifname))
            for i in xrange(natom):
                tokens = infile.next().split()
                type = 'Li' if tokens[1] == '1' else 'H'
                outfile.write('%s %s %s %s\n' % (type, tokens[5], tokens[6], tokens[7]))
                
    infile.close()
    outfile.close()
    
    return frames

def convert2vel(ifname, ofname):

    infile = open(ifname, 'r')          
    outfile = open(ofname, 'w')
    
    timestep = 0
    natom = 0
    frames = 0
    
    for line in infile:
        if line.find("ITEM: TIMESTEP") != -1:
            timestep = int(infile.next())
            time = float(timestep) *  0.0001 # in fs
            frames += 1
            outfile.write('%d %f\n' % (timestep, time))
        if line.find("ITEM: NUMBER OF ATOMS") != -1:
            natom = int(infile.next())
        if line.find("ITEM: ATOMS") != -1:
            # id type q spin eradius x y z vx vy vz ervel
            #  0    1 2    3       4 5 6 7  8  9 10    11
            for i in xrange(natom):
                tokens = infile.next().split()
                if tokens[1] == '2':
                    outfile.write('%s %s %s\n' % (tokens[8], tokens[9], tokens[10]))
                
    infile.close()
    outfile.close()

if __name__ == '__main__':


    description = 'eFF trajectory format conversion'
    epilog = 'Report bugs to feng@ugcs.caltech.edu'
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-i', type=string, help='input filename', required=True)
    parser.add_argument('-o', type=string, help='output filename', required=True)
    parser.add_argument('-f', type=int, help='number of frames to convert', default=0, required=False)
    args = parser.parse_args()



    if len(sys.argv) < 3:
        print 'incorrect usage'
        sys.exit(0)
        
    ifname = args.i
    ofname = args.o
    frames = args.f

    if ofname.endswith('vel'):
        convert2vel(ifname, ofname)    
    if ofname.endswith('xyz'):
        convert2xyz(ifname, ofname)

    print 'wrote %d frames to %s' % (frames, ofname)
    

