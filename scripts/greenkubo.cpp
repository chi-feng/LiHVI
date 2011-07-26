#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

struct atom {
	float q;
	float vx;
	float vy;
	float vz;
};

float GKConductivity(vector<atom> atoms, vector<atom> t0, int natom) {
	// $\sum_{i=1}^N\sum_{j=1}^N\langle q_iq_j\vec{v}_i(t)\cdot\vec{v}_j(0)\rangle$
	float J = 0;
	float dot = 0;
	for (int i = 0; i < natom; i++) { 
		for (int j = 0; j < natom; j++) {
			dot = atoms[i].vx * t0[j].vx + atoms[i].vy * t0[j].vy + atoms[i].vz * t0[j].vz;
			J += atoms[i].q * atoms[j].q * dot;
		}
	}
	return J;
}


int main(int argc, char *argv[]) {

	if (argc < 4) {
		cout << "Usage: GKConductivity ifname ofname #steps" << endl;
		return 1;
	}

	char *ifname = argv[1];
	char *ofname = argv[2];
	int timesteps = atoi(argv[3]);

	ifstream file(ifname);
	ofstream outfile (ofname);
  	
	// max number of characters in a line
	int maxchar = 512;

	vector<atom> t0; // atoms at t = 0
	int timestep = 0;
	int natom = 0;
	string line;

	for (int t = 0; t < timesteps; t++) {
		
		vector<atom> atoms;
		
		// skip 1 line to get to timestep
		getline(file, line);
		getline(file, line);
		timestep = atoi(line.c_str());
		printf("Timestep = %d\n", timestep);

		// skip 1 line to get natom
		getline(file, line);
		getline(file, line);
		natom = atoi(line.c_str());

		// skip 5 lines to get to data
		for (int skip = 0; skip < 5; skip++) {
			getline(file, line);
		}

		// ITEM: ATOMS id type q spin eradius x y z vx vy vz ervel 
		for (int i = 0; i < natom; i++) {
			
			int id, type, spin;
			float q, eradius, x, y, z, vx, vy, vz, ervel;
			
			file >> id >> type >> q >> spin >> eradius >> x >> y >> z >> vx >> vy >> vz >> ervel;
			
			atom a = {q, vx, vy, vz};
			atoms.push_back(a);
			
			// store t = 0 in t0_array
			if (timestep == 0) {
				t0 = atoms;
			}
		}

		// compute J for this timestep
		float J = GKConductivity(atom_array, t0_array, natom);
		outfile << timestep << ' ' << J << endl;

		// get rid of extra newline
		getline(file, line);
	}

	file.close();
	outfile.close();
}