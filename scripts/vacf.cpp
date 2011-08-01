#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


using namespace std;

struct atom {
    float v0[3];
    vector<float> peratomvacf;    
};

void progress(int n, int N) {
    int inc = N / 100;
    if (n % inc == 0) {
        printf("%5d /%5d [%3d%%]\r", n+1, N, 100 * n / N + 1);
		fflush(stdout);
    }
}

/* $\sum_{t=t_0}^{t_{n+1}} v(0)\cdot v(t_{n+1}) = v(0)\cdot v(t_{n+1}) + \sum_{t=t_0}^{t_n} v(0)\cdot v(t_n)$ */ 
void add_timestep(atom &a, float v[]) {
    float vacf = (a.peratomvacf.size() > 0) ? a.peratomvacf.back() : 0.0;
    vacf += a.v0[0] * v[0];
    vacf += a.v0[1] * v[1];
    vacf += a.v0[2] * v[2];
    a.peratomvacf.push_back(vacf);
}

vector<float> compute_ensemble_vacf(vector<atom> atoms) {
    int natom = atoms.size();
    int timesteps = atoms[0].peratomvacf.size();
    vector<float> vacf;
    for (int t = 0; t < timesteps; t++) {
        float sum = 0.0;
        for (int i = 0; i < natom; i++) {
            float v0 = 0.0;
            for (int j = 0; j < 3; j++) {
                v0 += atoms[i].v0[j] * atoms[i].v0[j];
            }
            sum += atoms[i].peratomvacf[t];
        }
        vacf.push_back(sum/(t+1));
    }
    return vacf;
}

int get_natom(char *ifname) {
	ifstream infile(ifname);
    string line;
	for (int i = 0; i < 4; i++) {
		getline(infile, line);
    }
	int natom = atoi(line.c_str());
    infile.close();
	return natom;
}

vector<int> get_timesteps(char *ifname, int natom) {
    vector<int> timesteps;
	ifstream infile(ifname);
    string line;
    while (!infile.eof()) {
        getline(infile, line);
        getline(infile, line);
        int timestep = atoi(line.c_str());
		timesteps.push_back(timestep);
		for (int i = 0; i < 7 + natom; i++) {
			getline(infile, line);
		}
	}
    infile.close();
	return timesteps;
}
    
int main(int argc, char *argv[]) {

	if (argc < 2) {
		cout << "Usage: vacf ifname.trj" << endl;
		return 1;
	}

	char *ifname = argv[1];
	int len = strlen(argv[1]);
	char ofname[len + 5];
	strcpy(ofname, ifname);
	strcat(ofname, ".vacf");

	ifstream infile(ifname);
	ofstream outfile(ofname);
 	
  	int natom = get_natom(ifname);
  	vector<int> timesteps = get_timesteps(ifname, natom);
  	  	
    vector<atom> atoms;
    for (int i = 0; i < natom; i++) {
        atom a;
        atoms.push_back(a);
    }
    
    cout << "Reading frames from " << ifname << endl;
    
    int M = timesteps.size();
    string line;
    for (int t = 0; t < M; t++) {
        progress(t, M);
		for (int i = 0; i < 9; i++) {
			getline(infile, line);
		}
		// ITEM: ATOMS id type q spin eradius x y z vx vy vz ervel 
		for (int i = 0; i < natom; i++) {
			int id, type, spin;
			float q, eradius, x, y, z, vx, vy, vz, ervel;
			infile >> id >> type >> q >> spin >> eradius >> x >> y >> z >> vx >> vy >> vz >> ervel;
            float v[] = {vx, vy, vz};
            if (timesteps[t] == 0) {
                atoms[i].v0[0] = vx;
                atoms[i].v0[1] = vy;
                atoms[i].v0[2] = vz;
            }
			add_timestep(atoms[i], v);
		}
		getline(infile, line);
    }
    
    cout << "\nComputing ensemble VACF...";
  	vector<float> vacf = compute_ensemble_vacf(atoms);
  	cout << "Done" << endl;
 	
  	cout << "Writing " << timesteps.size() << " frames to " << ofname << "...";
  	for (int i = 0; i < vacf.size(); i++) {
		outfile << timesteps[i] << ' ' << vacf[i] << endl;
  	}
  	cout << "Done" << endl;
  
	infile.close();
	outfile.close();
}
