import MDAnalysis as mda
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
from sys import argv

try:
    resname = str(argv[1]).upper()
except:
    print("Usage: $python pdb-survey_pyranose_sugar_dihedrals.py <resname>")
    raise

def pyranose_sugar_dihedrals_histogram(u, resname, resid, segid="SYSTEM", skip=10):
    u = u
    #print("...Scanning for backbone dihedrals")
    bb = u.select_atoms("segid %s and resname %s and resid %d and name C1 C2 C3 C4 C5 O5"%(segid, resname, resid))
    D = []
    for ts in u.trajectory[0::skip]:
        c1c2c3c4 = bb.select_atoms("name C1", "name C2", "name C3", "name C4")
        c2c3c4c5 = bb.select_atoms("name C2", "name C3", "name C4", "name C5")
        c3c4c5o5 = bb.select_atoms("name C3", "name C4", "name C5", "name O5")
        c4c5o5c1 = bb.select_atoms("name C4", "name C5", "name O5", "name C1")
        c5o5c1c2 = bb.select_atoms("name C5", "name O5", "name C1", "name C2")
        o5c1c2c3 = bb.select_atoms("name O5", "name C1", "name C2", "name C3")
        a1 = c1c2c3c4.dihedral.value() % 360
        a2 = c2c3c4c5.dihedral.value() % 360
        a3 = c3c4c5o5.dihedral.value() % 360
        a4 = c4c5o5c1.dihedral.value() % 360
        a5 = c5o5c1c2.dihedral.value() % 360
        a6 = o5c1c2c3.dihedral.value() % 360
        D.append((a1, a2, a3, a4, a5, a6))
    return D

pdbs = sorted(glob.glob("./pdbs-with-%s/*.pdb"%(resname)))
f = open("pyranose_sugar_dihedrals-%s.txt"%(resname),"w")
f.write("%10s %10s %10s %10s %10s %10s %10s %10s %10s\n"%('PDB', 'segid', 'resid', 'c1c2c3c4', 'c2c3c4c5', 'c3c4c5o5', 'c4c5o5c1', 'c5o5c1c2', 'o5c1c2c3'))
for i,pdb in enumerate(pdbs):
    #print(pdb[:-4])
    try:
        u = mda.Universe(pdb)
    except:
        continue
    sele1 = u.select_atoms('resname %s and name C1 C2 C3 C4 C5 O5'%(resname))
    for segid in sele1.segments.segids:
        sele2 = sele1.select_atoms('segid %s'%(segid))
        for resid in sele2.residues.resids:
            try:
                D = pyranose_sugar_dihedrals_histogram(u, resname, resid, segid=segid, skip=1)
                print(round(100*(i+1)/len(pdbs),2), pdb[:-4][-4:], segid, resid, D[0])
                f.write("%10s %10s %10s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n"%(pdb[:-4][-4:],segid,resid,D[0][0],D[0][1],D[0][2],D[0][3],D[0][4],D[0][5]))
            except:
                continue

f.close()
D = np.loadtxt("pyranose_sugar_dihedrals-%s.txt"%(resname), usecols=(3,4,5,6,7,8))
histogram = plt.figure(figsize=(9,15))
plt.subplots_adjust(hspace=0.5)
plt.subplot(611)
plt.title('c1-c2-c3-c4')
plt.xlim(0,360)
plt.hist(D[:,0], bins=180, density=True, range=(0,360))
plt.subplot(612)
plt.title('c2-c3-c4-c5')
plt.xlim(0,360)
plt.hist(D[:,1], bins=180, density=True, range=(0,360))
plt.subplot(613)
plt.title('c3-c4-c5-o5')
plt.xlim(0,360)
plt.hist(D[:,2], bins=180, density=True, range=(0,360))
plt.subplot(614)
plt.title('c4-c5-o5-c1')
plt.xlim(0,360)
plt.hist(D[:,3], bins=180, density=True, range=(0,360))
plt.subplot(615)
plt.title('c5-o5-c1-c2')
plt.xlim(0,360)
plt.hist(D[:,4], bins=180, density=True, range=(0,360))
plt.subplot(616)
plt.title('o5-c1-c2-c3')
plt.xlim(0,360)
plt.hist(D[:,5], bins=180, density=True, range=(0,360))
histogram.savefig("pyranose_sugar_dihedrals-%s.png"%(resname), dpi=300)
