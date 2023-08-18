

# https://zenn.dev/kaityo256/articles/md_initial_condition 


import numpy as np
import yaml
import random
import sys


def get_lattice_number(L, rho):
    m = np.ceil((L**3 * rho / 4.0)**(1.0 / 3.0))
    return int(m)


def make_fcc_pure(L, rho, Lb):
    m = get_lattice_number(L, rho)
    a = L / m
    ha = a * 0.5
    isin_box = lambda x,y,z: (0.0<=x<Lb[0]) and (0.0<=y<Lb[1]) and (0.0<=z<Lb[2])
    atoms = []
    for i in range(m**3):
        ix = i % m
        iy = (i // m) % m
        iz = i // (m * m)
        x = ix * a
        y = iy * a
        z = iz * a
        if isin_box(x,y,z): atoms.append((x, y, z))
        if isin_box(x+ha,y+ha,z): atoms.append((x + ha, y + ha, z))
        if isin_box(x+ha,y,z+ha): atoms.append((x + ha, y, z + ha))
        if isin_box(x,y+ha,z+ha): atoms.append((x, y + ha, z + ha))
    return atoms


def make_fcc_defect(L, rho, Lb):
    atoms = make_fcc_pure(L, rho, Lb)
    V = Lb[0]*Lb[1]*Lb[2]
    n = int(rho * V)
    return random.sample(atoms, n)

def scaling_momentum(T,vel):
    Kmean = (0.5/N)*np.sum(vel*vel)
    Tobs = Kmean/1.5
    alpha = np.sqrt(T/Tobs)
    vel *= alpha
    return

with open(sys.argv[1],"r+") as f:
    p = yaml.safe_load(f)
    Lb = p["system"]["box_size"]
    V = Lb[0]*Lb[1]*Lb[2]
    L = max(Lb)
    rho = float(p["system"]["number_density"])
    mole_size = int(p["molecule"]["atoms"])
    T   = float(p["thermostat"]["temperature"])

    pos = make_fcc_defect(L, rho, Lb)
    N = len(pos)
    N -= N%mole_size

    
    print(f"atoms={N}")
    print(f"input number_density={rho}")
    print(f"redefined number_density={N/V}")

    print(f"temperature={T}")

    vel = np.random.rand(N*3).reshape(N,3)-0.5
    vel -= np.mean(vel,axis=0)
    scaling_momentum(T,vel)

with open("init.dump","w") as file:
    print("ITEM: TIMESTEP",file=file)
    print(f"{0}",file=file)
    print("ITEM: NUMBER OF ATOMS",file=file)
    print(f"{N}",file=file)
    print("ITEM: BOX BOUNDS xx yy zz",file=file)
    print(f"0.0 {Lb[0]}",file=file)
    print(f"0.0 {Lb[1]}",file=file)
    print(f"0.0 {Lb[2]}",file=file)
    print("ITEM: ATOMS id mol type x y z vx vy vz radius",file=file)
    for i in range(N):
        a = pos[i]
        v = vel[i]
        # if i <= N//2: 
        molei = i//mole_size
        s=0
        #     s = 0 if i%mole_size == 0 else 1
        # else:
        #     molei = molei + 1
        #     s = 2
        print(f"{i} {molei} {s} {a[0]:.6f} {a[1]:.6f} {a[2]:.6f} {v[0]:.6f} {v[1]:.6f} {v[2]:.6f} {0.5}",file=file)
print("write init.dump")

