import numpy as np
import yaml 
import sys 

def read_data_in_a_step(f):
  _ = f.readline()
  STEP: int = int(f.readline())
  _ = f.readline()
  NUM_PARTICLE: int = int(f.readline())
  _ = f.readline()
  _, LX = map(float,f.readline().split())
  _, LY = map(float,f.readline().split())
  _, LZ = map(float,f.readline().split())
  _ = f.readline()

  DENSITY = NUM_PARTICLE/(LX*LY*LZ)
  
  NUM_MOLE = NUM_PARTICLE//PER_ATOMS
  q = np.empty((NUM_MOLE,PER_ATOMS,3))
  for i in range(NUM_MOLE):
    for j in range(PER_ATOMS):
      _,_,_,q[i,j,0],q[i,j,1],q[i,j,2],_,_,_,_ = map(float,f.readline().split())
  return STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q

def adjust_vector(dq,LX,LY,LZ):
  if dq[0] > +0.5*LX: dq[0]-=LX 
  if dq[0] < -0.5*LX: dq[0]+=LX 
  if dq[1] > +0.5*LY: dq[1]-=LY 
  if dq[1] < -0.5*LY: dq[1]+=LY 
  if dq[2] > +0.5*LZ: dq[2]-=LZ 
  if dq[2] < -0.5*LZ: dq[2]+=LZ 
  return dq


with open(sys.argv[1], "r") as f:
  PER_ATOMS = yaml.safe_load(f)["molecule"]["atoms"]

with open("LJrun.dump","r") as f, open("unwrap_xyz.dump","w") as fw:

  for k in range(0): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    print(f"step={STEP} skipped")

  USE_K = 100
  for k in range(USE_K): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    print(f"step={STEP} calculate")

    print("ITEM: TIMESTEP",file=fw)
    print(f"{STEP}",file=fw)
    print("ITEM: NUMBER OF ATOMS",file=fw)
    print(f"{NUM_MOLE*PER_ATOMS}",file=fw)
    print("ITEM: BOX BOUNDS xx yy zz",file=fw)
    print(f"0.0 {LX}",file=fw)
    print(f"0.0 {LY}",file=fw)
    print(f"0.0 {LZ}",file=fw)
    print("ITEM: ATOMS id mol type x y z vx vy vz radius",file=fw)

    cnt=0
    for i in range(NUM_MOLE):
      Rn = np.empty((PER_ATOMS,3))
      Rn[0] = q[i,0]
      Re = np.zeros(3)
      for j in range(PER_ATOMS-1):
        dq = adjust_vector(q[i,j+1] - q[i,j],LX,LY,LZ)
        Rn[j+1] = Rn[j] + dq
        Re += dq
      Rc = np.mean(Rn,axis=0)
      while Rc[0] > LX: 
        Rc[0] -= LX
        Rn[:,0] -= LX 
      while Rc[0] < 0.0: 
        Rc[0] += LX
        Rn[:,0] += LX 
      while Rc[1] > LY: 
        Rc[1] -= LY
        Rn[:,1] -= LY 
      while Rc[1] < 0.0: 
        Rc[1] += LY
        Rn[:,1] += LY 
      while Rc[2] > LZ: 
        Rc[2] -= LZ
        Rn[:,2] -= LZ 
      while Rc[2] < 0.0: 
        Rc[2] += LZ
        Rn[:,2] += LZ 

      for j in range(PER_ATOMS):
        print(f"{cnt} {i} {j} {Rn[j,0]:.6f} {Rn[j,1]:.6f} {Rn[j,2]:.6f} {0.0:.6f} {0.0:.6f} {0.0:.6f} {0.5}",file=fw)
        cnt+=1
