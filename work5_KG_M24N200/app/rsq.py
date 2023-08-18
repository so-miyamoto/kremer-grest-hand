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

with open("LJrun.dump","r") as f:

  for k in range(0): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"step={STEP} skipped")

  USE_K = 100
  for k in range(USE_K): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"step={STEP} calculate")
    Re_sq = []
    Re_qu = []
    Rg_sq = []
    for i in range(NUM_MOLE):
      Rn = np.empty((PER_ATOMS,3))
      Rn[0] = q[i,0]
      Re = np.zeros(3)
      for j in range(PER_ATOMS-1):
        dq = adjust_vector(q[i,j+1] - q[i,j],LX,LY,LZ)
        Rn[j+1] = Rn[j] + dq
        Re += dq
      Rc = np.mean(Rn,axis=0)
      Rg_sq_ = np.sum(np.mean((Rn - Rc)**2,axis=0))
      Re_sq.append(np.sum(Re**2))
      Re_qu.append(np.sum(Re**2)**2)
      Rg_sq.append(Rg_sq_)
      
    Re2 = np.mean(Re_sq) 
    Re4 = np.mean(Re_qu) 
    Rg2 = np.mean(Rg_sq)
    print("{:d} {:.4f} {:.4f} {:.4f} {:.4f}".format(STEP, Re2, Rg2, Re2/Rg2, Re4/Re2**2))
