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

  for k in range(80): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"# step={STEP} skipped")

  USE_K = 20
  Re2 = np.zeros(PER_ATOMS-1)
  for k in range(USE_K): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"# step={STEP} calculate")
    Re_sq = [ [] for i in range(PER_ATOMS-1) ]
    for i in range(NUM_MOLE):
      Re = np.zeros((PER_ATOMS-1,3))
      cnt = np.zeros(PER_ATOMS-1,dtype=np.int64)
      for j in range(PER_ATOMS-1):
        dq = adjust_vector(q[i,j+1] - q[i,j],LX,LY,LZ)
        Re[:] += dq
        cnt[:] += 1
        for bj in range(PER_ATOMS-1):
          if cnt[bj] == bj+1:
            Re_sq[bj].append(np.sum(Re[bj]**2))
            Re[bj] = 0.0
            cnt[bj] = 0
      
    Re2 += np.array([ np.mean(Re_sq_) for Re_sq_ in Re_sq ] )
    
  for j in range(PER_ATOMS-1):
    print(j+1,Re2[j]/USE_K,Re2[j]/USE_K/(j+1))
