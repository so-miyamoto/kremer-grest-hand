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
  while dq[0] > +0.5*LX: dq[0]-=LX 
  while dq[0] < -0.5*LX: dq[0]+=LX 
  while dq[1] > +0.5*LY: dq[1]-=LY 
  while dq[1] < -0.5*LY: dq[1]+=LY 
  while dq[2] > +0.5*LZ: dq[2]-=LZ 
  while dq[2] < -0.5*LZ: dq[2]+=LZ 
  return dq

def get_Xp(Rn):
  N = len(Rn)
  Xp = np.empty_like(Rn)
  for p in range(N):
    i = np.arange(N)
    Xp[p] = (np.cos((p*np.pi*i)/(N-1)).dot(Rn[i,:]))/N + (0.5/N * (Rn[0]+Rn[N-1]))
  return Xp

with open(sys.argv[1], "r") as f:
  PER_ATOMS = yaml.safe_load(f)["molecule"]["atoms"]

with open("LJrun.dump","r") as f, open("Xp.dump","w") as fw:
  for k in range(0): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"# step={STEP} skipped")


  USE_K = 100
  Xp = np.zeros((PER_ATOMS,3))
  for k in range(USE_K): 
    STEP,NUM_MOLE,LX,LY,LZ,DENSITY,q = read_data_in_a_step(f)
    # print(f"# step={STEP} calculate")
    print(f"writing Xp.dump STEP={STEP}")
    print("ITEM: TIMESTEP",file=fw)
    print("0",file=fw)
    print("ITEM: NUMBER OF ATOMS",file=fw)
    print(NUM_MOLE*PER_ATOMS,file=fw)
    print("ITEM: BOX BOUNDS xx yy zz",file=fw)
    print(f"0.0 {LX}",file=fw)
    print(f"0.0 {LY}",file=fw)
    print(f"0.0 {LZ}",file=fw)
    print("ITEM: ATOMS id mol x y z",file=fw)

    Xp_k = np.zeros((NUM_MOLE,PER_ATOMS,3))
    for i in range(NUM_MOLE):
      Rn = np.zeros_like(Xp_k[i])
      Rn[0] = q[i,0]
      for j in range(PER_ATOMS-1):
        dq = adjust_vector(q[i,j+1] - q[i,j],LX,LY,LZ)
        Rn[j+1] = Rn[j] + dq
      Xp_k[i] = get_Xp(Rn)

    if k == 0: Xp_pre = Xp_k
            
    for i in range(NUM_MOLE):    
      for j in range(PER_ATOMS):
        Xp_adjust = Xp_pre[i][j] + adjust_vector(Xp_k[i][j]-Xp_pre[i][j],LX,LY,LZ)
        print(i,j,*Xp_adjust,file=fw)

    Xp_pre = Xp_k
    print("done")
