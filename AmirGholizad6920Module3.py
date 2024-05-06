############# Author: Amir Gholizad #############
###### Submitted in partial fulfillment of ######
############ the requirements for ###############
############ the course CMSC 6920 ###############
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Data loading and Pre-processing
str = 'cyl_data_clean.csv'
ny = 199
nx = 449
U = np.loadtxt(str,delimiter=",")
print("Shape of the original data is: {}".format(U.shape))

# Plotting 4 random snapshots
fig, axs = plt.subplots(2,2)
axs[0,0].pcolormesh(np.reshape(U[:,0],(nx, ny)).T, cmap='RdBu', vmin=-3,vmax=3)
axs[0,0].set_xlabel("X")
axs[0,0].set_title("Snapshot 1")


axs[0,1].pcolormesh(np.reshape(U[:,10],(nx, ny)).T, cmap='RdBu', vmin=-3,vmax=3)
axs[0,1].set_xlabel("X")
axs[0,1].set_title("Snapshot 10")

axs[1,0].pcolormesh(np.reshape(U[:,100],(nx, ny)).T, cmap='RdBu', vmin=-3,vmax=3)
axs[1,0].set_xlabel("X")
axs[1,0].set_title("Snapshot 100")

axs[1,1].pcolormesh(np.reshape(U[:,150],(nx, ny)).T, cmap='RdBu', vmin=-3,vmax=3)
axs[1,1].set_xlabel("X")
axs[1,1].set_title("Snapshot 150")

fig.tight_layout()
plt.show()


###############################################
################ POD Modes ####################
###############################################

# Defining POD function
def POD(U, m):
    # SVD of Original data
    Phi, Sig, PsiT = np.linalg.svd(U, full_matrices=0)
    # Truncating the modes
    rPhi = Phi[:,:m]
    rSig = Sig[:m]
    rPsiT = PsiT[:m,:]
    # Refurbished U
    rU = np.dot(rPhi * rSig, rPsiT)
    return rU, rPhi, rSig, rPsiT


# Getting the SVD results for Original data
U0, Phi0, Sig0, PsiT0 = POD(U, 151)

# Computing the Energy for Original data
Energy0 = Sig0[0]**2/np.sum(Sig0**2)

# Creating a list of different snapshot choices
mlist = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 150]


# Computing ERRORs (1 - E_i/E_0) for different snapshot choices
ERR = []
rEnergylist = []
for m in mlist:
    rU, rPhi, rSig, rPsiT = POD(U, m)
    rEnergy = rSig[0]**2/np.sum(rSig**2)
    rEnergylist.append(rEnergy)
    ERR.append(abs(1-(rEnergy/Energy0)))

# Normalizing the errors
ERR = ERR/np.max(ERR)
print("The choosen values for m are: {}".format(mlist))
print("And the ERRORs are: {}".format(ERR))

# Plotting Number of snapshots VS ERROR levels
plt.stem(mlist, ERR, label='Number of snapshots VS ERROR levels')
plt.xlabel('Truncation Number')
plt.ylabel('ERROR')
plt.show()



# Caculate the Energy and Cumulative Energy
TEnergy0 = Sig0**2/np.sum(Sig0**2)
CumEnergy = np.cumsum(TEnergy0)


# Plotting Energy vs Number of snapshots
plt.plot(TEnergy0)
plt.ylabel('Energy')
plt.xlabel('Modes')
plt.show()


# Plotting Cumulative Energy
plt.plot(CumEnergy)
plt.ylabel('Cumulative Energy')
plt.xlabel('Modes')
plt.show()

# Getting the Truncation Number using the Cumulative Energy
TruncationNumber = len([e for e in CumEnergy if e < 0.99])
print("The number of modes that hold the 99/100 of the energy is: {}".format(TruncationNumber))


# Getting the truncated version of U
rU, rPhi, rSig, rPsiT = POD(U, TruncationNumber)
plt.plot(U, rU)
plt.ylabel('rU')
plt.xlabel('U')
plt.show()
print("The Correlation between original and refurbished data is: \n {}".format(np.corrcoef(rU[:,0], U[:,0])))
###############################################
################## DMD Modes ##################
###############################################
# Creating the U2 = A*U1 equation
U1 = U[:,:-1]
U2 =U[:,1:]


# Getting POD modes for U1
rU1, rPhi1, rSig1, rPsiT1 = POD(U1, TruncationNumber)


# Solving equation Atilde = rPhi.T * U2 * rPsiT.T/rSig for Atilde
At = np.linalg.multi_dot([rPhi1.T, U2, rPsiT1.T/rSig1])
Omega, R = np.linalg.eig(At)
print("Shapes of Omega and R are: {}, {}".format(Omega.shape, R.shape))


# computing dynamic modes = PHI_dmd = eigenvectors of A
PHI = np.linalg.multi_dot([U2, rPsiT1.T/rSig1, R])
# computing projected dynamic modes (PHI_dmd)
pPHI = np.dot(rPhi1, R)


# Plotting the first 6 DMD modes (PHI)
fig, axs = plt.subplots(2,3)
axs[0,0].pcolormesh(np.reshape(PHI[:,0].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[0,0].set_xlabel("X")
axs[0,0].set_title("DMD mode 1")

axs[0,1].pcolormesh(np.reshape(PHI[:,1].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[0,1].set_xlabel("X")
axs[0,1].set_title("DMD mode 2")

axs[0,2].pcolormesh(np.reshape(PHI[:,2].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[0,2].set_xlabel("X")
axs[0,2].set_title("DMD mode 3")

axs[1,0].pcolormesh(np.reshape(PHI[:,3].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[1,0].set_xlabel("X")
axs[1,0].set_title("DMD mode 4")

axs[1,1].pcolormesh(np.reshape(PHI[:,4].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[1,1].set_xlabel("X")
axs[1,1].set_title("DMD mode 5")

axs[1,2].pcolormesh(np.reshape(PHI[:,5].real,(nx, ny)).T, cmap='RdBu', vmin=-0.002,vmax=0.002)
axs[1,2].set_xlabel("X")
axs[1,2].set_title("DMD mode 6")
plt.show()


# Comparing middle PHI and projected PHI mode
fig, axs = plt.subplots(2,1)
axs[0].pcolormesh(np.reshape(pPHI[:,0].real,(nx, ny)).T, cmap='RdBu', vmin=-0.005,vmax=0.005)
axs[0].set_title("Projected PHI mode 1")
axs[1].pcolormesh(np.reshape(PHI[:,0].real,(nx, ny)).T, cmap='RdBu', vmin=-0.005,vmax=0.005)
axs[1].set_title("PHI mode 1")
plt.show()


# Sorting DMD eigen values
idx = Omega.argsort()[::-1]
# strongest DMD modes
print("The 5 strongest DMD modes are: {}".format(idx[:5]))
print("And strongest projected DMD mode is: {}".format(pPHI[:, idx[0]]))

# Assessing convergence of the DMD modes
theta = np.linspace(0,2*np.pi, 101)
circ = np.array([np.cos(theta) ,np.sin(theta)])
plt.figure(figsize=(5,5))
plt.plot(circ[0,:],circ[1,:])
plt.plot(Omega.real, Omega.imag, '.')
plt.title("Convergence of the eigenvalues")
plt.show()



# Reconstructing the data
u0 = U[:,0]
b = np.linalg.lstsq(pPHI, u0, rcond=None)[0] # pPHI*b = u0
print("The shape of the vector b is: {}".format(b.shape))
dt = 0.125
Lambda = np.log(Omega)/dt
for s in [4,99]:
    t = s*dt
    V = np.dot(pPHI, b*np.exp(Lambda*t)).real
    print("The correlattion matrix for V and the {}th mode of U is: \n {}".format(s,np.corrcoef(V,U[:,s])))
    fig, axs = plt.subplots(1,2)
    axs[0].pcolormesh(np.reshape(V.real,(nx, ny)).T, cmap='RdBu', vmin=-2,vmax=2)
    axs[1].pcolormesh(np.reshape(U[:,s].real,(nx, ny)).T, cmap='RdBu', vmin=-2,vmax=2)
    plt.show()






############################# Bonus Task #############################
################## Animating POD and DMD modes ##################

### Animating the original and truncated data
i = 0
for F in [U, rU]:

    # Create the initial plot
    fig, axs = plt.subplots(1,1)
    mesh = axs.pcolormesh(np.reshape(F[:,0],(nx, ny)).T, cmap='RdBu', vmin=-3, vmax=3)
    
    def update(frame):
        mesh.set_array(np.reshape(F[:,frame],(nx, ny)).T)
        return mesh,
    
    # Creating and saving the animation
    ani = FuncAnimation(fig, update, frames=len(U), interval=50, blit=True)
    # if i == 0:
    #     ani.save('U_animation.mp4')
    #     i =+ 1
    # else:
    #     ani.save('rU_animation.mp4')
    # Show the animation
    plt.show()


### Animating the DMD modes
# Create the initial plot
fig, axs = plt.subplots(1,1)
mesh = axs.pcolormesh(np.reshape(np.real(PHI[:,0]), (nx, ny)).T, cmap='RdBu', vmin=-0.002, vmax=0.002)


def update(frame):
    mesh.set_array(np.reshape(np.real(PHI[:,frame]), (nx, ny)).T)
    return mesh,

# Creating and saving the animation
ani = FuncAnimation(fig, update, frames=PHI.shape[1], interval=50, blit=True)
# ani.save('dmd_animation.mp4')
plt.show()


### Animating the predicted V snapshots
Vlist = [] 

for s in range(1,100): 
    t = s * dt
    V = np.dot(pPHI, b * np.exp(Lambda * t)).real
    Vlist.append(V)

VMatrix = np.column_stack(Vlist)

# Creating the initial plot
fig, axs = plt.subplots(1,1)
mesh = axs.pcolormesh(np.reshape(VMatrix[:,0],(nx, ny)).T, cmap='RdBu', vmin=-3, vmax=3)
    
def update(frame):
    mesh.set_array(np.reshape(VMatrix[:,frame],(nx, ny)).T)
    return mesh,
    
# Creating and saving the animation
ani = FuncAnimation(fig, update, frames=len(VMatrix-1), interval=50, blit=True)
# ani.save('V_animation.mp4')    
plt.show()


