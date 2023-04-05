import pickle
import numpy as np
from timeit import default_timer as timer
np.random.seed(0)



size = 300 #no.of lattice points is size*size

#initial number of molecules on the cell surface
init_Bh = 1500 # horizontal bcr init_Bh
init_Bv = 1500 # vertical bcr init_Bv
init_Ag = 1000 
init_Lfa = 3000
init_Icam = 2000


#Numbers allotted to each molecule
Bh_h_id = 1 # horizontal bcr head Bh_h_id
Bh_b_id = 2 # horizontal bcr body Bh_b_id 
Bh_t_id = 3 # horizontal bcr tail Bh_t_id
Bv_h_id = 4 # vertical bcr head Bv_h_id
Bv_b_id = 5 # vertical bcr body Bv_h_id
Bv_t_id = 6 # vertical bcr tail Bv_a_id
Ag_id = 7 # Ag_id
Lfa_id = 8 # Lfa_id
Icam_id = 9 # Icam_id

Bh_h_hA_id = 10# head of hori bcr bound to ag at head
Bh_b_hA_id = 11
Bh_t_hA_id = 12

Bv_h_hA_id = 13
Bv_b_hA_id = 14
Bv_t_hA_id = 15

Bh_h_tA_id = 16# BhhA_h,
Bh_b_tA_id = 17
Bh_t_tA_id = 18

Bv_h_tA_id = 19
Bv_b_tA_id = 20
Bv_t_tA_id = 21# BvtA

Bh_h_A2_id = 22 # BhA2_h_id
Bh_b_A2_id = 23
Bh_t_A2_id = 24
Bv_h_A2_id = 25# BvA2_h_id
Bv_b_A2_id = 26
Bv_t_A2_id = 27
LI_id = 28 # LI


#lists of molecules indices
Bcr_h = []
Bcr_v = []
Ag = []
Lfa = []
Icam = []
Bh_h_A = [] # BhhA_h,
Bh_t_A = []  
Bv_h_A = []
Bv_t_A = [] # BvtA
Bh_A2 = [] # BhA2_h_id
Bv_A2 = [] # BvA2_h_id
LI = [] # LI


"""
Bcr_h_num = len(Bcr_h) 
len(Bcr_v) = len(Bcr_v) 
Ag_num = len(Ag) 
len(Lfa) = len(Lfa)
len(Icam) = len(Icam) 
len(Bh_h_A) = len(Bh_h_A) 
len(Bh_t_A) = len(Bh_t_A) 
len(Bv_h_A) = len(Bv_h_A)
len(Bv_t_A) = len(Bv_t_A)
len(Bh_A2) = len(Bh_A2) 
len(Bv_A2) = len(Bv_A2)
len(Bh_A2) = len(LI)
"""

#p_diff_free = 0.8
#p_diff_complex = 0.1

#p_asso = 0.8
#p_diss = 0.3



#just to make sure everything runs/happens
p_diff_free = 1
p_diff_complex = 1
p_asso = 1
p_diss = 1
p_on_LI = 1#eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
p_off_LI = 1
p_on_BA = 1
p_off_BA = 1


z0 = 40 * (10**-9)  #Least distance between molecules(confirm)
x0 = size//2  #x-coordinate of center of cell/bilayer
y0 = size//2  #x-coordinate of center of cell/bilayer
Rb = 6 * (10**-6)  #radius of B lymphocytes

startie = timer()
surf = np.zeros((2,size,size),dtype = int)#initialising matrix
#print(surf)


#creating horizontal bcr molecules
ink = 0
while ink < init_Bh:
	row = np.random.randint(0,size-1)
	col = np.random.randint(1,size-2)
	
	#checking if the cells are avalable fr bcr to occupy
	if surf[1][row][col] == 0 and surf[1][row][col - 1] == 0 and surf[1][row][col + 1] == 0:
		
		#changing value in the matrix
		surf[1][row][col] = Bh_b_id
		surf[1][row][col - 1] = Bh_h_id
		surf[1][row][col + 1] = Bh_t_id
		
		#adding indices to list hbcr
		Bcr_h.append((1,row,col))

		ink += 1



#creating vertical bcr molecules ---not checked-----
ink = 0
while ink < init_Bv:
	row = np.random.randint(1,size-2)
	col = np.random.randint(0,size-1)
	
	#checking if the cells are avalable fr bcr to occupy
	if surf[1][row][col] == 0 and surf[1][row + 1][col] == 0 and surf[1][row - 1][col] == 0:
		
		#changing value in the matrix
		surf[1][row][col] = Bv_b_id
		surf[1][row - 1][col] = Bv_h_id
		surf[1][row + 1][col] = Bv_t_id
		
		#adding indices to list vbcr
		Bcr_v.append((1,row,col))

		ink += 1

 	
	


#creating ag molecules
ink = 0
while ink < init_Ag:

	row = np.random.randint(0,size-1)
	col = np.random.randint(0,size-1)
	
	if surf[0][row][col] == 0:

		surf[0][row][col] = Ag_id
		Ag.append((0,row,col))
		ink += 1

#creating lfa molecules
ink = 0
while ink < init_Lfa:

	row = np.random.randint(0,size-1)
	col = np.random.randint(0,size-1)

	if surf[1][row][col] == 0:
		surf[1][row][col] = Lfa_id
		Lfa.append((1,row,col))
		ink += 1

#creating icam molecules
ink = 0
while ink < init_Icam:

	row = np.random.randint(0,size-1)
	col = np.random.randint(0,size-1)

	if surf[0][row][col] == 0:
		surf[0][row][col] = Icam_id
		Icam.append((0,row,col))
		ink += 1

stopie = timer()
#print(stopie - startie)

#print(surf)
#print(len(Bcr_h) , len(Bcr_v) , len(Ag) , len(Lfa) , len(Icam) , len(Bh_h_A) , len(Bh_t_A) , len(Bv_h_A) , len(Bv_t_A) , len(Bh_A2) , len(Bv_A2) , len(Bh_A2))
#print(np.count_nonzero(surf))
#print(np.nonzero(surf))


#the total number of molecules in the matrix
nom = len(Bcr_h) + len(Bcr_v) + len(Ag) + len(Lfa) + len(Icam) + len(Bh_h_A) + len(Bh_t_A) + len(Bv_h_A) + len(Bv_t_A) + len(Bh_A2) + len(Bv_A2) + len(Bh_A2) 
#nom = len(Bcr_h) + len(Bcr_v) +len(Ag) + len(Lfa) + len(Icam) + len(Bh_h_A) + len(Bh_t_A) + len(Bv_h_A) + len(Bv_t_A) + len(Bh_A2) + len(Bv_A2) + len(LI)
#print(nom)
"""
#pickling 


with open('surf_{}'.format(size),'wb') as fil:
	pickle.dump(surf,fil)


data = [Bcr_h , Bcr_v , Ag , Lfa , Icam , surf] #fill in

def save(name):
	with open(name,'wb') as fil:
		pickle.dump(data,fil)

def load(name):
	with open(name,'rb') as fil:
		data = pickle.load(fil)

	Bcr_h = data[0]
	Bcr_v = data[1]
	Ag = [2]
	Lfa = [3]
	Icam = [4]
	surf = data[5]
	#fill in

#uncomment to save/load

save('surf_{}'.format(size))
#load('surf_{}'.format(size))

"""



     

#diffusion Functiona        
                        
def LIA_diffusion(coordinate,mol_list):
        
        x,y,z = coordinate

        
        if np.random.random() < p_diff_free: #cehcking if diffusion will occur:

                #direction = np.random.random()
                direction = 0.20

                #diffusing upward  
                if direction < 0.25:

                        ##print("up")
                       
                        if y - 1 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z] == 0:
                                        
                                        surf[x][y-1][z] = surf[x][y][z]
                                        surf[x][y][z] = 0
                                        mol_list.remove((x,y,z))
                                        mol_list.append((x,y - 1,z))
                                
                                        ##print((x,y,z),(x,y - 1,z))


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z] == 0:
                                        
                                        surf[x][y+1][z] = surf[x][y][z]
                                        surf[x][y][z] = 0
                                        mol_list.remove((x,y,z))
                                        
                                        mol_list.append((x,y + 1,z))

                                        ##print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0:
                                        
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = 0
                                        mol_list.remove((x,y,z))
                                        mol_list.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0:
                                        
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = 0
                                        mol_list.remove((x,y,z))
                                        mol_list.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))





def Bcr_v_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_free: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")
                       
                        if y - 2 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 2][z] == 0:
                                        
                                        surf[x][y - 2][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = 0
                                        
                                        Bcr_v.remove((x,y,z))
                                        Bcr_v.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 2 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 2][z] == 0:
                                        
                                        surf[x][y + 2][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = 0
                
                                        Bcr_v.remove((x,y,z))
                                        Bcr_v.append((x,y + 1,z))

                                        #print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0 and surf[x][y-1][z + 1] == 0 and surf[x][y+1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y + 1][z]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y - 1][z + 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0
                                        
                                        Bcr_v.remove((x,y,z))
                                        Bcr_v.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0 and surf[x][y - 1][z - 1] == 0 and surf[x][y +1][z -1] == 0:
                                         
                                        surf[x][y + 1][z - 1] = surf[x][y + 1][z]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        
                                        Bcr_v.remove((x,y,z))
                                        Bcr_v.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))


                
        

#print(Lfa[0])
#print(np.random.choice(np.array(Lfa)))
#zio = np.random.choice(Lfa)
"""
a1 = 1
a2 = size - 3
a3 = 2

b2 = a2 + 1
b3 = a3
c2 = a2 - 1
c3 = a3

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

Bcr_v.append((a1,a2,a3))

while (a1,1,b3) not in Bcr_v:
        Bcr_v_diffusion(Bcr_v[0],Bcr_v,1)
        print(surf)
        print(np.count_nonzero(surf))
-
"""                        

def Bcr_h_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_free: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")

                        if y - 1 >= 0:#checks if y is the topmost row


                       

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z - 1] == 0 and surf[x][y-1][z] == 0 and surf[x][y - 1][z + 1] == 0:
                                        
                                        surf[x][y - 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0
                                        
                                        Bcr_h.remove((x,y,z))
                                        Bcr_h.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))
                       
                       

                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row


                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z - 1] == 0 and surf[x][y + 1][z] == 0 and surf[x][y + 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y + 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0
                                        
                                        Bcr_h.remove((x,y,z))
                                        Bcr_h.append((x,y + 1,z))

                                        #print((x,y,z),(x,y + 1,z))

                                
                #diffusing to the right
                elif direction < 0.75:

                        ##print("right")
                       
                        if z + 2 < size:#checks if it is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 2] == 0:
                                        
                                        surf[x][y][z + 2] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                
                                        Bcr_h.remove((x,y,z))
                                        Bcr_h.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z + 1))


                

                           
                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 2 >= 0:#checks if y is the left-most column
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 2] == 0:
                                        
                                        surf[x][y][z - 2] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = 0
                
                                        Bcr_h.remove((x,y,z))
                                        Bcr_h.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z - 1))

"""                                
a1 = 1
a2 = size - 2
a3 = 3

b2 = a2 
b3 = a3 + 1
c2 = a2
c3 = a3 - 1

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

Bcr_h.append((a1,a2,a3))
print(surf)

while (a1,1,b3) not in Bcr_h:
        Bcr_h_diffusion(Bcr_h[0],Bcr_h,1)
        print(surf)
        print(np.count_nonzero(surf))
"""

                        
def LI_diffusion(coordinate):
        
        x,y,z = coordinate

        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20

                #diffusing upward  
                if direction < 0.25:

                        #print("up")
                       
                        if y - 1 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z] == 0 and surf[0][y - 1][z] == 0:
                                        
                                        surf[x][y-1][z] = surf[x][y][z]
                                        surf[x][y][z] = 0

                                        surf[0][y-1][z] = surf[0][y][z]
                                        surf[0][y][z] = 0
                        
                                        
                                        LI.remove((x,y,z))
                                        LI.append((x,y - 1,z))
                                
                                        #print((x,y,z),(x,y - 1,z))

                                        


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z] == 0 and surf[0][y + 1][z] == 0:
                                        
                                        surf[x][y+1][z] = surf[x][y][z]
                                        surf[x][y][z] = 0

                                        surf[0][y+1][z] = surf[0][y][z]
                                        surf[0][y][z] = 0
                                        
                                        LI.remove((x,y,z))
                                        LI.append((x,y + 1,z))

                                        #print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0 and surf[0][y][z + 1] == 0:
                                        
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = 0

                                        surf[0][y][z + 1] = surf[0][y][z]
                                        surf[0][y][z] = 0
                                        
                                        LI.remove((x,y,z))
                                        LI.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0 and surf[0][y][z - 1] == 0:
                                        
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = 0

                                        surf[0][y][z - 1] = surf[0][y][z]
                                        surf[0][y][z] = 0
                                        
                                        LI.remove((x,y,z))
                                        LI.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))

"""

a1,a2,a3 = 1,size - 1, size - 1

surf[a1][a2][a3] = 1000
surf[0][a2][a3] = 1000
LI.append((a1,a2,a3))
while (a1,1,a3) not in LI:
        LI_diffusion(LI[0],LI,1)
        print(surf)
        print(np.count_nonzero(surf))

"""
"""
a1 = 1
a2 = size - 2
a3 = 3

b2 = a2 
b3 = a3 + 1
c2 = a2
c3 = a3 - 1

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

surf[1 - a1][a2][a3] = 1000
surf[1 - a1][b2][b3] = 2000
surf[1 - a1][c2][c3] = 3000

Bcr_h.append((a1,a2,a3))
print(surf)

while (a1,1,b3) not in Bcr_h:
        Bcr_h_diffusion(Bcr_h[0],Bcr_h,1)
        print(surf)
        print(np.count_nonzero(surf))
"""


                        
def Bv_A2_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")
                       
                        if y - 2 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 2][z] == 0 and surf[0][y - 2][z] == 0 and surf[0][y][z] == 0 :
                                        
                                        surf[x][y - 2][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = 0

                                        surf[0][y][z] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                        surf[0][y - 2][z] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        
                                        
                                        Bv_A2.remove((x,y,z))
                                        Bv_A2.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 2 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 2][z] == 0 and surf[0][y][z] == 0 and surf[0][y + 2][z] == 0 :
                                        
                                        surf[x][y + 2][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = 0

                                        surf[0][y][z] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        surf[0][y + 2][z] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                        
                
                                        Bv_A2.remove((x,y,z))
                                        Bv_A2.append((x,y + 1,z))

                                        #print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0 and surf[x][y-1][z + 1] == 0 and surf[x][y+1][z + 1] == 0 and surf[0][y - 1][z + 1] == 0 and surf[0][y + 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y + 1][z]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y - 1][z + 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y - 1][z + 1] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        surf[0][y + 1][z + 1] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                        
                                        Bv_A2.remove((x,y,z))
                                        Bv_A2.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0 and surf[x][y - 1][z - 1] == 0 and surf[x][y +1][z -1] == 0 and surf[0][y - 1][z - 1] == 0 and surf[0][y + 1][z - 1] == 0:
                                         
                                        surf[x][y + 1][z - 1] = surf[x][y + 1][z]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y - 1][z - 1] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        surf[0][y + 1][z - 1] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0

                                        
                                        Bv_A2.remove((x,y,z))
                                        Bv_A2.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))

"""                         
a1 = 1
a2 = size - 3
a3 = 2

b2 = a2 + 1
b3 = a3
c2 = a2 - 1
c3 = a3

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

surf[0][b2][a3] = 1000
surf[0][c2][c3] = 3000

Bv_A2.append((a1,a2,a3))
print(surf)
while (a1,1,b3) not in Bv_A2:
        Bv_A2_diffusion(Bv_A2[0],Bv_A2,1)
        print(surf)
        print(np.count_nonzero(surf))        

"""

	
def Bh_A2_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")

                        if y - 1 >= 0:#checks if y is the topmost row


                       

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z - 1] == 0 and surf[x][y-1][z] == 0 and surf[x][y - 1][z + 1] == 0 and surf[0][y - 1][z - 1] == 0 and surf[0][y - 1][z + 1] == 0:
                                        
                                        surf[x][y - 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y - 1][z + 1] = surf[0][y][z + 1]
                                        surf[0][y - 1][z - 1] = surf[0][y][z - 1]
                                        surf[0][y][z + 1] = 0
                                        surf[0][y][z - 1] = 0                                        
                                        
                                        Bh_A2.remove((x,y,z))
                                        Bh_A2.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))
                       
                       

                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row


                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z - 1] == 0 and surf[x][y + 1][z] == 0 and surf[x][y + 1][z + 1] == 0 and surf[0][y + 1][z - 1] == 0 and surf[0][y + 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y + 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y + 1][z + 1] = surf[0][y][z + 1]
                                        surf[0][y + 1][z - 1] = surf[0][y][z - 1]
                                        surf[0][y][z + 1] = 0
                                        surf[0][y][z - 1] = 0 
                                        
                                        Bh_A2.remove((x,y,z))
                                        Bh_A2.append((x,y + 1,z))

                                        #print((x,y,z),(x,y + 1,z))

                                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 2 < size:#checks if it is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 2] == 0 and surf[0][y][z + 2] == 0 and surf[0][y][z] == 0:
                                        
                                        surf[x][y][z + 2] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0

                                        surf[0][y][z + 2] = surf[0][y][z + 1]
                                        surf[0][y][z] = surf[0][y][z - 1]
                                        surf[0][y][z + 1] = 0
                                        surf[0][y][z - 1] = 0
                                        
                
                                        Bh_A2.remove((x,y,z))
                                        Bh_A2.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z + 1))


                

                           
                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 2 >= 0:#checks if y is the left-most column
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 2] == 0 and surf[0][y][z - 2] == 0 and surf[0][y][z] == 0:
                                        
                                        surf[x][y][z - 2] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = 0

                                        surf[0][y][z - 2] = surf[0][y][z - 1]
                                        surf[0][y][z] = surf[0][y][z + 1]
                                        surf[0][y][z - 1] = 0
                                        surf[0][y][z + 1] = 0
                
                                        Bh_A2.remove((x,y,z))
                                        Bh_A2.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z - 1))



"""                               
a1 = 1
a2 = size - 2
a3 = 3

b2 = a2 
b3 = a3 + 1
c2 = a2
c3 = a3 - 1

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

surf[0][b2][b3] = 1000
surf[0][c2][c3] = 3000

Bcr_h.append((a1,a2,a3))
print(surf)

while (a1,1,b3) not in Bcr_h:
        Bcr_h_diffusion(Bcr_h[0],Bcr_h,1)
        print(surf)
        print(np.count_nonzero(surf))

"""


def Bv_h_A_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")
                       
                        if y - 2 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 2][z] == 0 and surf[0][y - 2][z] == 0:
                                        
                                        surf[x][y - 2][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = 0

                                        surf[0][y - 2][z] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        
                                        
                                        Bv_h_A.remove((x,y,z))
                                        Bv_h_A.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 2 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 2][z] == 0 and surf[0][y][z] == 0:
                                        
                                        surf[x][y + 2][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = 0

                                        surf[0][y][z] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        
                
                                        Bv_h_A.remove((x,y,z))
                                        Bv_h_A.append((x,y + 1,z))

                                        #print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0 and surf[x][y-1][z + 1] == 0 and surf[x][y+1][z + 1] == 0 and surf[0][y - 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y + 1][z]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y - 1][z + 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y - 1][z + 1] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                                                                
                                        Bv_h_A.remove((x,y,z))
                                        Bv_h_A.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0 and surf[x][y - 1][z - 1] == 0 and surf[x][y +1][z -1] == 0 and surf[0][y - 1][z - 1] == 0:
                                        
                                        surf[x][y + 1][z - 1] = surf[x][y + 1][z]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y - 1][z - 1] = surf[0][y - 1][z]
                                        surf[0][y - 1][z] = 0
                                        
                                        
                                        Bv_h_A.remove((x,y,z))
                                        Bv_h_A.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))
"""
a1 = 1
a2 = size - 3
a3 = 2

b2 = a2 + 1
b3 = a3
c2 = a2 - 1
c3 = a3

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

#surf[0][b2][a3] = 1000
surf[0][c2][c3] = 3000

Bv_A2.append((a1,a2,a3))
print(surf)
while (a1,1,b3) not in Bv_A2:
        Bv_h_A_diffusion(Bv_A2[0],Bv_A2,1)
        print(surf)
        print(np.count_nonzero(surf))   
"""


def Bv_t_A_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")
                       
                        if y - 2 >= 0:#checks if y is the topmost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 2][z] == 0 and surf[0][y][z] == 0 :
                                        
                                        surf[x][y - 2][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = 0

                                        surf[0][y][z] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                                                                
                                        
                                        Bv_t_A.remove((x,y,z))
                                        Bv_t_A.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))


                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 2 < size:#checks if y is the bottommost row

                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 2][z] == 0 and surf[0][y + 2][z] == 0 :
                                        
                                        surf[x][y + 2][z] = surf[x][y + 1][z]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y - 1][z]
                                        surf[x][y - 1][z] = 0

                                        surf[0][y + 2][z] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                        
                
                                        Bv_t_A.remove((x,y,z))
                                        Bv_t_A.append((x,y + 1,z))

                                        #print((x,y,z),(x,y+1,z))


                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 1 < size:#checks if y is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 1] == 0 and surf[x][y-1][z + 1] == 0 and surf[x][y+1][z + 1] == 0 and surf[0][y + 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y + 1][z]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y - 1][z + 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y + 1][z + 1] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0
                                        
                                        Bv_t_A.remove((x,y,z))
                                        Bv_t_A.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z+1))


                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 1 >= 0:#checks if y is the left-most column
                        
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 1] == 0 and surf[x][y - 1][z - 1] == 0 and surf[x][y +1][z -1] == 0 and surf[0][y + 1][z - 1] == 0:
                                         
                                        surf[x][y + 1][z - 1] = surf[x][y + 1][z]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y - 1][z]
                                        surf[x][y + 1][z] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y - 1][z] = 0

                                        surf[0][y + 1][z - 1] = surf[0][y + 1][z]
                                        surf[0][y + 1][z] = 0

                                        
                                        Bv_t_A.remove((x,y,z))
                                        Bv_t_A.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z-1))

"""                       
a1 = 1
a2 = size - 3
a3 = 2

b2 = a2 + 1
b3 = a3
c2 = a2 - 1
c3 = a3

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

surf[0][b2][a3] = 1000
#surf[0][c2][c3] = 3000

Bv_A2.append((a1,a2,a3))
print(surf)
while (a1,1,b3) not in Bv_A2:
        Bv_A2_diffusion(Bv_A2[0],Bv_A2,1)
        print(surf)
        print(np.count_nonzero(surf))        
"""

def Bh_h_A_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")

                        if y - 1 >= 0:#checks if y is the topmost row


                       

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z - 1] == 0 and surf[x][y-1][z] == 0 and surf[x][y - 1][z + 1] == 0 and surf[0][y - 1][z - 1] == 0:
                                        
                                        surf[x][y - 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y - 1][z - 1] = surf[0][y][z - 1]
                                        surf[0][y][z - 1] = 0                                        
                                        
                                        Bh_h_A.remove((x,y,z))
                                        Bh_h_A.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))
                       
                       

                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row


                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z - 1] == 0 and surf[x][y + 1][z] == 0 and surf[x][y + 1][z + 1] == 0 and surf[0][y + 1][z - 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y + 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y + 1][z - 1] = surf[0][y][z - 1]
                                        surf[0][y][z - 1] = 0 
                                        
                                        Bh_h_A.remove((x,y,z))
                                        Bh_h_A.append((x,y + 1,z))

                                        #print((x,y,z),(x,y + 1,z))

                                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 2 < size:#checks if it is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 2] == 0 and surf[0][y][z] == 0:
                                        
                                        surf[x][y][z + 2] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0

                                        surf[0][y][z] = surf[0][y][z - 1]
                                        surf[0][y][z - 1] = 0
                                        
                
                                        Bh_h_A.remove((x,y,z))
                                        Bh_h_A.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z + 1))


                

                           
                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 2 >= 0:#checks if y is the left-most column
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 2] == 0 and surf[0][y][z - 2] == 0:
                                        
                                        surf[x][y][z - 2] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = 0

                                        surf[0][y][z - 2] = surf[0][y][z - 1]
                                        surf[0][y][z - 1] = 0
                
                                        Bh_h_A.remove((x,y,z))
                                        Bh_h_A.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z - 1))



"""                               
a1 = 1
a2 = size - 2
a3 = 3

b2 = a2 
b3 = a3 + 1
c2 = a2
c3 = a3 - 1

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

#surf[0][b2][b3] = 1000
surf[0][c2][c3] = 3000

Bcr_h.append((a1,a2,a3))
print(surf)

while (a1,1,b3) not in Bcr_h:
        Bh_h_A_diffusion(Bcr_h[0],Bcr_h,1)
        print(surf)
        print(np.count_nonzero(surf))
"""

def Bh_t_A_diffusion(coordinate):
        
        x,y,z = coordinate
        
        if np.random.random() < p_diff_complex: #cehcking if diffusion will occur:

                direction = np.random.random()
                #direction = 0.20
                
                #diffusing upward  
                if direction < 0.25:

                        #print("up")

                        if y - 1 >= 0:#checks if y is the topmost row


                       

                                #check if free space available on the molecule's side of surface
                                if surf[x][y - 1][z - 1] == 0 and surf[x][y-1][z] == 0 and surf[x][y - 1][z + 1] == 0 and surf[0][y - 1][z + 1] == 0:
                                        
                                        surf[x][y - 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y - 1][z] = surf[x][y][z]
                                        surf[x][y - 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y - 1][z + 1] = surf[0][y][z + 1]
                                        surf[0][y][z + 1] = 0                                     
                                        
                                        Bh_t_A.remove((x,y,z))
                                        Bh_t_A.append((x,y - 1,z))

                                        #print((x,y,z),(x,y - 1,z))
                       
                       

                #diffusing downward 
                elif direction < 0.5:

                        #print("down")
                       
                        if y + 1 < size:#checks if y is the bottommost row


                                #check if free space available on the molecule's side of surface
                                if surf[x][y + 1][z - 1] == 0 and surf[x][y + 1][z] == 0 and surf[x][y + 1][z + 1] == 0 and surf[0][y + 1][z + 1] == 0:
                                        
                                        surf[x][y + 1][z + 1] = surf[x][y][z + 1]
                                        surf[x][y + 1][z] = surf[x][y][z]
                                        surf[x][y + 1][z - 1] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0
                                        surf[x][y][z] = 0
                                        surf[x][y][z + 1] = 0

                                        surf[0][y + 1][z + 1] = surf[0][y][z + 1]
                                        surf[0][y][z + 1] = 0
                                        
                                        Bh_t_A.remove((x,y,z))
                                        Bh_t_A.append((x,y + 1,z))

                                        #print((x,y,z),(x,y + 1,z))

                                
                #diffusing to the right
                elif direction < 0.75:

                        #print("right")
                       
                        if z + 2 < size:#checks if it is the right-most column

                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z + 2] == 0 and surf[0][y][z + 2] == 0:
                                        
                                        surf[x][y][z + 2] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = 0

                                        surf[0][y][z + 2] = surf[0][y][z + 1]
                                        surf[0][y][z + 1] = 0
                                        
                
                                        Bh_t_A.remove((x,y,z))
                                        Bh_t_A.append((x,y,z + 1))

                                        #print((x,y,z),(x,y,z + 1))


                

                           
                #diffusing to the left
                else:

                        #print("left")
                       
                        if z - 2 >= 0:#checks if y is the left-most column
                                #check if free space available on the molecule's side of surface
                                if surf[x][y][z - 2] == 0 and surf[0][y][z] == 0:
                                        
                                        surf[x][y][z - 2] = surf[x][y][z - 1]
                                        surf[x][y][z - 1] = surf[x][y][z]
                                        surf[x][y][z] = surf[x][y][z + 1]
                                        surf[x][y][z + 1] = 0

                                        surf[0][y][z] = surf[0][y][z + 1]
                                        surf[0][y][z + 1] = 0
                
                                        Bh_t_A.remove((x,y,z))
                                        Bh_t_A.append((x,y,z - 1))

                                        #print((x,y,z),(x,y,z - 1))



"""                             
a1 = 1
a2 = size - 2
a3 = 3

b2 = a2 
b3 = a3 + 1
c2 = a2
c3 = a3 - 1

surf[a1][a2][a3] = 1000
surf[a1][b2][b3] = 2000
surf[a1][c2][c3] = 3000

surf[0][b2][b3] = 1000
#surf[0][c2][c3] = 3000

Bcr_h.append((a1,a2,a3))
print(surf)

while (a1,1,b3) not in Bcr_h:
        Bh_t_A_diffusion(Bcr_h[0],Bcr_h,1)
        print(surf)
        print(np.count_nonzero(surf))
"""


"""
remove the parts btween #... and #... in all reaction functions
"""



#reaction Functions------



def LI_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_off_LI:

                LI.remove((x,y,z))
                Lfa.append((x,y,z))
                Icam.append((0,y,z))

                #....
                #print(x,y,z)
                surf[x][y][z] = Lfa_id
                surf[0][y][z] = Icam_id
                #...
"""
a1,a2,a3 = 0,3,4
b1,b2,b3 = 1,a2,a3

surf[a1][a2][a3] = 2000
surf[b1][b2][b3] = 2000

LI.append((a1,a2,a3))
Lfa_id = 1000
Icam_id = 1000

print(surf)
LI_reaction(LI[0],1)
print(surf)
print(Lfa,Icam,LI)
"""

def Lfa_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_on_LI:

                if surf[0][y][z] == Icam_id:

                        LI.append((1,y,z))
                        Lfa.remove((1,y,z))
                        Icam.remove((0,y,z))

                        #...
                        surf[1][y][z] = LI_id
                        surf[0][y][z] = LI_id
                        #...

"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = 0,a2,a3

surf[a1][a2][a3] = 1000
surf[b1][b2][b3] = 1000

Lfa.append((a1,a2,a3))
Icam.append((b1,b2,b3))
Lfa_id = 1000
Icam_id = 1000
LI_id = 2000

print(surf)
L_reaction(Lfa[0],1)
print(surf)
print(Lfa,Icam,LI)
"""

def Icam_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_on_LI:

                if surf[0][y][z] == Lfa_id:

                        LI.append((1,y,z))
                        Lfa.remove((1,y,z))
                        Icam.remove((0,y,z))

                        #...
                        surf[1][y][z] = LI_id
                        surf[0][y][z] = LI_id
                        #...

"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = 0,a2,a3

surf[a1][a2][a3] = 1000
surf[b1][b2][b3] = 1000

Lfa.append((a1,a2,a3))
Icam.append((b1,b2,b3))
Lfa_id = 1000
Icam_id = 1000
LI_id = 2000

print(surf)
L_reaction(Icam[0],1)
print(Lfa,Icam,LI)
print(surf)
"""

def Bcr_h_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_on_BA:

                #if random.random() < 0.5:#head is chosen
                if 1>2:

                        if surf[0][y][z - 1] == Ag_id:

                                Bh_h_A.append((x,y,z))#Bh_h_A
                                Bcr_h.remove((x,y,z))
                                Ag.remove((0,y,z-1))

                                #...
                                surf[x][y][z - 1] = Bh_h_hA_id
                                surf[x][y][z] = Bh_b_hA_id
                                surf[x][y][z + 1] = Bh_t_hA_id
                                surf[0][y][z - 1] = Bh_h_hA_id
                                #...

                else:

                        if surf[0][y][z + 1] == Ag_id:#tail is chosen

                                Bh_t_A.append((x,y,z))
                                Bcr_h.remove((x,y,z))
                                Ag.remove((0,y,z + 1))

                                #...
                                surf[x][y][z - 1] = Bh_h_tA_id
                                surf[x][y][z] = Bh_b_tA_id
                                surf[x][y][z + 1] = Bh_t_tA_id
                                surf[0][y][z + 1] = Bh_t_tA_id
                                #...
"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
bh_h_A_id = 10000
bh_t_A_id = 1000000

surf[a1][a2][a3] = 300
surf[b1][b2][b3] = 200
surf[c1][c2][c3] = 400
surf[0][b2][b3] = surf[0][c2][c3] = Ag_id

Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bcr_h.append((a1,a2,a3))

print(surf)

Bcr_h_reaction(Bcr_h[0],1)

print(Bcr_h,Ag,Bh_h_A,Bh_t_A)                        

print(surf)
                        
"""
                        
def Bcr_v_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_on_BA:

                if np.random.random() < 0.5:#head is chosen
                #if 1<2:

                        if surf[0][y - 1][z] == Ag_id:

                                Bv_h_A.append((x,y,z))#Bh_h_A
                                Bcr_v.remove((x,y,z))
                                Ag.remove((0,y - 1,z))

                                #...
                                surf[x][y - 1][z] = Bv_h_hA_id
                                surf[x][y][z] = Bv_b_hA_id
                                surf[x][y + 1][z] = Bv_t_hA_id
                                surf[0][y - 1][z] = Bv_h_hA_id
                                #...

                else:

                        if surf[0][y + 1][z] == Ag_id:#tail is chosen

                                Bv_t_A.append((x,y,z))
                                Bcr_v.remove((x,y,z))
                                Ag.remove((0,y + 1,z))

                                #...
                                surf[x][y - 1][z] = Bv_h_tA_id
                                surf[x][y][z] = Bv_b_tA_id
                                surf[x][y + 1][z] = Bv_t_tA_id
                                surf[0][y + 1][z] = Bv_t_tA_id
                                #...                 


"""        
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
bv_h_A_id = 10000
bv_t_A_id = 100000

surf[a1][a2][a3] = 300
surf[b1][b2][b3] = 200
surf[c1][c2][c3] = 400
surf[0][b2][b3] = surf[0][c2][c3] = Ag_id

Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bcr_v.append((a1,a2,a3))

print(surf)

Bcr_v_reaction(Bcr_v[0],1)

print(Bcr_v,Ag,Bv_h_A,Bv_t_A)                        

print(surf)
"""        

def Bh_A2_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_off_BA:

                if np.random.random() < 0.5:#head is chosen

                        Bh_A2.remove((x,y,z))
                        Bh_t_A.append((x,y,z))
                        Ag.append((0,y,z-1))

                        #...
                        surf[x][y][z - 1] = Bh_h_tA_id
                        surf[x][y][z] = Bh_b_tA_id
                        surf[x][y][z + 1] = Bh_t_tA_id
                        surf[0][y][z + 1] = Bh_t_tA_id
                        surf[0][y][z - 1] = Ag_id
                        
                        #...


                else:

                        Bh_A2.remove((x,y,z))
                        Bh_h_A.append((x,y,z))
                        Ag.append((0,y,z + 1))

                        #...
                        surf[x][y][z - 1] = Bh_h_hA_id
                        surf[x][y][z] = Bh_b_hA_id
                        surf[x][y][z + 1] = Bh_t_hA_id
                        surf[0][y][z - 1] = Bh_h_hA_id
                        surf[0][y][z + 1] = Ag_id
                        #...


"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
bh_h_A_id = 10000
bh_t_A_id = 1000000
Bh_A2_id = 3000

surf[a1][a2][a3] = Bh_A2_id
surf[b1][b2][b3] = Bh_A2_id
surf[c1][c2][c3] = Bh_A2_id
surf[0][b2][b3] = surf[0][c2][c3] = Bh_A2_id

#Ag.append((0,b2,b3))
#Ag.append((0,c2,c3))
Bh_A2.append((a1,a2,a3))

print(surf)

Bh_A2_reaction(Bh_A2[0],1)

print(Bh_A2,Ag,Bh_h_A,Bh_t_A)                        

print(surf)                        
"""


def Bv_A2_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_off_BA:

                if np.random.random() < 0.5:#head is chosen

                        Bv_A2.remove((x,y,z))#Bh_h_A
                        Bv_t_A.append((x,y,z))
                        Ag.append((0,y - 1,z))

                        #...
                        surf[x][y - 1][z] = Bv_h_tA_id
                        surf[x][y][z] = Bv_b_tA_id
                        surf[x][y + 1][z] = Bv_t_tA_id
                        surf[0][y + 1][z] = Bv_t_tA_id
                        surf[0][y - 1][z] = Ag_id
                        
                        #...


                else:

                        Bv_A2.remove((x,y,z))
                        Bv_h_A.append((x,y,z))
                        Ag.append((0,y + 1,z))

                        #...
                        surf[x][y - 1][z] = Bv_h_hA_id
                        surf[x][y][z] = Bv_b_hA_id
                        surf[x][y + 1][z] = Bv_t_hA_id
                        surf[0][y - 1][z] = Bv_h_hA_id
                        surf[0][y + 1][z] = Ag_id
                        #...


"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
bv_h_A_id = 10000
bv_t_A_id = 1000000
Bv_A2_id = 3000

surf[a1][a2][a3] = Bv_A2_id
surf[b1][b2][b3] = Bv_A2_id
surf[c1][c2][c3] = Bv_A2_id
surf[0][b2][b3] = surf[0][c2][c3] = Bv_A2_id

#Ag.append((0,b2,b3))
#Ag.append((0,c2,c3))
Bv_A2.append((a1,a2,a3))

print(surf)

Bv_A2_reaction(Bv_A2[0],1)

print(Bv_A2,Ag,Bv_h_A,Bv_t_A)                        

print(surf)                        
"""

def Bh_h_A_reaction(coordinate):

        x,y,z = coordinate
        
        if np.random.random() < 0.5:#head is chosen, and therefore dissociation occurs

                if np.random.random() < p_off_BA:

                        Bh_h_A.remove((x,y,z))#Bh_h_A
                        Bcr_h.append((x,y,z))
                        Ag.append((0,y,z-1))

                        #...
                        surf[x][y][z - 1] = Bh_h_id
                        surf[x][y][z] = Bh_b_id
                        surf[x][y][z + 1] = Bh_t_id
                        #surf[0][y][z + 1] = Bh_h_id
                        surf[0][y][z - 1] = Ag_id
                        #...

        else: #tail is chosen and therefore reaction occurs

                if np.random.random() < p_on_BA:

                        if surf[0][y][z + 1] == Ag_id:#tail is chosen
                      

                                Bh_A2.append((x,y,z))
                                Bh_h_A.remove((x,y,z))
                                Ag.remove((0,y,z + 1))

                                        #...
                                surf[x][y][z - 1] = Bh_h_A2_id
                                surf[x][y][z] = Bh_b_A2_id
                                surf[x][y][z + 1] = Bh_t_A2_id
                                surf[0][y][z + 1] = Bh_t_A2_id
                                surf[0][y][z - 1] = Bh_h_A2_id
                                        #...


"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
Bh_h_id = 200
Bh_h_A_id = 10000
#bh_t_A_id = 20000
Bh_A2_id = 3000

surf[a1][a2][a3] = Bh_h_A_id
surf[b1][b2][b3] = Bh_h_A_id
surf[c1][c2][c3] = Bh_h_A_id
surf[0][b2][b3] = Bh_h_A_id
surf[0][c2][c3] = Ag_id
#Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bh_h_A.append((a1,a2,a3))

print(surf)

Bh_h_A_reaction(Bh_h_A[0],1,1)

print(Ag,Bh_h_A,Bh_A2,Bcr_h)                        

print(surf)                        
"""               
                        

def Bh_t_A_reaction(coordinate):

        x,y,z = coordinate
        
        if np.random.random() < 0.5:#head is chosen, and therefore dissociation occurs

                if np.random.random() < p_on_BA:

                        if surf[0][y][z - 1] == Ag_id:
                      

                                Bh_A2.append((x,y,z))
                                Bh_t_A.remove((x,y,z))
                                Ag.remove((0,y,z - 1))

                                        #...
                                surf[x][y][z - 1] = Bh_h_A2_id
                                surf[x][y][z] = Bh_b_A2_id
                                surf[x][y][z + 1] = Bh_t_A2_id
                                surf[0][y][z + 1] = Bh_t_A2_id
                                surf[0][y][z - 1] = Bh_h_A2_id
                                        #...


        else:

                if np.random.random() < p_off_BA:

                        Bh_t_A.remove((x,y,z))
                        Bcr_h.append((x,y,z))
                        Ag.append((0,y,z + 1))

                        #...
                        surf[x][y][z - 1] = Bh_h_id
                        surf[x][y][z] = Bh_b_id
                        surf[x][y][z + 1] = Bh_t_id
                        #surf[0][y][z + 1] = Bh_h_id
                        surf[0][y][z + 1] = Ag_id
                        #...
"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
Bh_h_id = 200
Bh_t_A_id = 10000
Bh_t_A_id = 20000
Bh_A2_id = 3000

surf[a1][a2][a3] = Bh_t_A_id
surf[b1][b2][b3] = Bh_t_A_id
surf[c1][c2][c3] = Bh_t_A_id
surf[0][b2][b3] = Ag_id
surf[0][c2][c3] = Bh_t_A_id
#Ag.append((0,b2,b3))
Ag.append((0,b2,b3))
Bh_t_A.append((a1,a2,a3))

print(surf)

Bh_t_A_reaction(Bh_t_A[0],1,1)

print(Ag,Bh_t_A,Bh_A2,Bcr_h)                        

print(surf)               
"""


def Bv_h_A_reaction(coordinate):

        x,y,z = coordinate
        
        if np.random.random() < 0.5:#head is chosen, and therefore dissociation occurs

                if np.random.random() < p_off_BA:

                        Bv_h_A.remove((x,y,z))#Bh_h_A
                        Bcr_v.append((x,y,z))
                        Ag.append((0,y-1,z))

                        #...
                        surf[x][y - 1][z] = Bv_h_id
                        surf[x][y][z] = Bv_b_id
                        surf[x][y + 1][z] = Bv_t_id
                        #surf[0][y][z + 1] = Bh_h_id
                        surf[0][y - 1][z] = Ag_id
                        #...

        else: #tail is chosen and therefore reaction occurs

                if np.random.random() < p_on_BA:

                        if surf[0][y + 1][z] == Ag_id:#tail is chosen
                      

                                Bv_A2.append((x,y,z))
                                Bv_h_A.remove((x,y,z))
                                Ag.remove((0,y + 1,z))

                                #...
                                surf[x][y - 1][z] = Bv_h_A2_id
                                surf[x][y][z] = Bv_b_A2_id
                                surf[x][y + 1][z] = Bv_t_A2_id
                                surf[0][y + 1][z] = Bv_t_A2_id
                                surf[0][y - 1][z] = Bv_h_A2_id
                                #...

"""        
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
Bcr_v_id = 200
Bv_h_A_id = 10000
#bh_t_A_id = 20000
Bv_A2_id = 3000

surf[a1][a2][a3] = Bv_h_A_id
surf[b1][b2][b3] = Bv_h_A_id
surf[c1][c2][c3] = Bv_h_A_id
surf[0][b2][b3] = Bv_h_A_id
surf[0][c2][c3] = Ag_id
#Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bv_h_A.append((a1,a2,a3))

print(surf)

Bv_h_A_reaction(Bv_h_A[0],1,1)

print(Ag,Bv_h_A,Bv_A2,Bcr_v)                        

print(surf)                   
"""


def Bv_t_A_reaction(coordinate):

        x,y,z = coordinate
        
        if np.random.random() < 0.5:#head is chosen, and therefore dissociation occurs

                if np.random.random() < p_on_BA:

                        if surf[0][y - 1][z] == Ag_id:
                      

                                Bv_A2.append((x,y,z))
                                Bv_t_A.remove((x,y,z))
                                Ag.remove((0,y - 1,z))

                                #...
                                surf[x][y - 1][z] = Bv_h_A2_id
                                surf[x][y][z] = Bv_b_A2_id
                                surf[x][y + 1][z] = Bv_t_A2_id
                                surf[0][y + 1][z] = Bv_t_A2_id
                                surf[0][y - 1][z] = Bv_h_A2_id
                                #...


        else:

                if np.random.random() < p_off_BA:

                        Bv_t_A.remove((x,y,z))
                        Bcr_v.append((x,y,z))
                        Ag.append((0,y + 1,z))

                        #...
                        surf[x][y - 1][z] = Bv_h_id
                        surf[x][y][z] = Bv_b_id
                        surf[x][y + 1][z] = Bv_t_id
                        #surf[0][y][z + 1] = Bh_h_id
                        surf[0][y + 1][z] = Ag_id
                        #...

"""
a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
Bcr_v_id = 200
Bv_t_A_id = 10000
Bv_t_A_id = 20000
Bv_A2_id = 3000

surf[a1][a2][a3] = Bv_t_A_id
surf[b1][b2][b3] = Bv_t_A_id
surf[c1][c2][c3] = Bv_t_A_id
surf[0][b2][b3] = Ag_id
surf[0][c2][c3] = Bv_t_A_id
#Ag.append((0,b2,b3))
Ag.append((0,b2,b3))
Bv_t_A.append((a1,a2,a3))

print(surf)

Bv_t_A_reaction(Bv_t_A[0],1,1)

print(Ag,Bv_t_A,Bv_A2,Bcr_v)                        

print(surf)               
"""

                
def Ag_reaction(coordinate):

        x,y,z = coordinate

        if np.random.random() < p_on_BA:

                if surf[1][y][z] != 0:

                        if surf[1][y][z] == Bh_h_id:

                                #print(surf[1][y][z])
                                
                                Bh_h_A.append((1,y,z+1))
                                Bcr_h.remove((1,y,z+1))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bh_h_hA_id
                                surf[1][y][z+1] = Bh_b_hA_id
                                surf[1][y][z+2] = Bh_t_hA_id
                                surf[0][y][z] = Bh_h_hA_id
                                

                        elif surf[1][y][z] == Bh_t_id:

                                Bh_t_A.append((1,y,z-1))
                                Bcr_h.remove((1,y,z-1))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bh_t_tA_id
                                surf[1][y][z-1] = Bh_b_tA_id
                                surf[1][y][z-2] = Bh_h_tA_id
                                surf[0][y][z] = Bh_t_tA_id

                                

                        elif surf[1][y][z] == Bv_h_id:

                                Bv_h_A.append((1,y+1,z))
                                Bcr_v.remove((1,y+1,z))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bv_h_hA_id
                                surf[1][y + 1][z] = Bv_b_hA_id
                                surf[1][y + 2][z] = Bv_t_hA_id
                                surf[0][y][z] = Bv_h_hA_id


                        elif surf[1][y][z] == Bv_t_id:

                                Bv_t_A.append((1,y - 1,z))
                                Bcr_v.remove((1,y - 1,z))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bv_t_tA_id
                                surf[1][y - 1][z] = Bv_b_tA_id
                                surf[1][y - 2][z] = Bv_h_tA_id
                                surf[0][y][z] = Bv_t_tA_id


                        elif surf[1][y][z] == Bh_t_hA_id:


                                Bh_A2.append((1,y,z-1))
                                Bh_h_A.remove((1,y,z-1))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bh_t_A2_id
                                surf[1][y][z-1] = Bh_b_A2_id
                                surf[1][y][z-2] = Bh_h_A2_id
                                surf[0][y][z] = Bh_t_A2_id
                                surf[0][y][z-2] = Bh_h_A2_id


                        elif surf[1][y][z] == Bh_h_tA_id:


                                Bh_A2.append((1,y,z + 1))
                                Bh_t_A.remove((1,y,z + 1))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bh_h_A2_id
                                surf[1][y][z + 1] = Bh_b_A2_id
                                surf[1][y][z + 2] = Bh_t_A2_id
                                surf[0][y][z] = Bh_h_A2_id
                                surf[0][y][z + 2] = Bh_t_A2_id



                        elif surf[1][y][z] == Bv_t_hA_id:


                                Bv_A2.append((1,y - 1,z))
                                Bv_h_A.remove((1,y - 1,z))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bv_t_A2_id
                                surf[1][y - 1][z] = Bv_b_A2_id
                                surf[1][y - 2][z] = Bv_h_A2_id
                                surf[0][y][z] = Bv_t_A2_id
                                surf[0][y - 2][z] = Bv_h_A2_id


                        elif surf[1][y][z] == Bv_h_tA_id:


                                Bv_A2.append((1,y + 1,z))
                                Bv_t_A.remove((1,y + 1,z))
                                Ag.remove((x,y,z))

                                surf[1][y][z] = Bv_h_A2_id
                                surf[1][y + 1][z] = Bv_b-A2_id
                                surf[1][y + 2][z] = Bv_t_A2_id
                                surf[0][y][z] = Bv_h_A2_id
                                surf[0][y + 2][z] = Bv_t_A2_id



"""
a1,a2,a3 = 1,3,3
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
Bcr_v_id = 200
Bv_t_A_id = 10000
Bv_t_A_id = 20000
Bv_A2_id = 3000

surf[a1][a2][a3] = Bv_t_A_id
surf[b1][b2][b3] = Bv_t_A_id
surf[c1][c2][c3] = Bv_t_A_id
surf[0][b2][b3] = Ag_id
surf[0][c2][c3] = Bv_t_A_id
#Ag.append((0,b2,b3))
Ag.append((0,b2,b3))
Bv_t_A.append((a1,a2,a3))

print(surf)

Ag_reaction(Ag[0],1)

print(Ag,Bv_t_A,Bv_A2,Bcr_v)                        

print(surf)           

a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
Bcr_v_id = 200
Bv_h_A_id = 10000
#bh_t_A_id = 20000
Bv_A2_id = 3000

surf[a1][a2][a3] = Bv_h_A_id
surf[b1][b2][b3] = Bv_h_A_id
surf[c1][c2][c3] = Bv_h_A_id
surf[0][b2][b3] = Bv_h_A_id
surf[0][c2][c3] = Ag_id
#Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bv_h_A.append((a1,a2,a3))

print(surf)
Ag_reaction(Ag[0],1)

print(Ag,Bv_h_A,Bv_A2,Bcr_v)                        

print(surf)                   

a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
Bh_h_id = 200
Bh_t_A_id = 10000
Bh_t_A_id = 20000
Bh_A2_id = 3000

surf[a1][a2][a3] = Bh_t_A_id
surf[b1][b2][b3] = Bh_t_A_id
surf[c1][c2][c3] = Bh_t_A_id
surf[0][b2][b3] = Ag_id
surf[0][c2][c3] = Bh_t_A_id
#Ag.append((0,b2,b3))
Ag.append((0,b2,b3))
Bh_t_A.append((a1,a2,a3))

print(surf)

Ag_reaction(Ag[0],1)

print(Ag,Bh_t_A,Bh_A2,Bcr_h)                        

print(surf)               


a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2,a3-1
c1,c2,c3 = a1,a2,a3+1

Ag_id = 100
Bh_h_id = 200
Bh_h_A_id = 10000
#bh_t_A_id = 20000
Bh_A2_id = 3000

surf[a1][a2][a3] = Bh_h_A_id
surf[b1][b2][b3] = Bh_h_A_id
surf[c1][c2][c3] = Bh_h_A_id
surf[0][b2][b3] = Bh_h_A_id
surf[0][c2][c3] = Ag_id
#Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bh_h_A.append((a1,a2,a3))

print(surf)

Ag_reaction(Ag[0],1)

print(Ag,Bh_h_A,Bh_A2,Bcr_h)                        

print(surf)                        


a1,a2,a3 = 1,3,4
b1,b2,b3 = a1,a2-1,a3
c1,c2,c3 = a1,a2+1,a3

Ag_id = 100
Bv_h_id = 1000
Bv_t_id = 2000
Bv_h_A_id = 222
Bv_t_A_id = 333
Bv_A2_id = 9999
surf[a1][a2][a3] = 1000
surf[b1][b2][b3] = Bv_h_id
surf[c1][c2][c3] = Bv_t_id
surf[0][b2][b3] = surf[0][c2][c3] = Ag_id

Ag.append((0,b2,b3))
Ag.append((0,c2,c3))
Bcr_v.append((a1,a2,a3))

print(surf)

Ag_reaction(Ag[1],1)
Ag_reaction(Ag[0],1)
print(Bcr_v,Ag,Bv_h_A,Bv_t_A)                        

print(surf)
"""


#Height Function
def height(x1,y1):

        return(z0 + Rb - (Rb**2 - ((x1 - x0)**2 + (y1 - y0)**2))**(1/2))


#calculating radial distance of each molecule
def radial_distance(x2,y2):

    return(((x0 - x2)**2 + (y0 - y2)**2)**(1/2))


#To find average radial distance of a molecule-type from the center
def synform(mol_list):

    if len(mol_list) != 0:

        Av_radial_distance = 0
    
        for mol in mol_list:

            Av_radial_distance += radial_distance(mol[1],mol[2])    

        return(Av_radial_distance/len(mol_list))

    else:

        return 0


synba = []
synli = []




def scanner():

#head=body=tail

    ag_s = 0
    lfa_s = 0
    icam_s = 0
    li_s = 0
    
    bcrh_s = (0,0,0)
    bcrv_s = (0,0,0)
    bahh_s = (0,0,0)
    bahv_s = (0,0,0)
    bath_s = (0,0,0)
    batv_s = (0,0,0)
    ba2h_s = (0,0,0)
    ba2v_s = (0,0,0)

    for i in range(size):
        for j in range(size):
            for k in range(2):

                if surf[k][i][j] == Bh_h_id:
                    heads += 1

                elif surf[k][i][j] == Bh_b_id:
                    bodies += 1

                elif surf[k][i][j] == Bh_t_id:
                    tails += 1

    return (heads,bodies,tails)




#Frameee

start = timer()        

"""
To check:

1. Order of increasing probability for cumuprob

"""

for time in range(100):

        print(time)
        #total number of molecules in the matrix
        #nom = Bcr_h_len + Bcr_v_len + Ag_len + Lfa_len + Icam_len + Bh_h_A_len + Bh_t_A_len + Bv_h_A_len + Bv_t_A_len + Bh_A2_len + Bv_A2_len + LI_len 
        nom  = len(Bcr_h) + len(Bcr_v) +len(Ag) + len(Lfa) + len(Icam) + len(Bh_h_A) + len(Bh_t_A) + len(Bv_h_A) + len(Bv_t_A) + len(Bh_A2) + len(Bv_A2) + len(LI)

        #iterating for nom times
        for mol in range(nom):

                moltype_checker = nom*np.random.random()

                cumu_prob = nom - len(Lfa)
                #cumu_prob = len(Lfa)

                if moltype_checker > cumu_prob:

                        molecule = Lfa[np.random.randint(0,len(Lfa) - 1)]
                        #molecule = np.random.choice(Lfa) #Change np.randomchoice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                LIA_diffusion(molecule,Lfa)

                        else:

                                Lfa_reaction(molecule)
                        
                #cumu_prob -= len(Icam)


                elif moltype_checker > cumu_prob - len(Icam):

                        molecule = Icam[np.random.randint(0,len(Icam) - 1)]
                        
                        #molecule = np.random.choice(Icam) #Change random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                LIA_diffusion(molecule,Icam)

                        else:

                                Icam_reaction(molecule)


                
                #cumu_prob -= len(Bcr_v)
                

                elif moltype_checker > cumu_prob - len(Icam) - len(Bcr_v):

                        molecule = Bcr_v[np.random.randint(0,len(Bcr_v) - 1)]

                        #molecule = np.random.choice(Bcr_v) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bcr_v_diffusion(molecule)

                        else:

                                Bcr_v_reaction(molecule)


                #cumu_prob -= len(Bcr_h)
                

                elif moltype_checker > cumu_prob - len(Icam) - len(Bcr_v) - len(Bcr_h):

                        molecule = Bcr_h[np.random.randint(0,len(Bcr_h) - 1)]

                        #molecule = np.random.choice(Bcr_h) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bcr_h_diffusion(molecule)

                        else:

                                Bcr_h_reaction(molecule)

                #cumu_prob -= len(Ag)
                

                elif moltype_checker > cumu_prob - len(Icam) - len(Bcr_v) - len(Bcr_h) - len(Ag):

                        molecule = Ag[np.random.randint(0,len(Ag) - 1)]

                        #molecule = np.random.choice(Ag) #Change random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                LIA_diffusion(molecule,Ag)

                        else:

                                Ag_reaction(molecule)

                #cumu_prob -= len(Bh_A2)
                

                elif moltype_checker > cumu_prob - len(Icam) - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI):

                        molecule = LI[np.random.randint(0,len(LI) - 1)]

                        #molecule = np.random.choice(LI) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                LI_diffusion(molecule)

                        else:

                                LI_reaction(molecule)

                #cumu_prob -= len(Bh_h_A)
                

                elif moltype_checker > cumu_prob - len(Icam) - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A):

                        molecule = Bh_h_A[np.random.randint(0,len(Bh_h_A) - 1)]

                        #molecule = np.random.choice(Bh_h_A) #Change np.randomchoice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bh_h_A_diffusion(molecule)

                        else:

                                Bh_h_A_reaction(molecule)


                #cumu_prob -= len(Bh_t_A)
                

                elif moltype_checker > cumu_prob - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A) - len(Bh_t_A):

                        molecule = Bh_t_A[np.random.randint(0,len(Bh_t_A) - 1)]

                        #molecule = np.random.choice(Bh_t_A) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bh_t_A_diffusion(molecule)

                        else:

                                Bh_t_A_reaction(molecule)

                #cumu_prob -= len(Bv_h_A)
                

                elif moltype_checker > cumu_prob - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A) - len(Bh_t_A) - len(Bv_h_A):

                        molecule = Bv_h_A[np.random.randint(0,len(Bv_h_A) - 1)]

                        #molecule = np.random.choice(Bv_h_A) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bv_h_A_diffusion(molecule)

                        else:

                                Bv_h_A_reaction(molecule)

                #cumu_prob -= len(Bv_t_A)
                

                elif moltype_checker > cumu_prob - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A) - len(Bh_t_A) - len(Bv_h_A) - len(Bv_t_A):

                        molecule = Bv_t_A[np.random.randint(0,len(Bv_t_A) - 1)]

                        #molecule = np.random.choice(Bv_t_A) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bv_t_A_diffusion(molecule)

                        else:

                                Bv_t_A_reaction(molecule)

                #cumu_prob -= len(Bh_A2)
                

                elif moltype_checker > cumu_prob - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A) - len(Bh_t_A) - len(Bv_h_A) - len(Bv_t_A) - len(Bh_A2):

                        molecule = Bh_A2[np.random.randint(0,len(Bh_A2) - 1)]

                        #molecule = np.random.choice(Bh_A2) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bh_A2_diffusion(molecule)

                        else:

                                Bh_A2_reaction(molecule)


                elif moltype_checker > cumu_prob - len(Bcr_v) - len(Bcr_h) - len(Ag) - len(LI) - len(Bh_h_A) - len(Bh_t_A) - len(Bv_h_A) - len(Bv_t_A) - len(Bh_A2) - len(Bv_A2):

                        molecule = Bv_A2[np.random.randint(0,len(Bv_A2) - 1)]
                        
                        #molecule = np.random.choice(Bv_A2) #Change np.random.choice into something elseeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

                        if np.random.random() > 0.5:

                                Bv_A2_diffusion(molecule)

                        else:

                                Bv_A2_reaction(molecule)

        

        synba.append(synform(Bv_A2+Bh_A2+Bv_t_A+Bv_h_A+Bh_t_A+Bh_h_A))
        synli.append(synform(LI))


stop = timer()
print(stop - start)

#print(surf)
#print(np.count_nonzero(surf))
"""
for lisi in lis:

        for mol in lisi:

                x,y,z = mol
                if surf[i][j][k] == 0:

                        print(mol,lisi)
        
"""
