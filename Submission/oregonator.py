import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import time as tm

class Elementary_reaction(object):
    def __init__(self,reactants,products, k):
        "Assumes reactants and products being list of manes of species"
        self.reactants=reactants
        self.products=products
        self.k=k
    def get_reactants(self):
        return self.reactants
    def get_products(self):
        return self.products
    def get_species(self):
        return self.reactants+self.products
    def get_rate_const(self):
        return self.k
    def set_rate_const(self,new_k):
        self.k=new_k

class Reactions(object):
    "Dictionary with each cpesies being a key to it's updated concentration and also contains the reactions"
    def __init__(self):
        self.species={}
        self.rx=[]
    def add_species(self,SP):
        """Add a species to the system"""
        if SP in self.species:
            raise ValueError('Duplicate species')
        else:
            self.species[SP]=0
    def get_species(self):
        """returns a dictionary with the names of the species and their concentration"""
        return self.species
    def set_concentration(self,SP,conc):
        """sets the concentration of the species SP"""
        self.species[SP]=conc
    def add_elem_rx(self,elem_rx):
        """Adds an elementary reaction under every SP that participates in it"""
        for SP in elem_rx.get_reactants():
            if not (SP in self.species):
                raise ValueError("Elementary reaction involves unknown species "+SP)
        for SP in elem_rx.get_products():
            if not (SP in self.species):
                raise ValueError("Elementary reaction involves unknown species "+SP)
        self.rx.append(elem_rx)
    def get_elem_rx(self):
        """Returns a list of the elementary reactions"""
        return self.rx
    
    def rate_of_change(self,time_step):
        """Euler numerical integrator"""
        i=0
        prod=[]
        for xx in self.rx:
            prod.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                prod[i]=prod[i]*self.species[yy]
            i+=1
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                self.species[sp]-=time_step*prod[j]
            for sp in self.rx[j].get_products():
                self.species[sp]+=time_step*prod[j]        

    def rate_of_change_RK4(self,time_step):
        i=0
        prod=[]
        for xx in self.rx:
            prod.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                prod[i]=prod[i]*self.species[yy]
            i+=1
        d1={}
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                if sp not in d1:
                    d1[sp]=0
                d1[sp]-=time_step*prod[j]
            for sp in self.rx[j].get_products():
                if sp not in d1:
                    d1[sp]=0
                d1[sp]+=time_step*prod[j]
        i=0
        prod2=[]
        for xx in self.rx:
            prod2.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                prod2[i]=prod2[i]*(self.species[yy]+0.5*d1[yy])
            i+=1
        d2={}
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                if sp not in d2:
                    d2[sp]=0
                d2[sp]-=time_step*prod2[j]
            for sp in self.rx[j].get_products():
                if sp not in d2:
                    d2[sp]=0
                d2[sp]+=time_step*prod2[j]
        i=0
        prod3=[]
        for xx in self.rx:
            prod3.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                prod3[i]=prod3[i]*(self.species[yy]+0.5*d2[yy])
            i+=1
        d3={}
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                if sp not in d3:
                    d3[sp]=0
                d3[sp]-=time_step*prod3[j]
            for sp in self.rx[j].get_products():
                if sp not in d3:
                    d3[sp]=0
                d3[sp]+=time_step*prod3[j]
        i=0
        prod4=[]
        for xx in self.rx:
            prod4.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                prod4[i]=prod4[i]*(self.species[yy]+d3[yy])
            i+=1
        d4={}
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                if sp not in d4:
                    d4[sp]=0
                d4[sp]-=time_step*prod4[j]
            for sp in self.rx[j].get_products():
                if sp not in d4:
                    d4[sp]=0
                d4[sp]+=time_step*prod4[j]
        for sp in self.species:
            self.species[sp]+=d1[sp]/6+d2[sp]/3+d3[sp]/3+d4[sp]/6


      
def specify_system(file,rx_system_name):
    "builds the reaction system from an appropriate file"
    f=open(file,"r")
    linenumb=0
    for line in f:
        if linenumb==0:
            l=line.split()
            for spec in l:
                rx_system_name.add_species(spec)
            linenumb+=1
        else:
            l=line.split()
            reacta=[]
            produ=[]
            i=0
            while l[i]!="-":
                reacta.append(l[i])
                i+=1
            for j in range(i+1,len(l)-2):
                produ.append(l[j])
            const=float(l[-1])
            elemen_rx=Elementary_reaction(reacta,produ,const)
            rx_system_name.add_elem_rx(elemen_rx)

class Line_of_cells(object):
    """Contains a 1D line of Reactions objects as cells"""
    def __init__(self):
        self.cells=[]
    def add_cell(self,cell):
        """Add a new cell at the end of the line"""
        self.cells.append(cell)
    def diffusion(self,cell_one,cell_two,cell_three,time_step,diffusion_coef,delta_x):
        "expects to reactions object having the same species present and diffuses between them"
        for sp in cell_one.get_species():
            lapl=(cell_one.species[sp]+cell_three.species[sp]-2*cell_two.species[sp])/delta_x**2
            cell_one.species[sp]+=time_step*diffusion_coef*lapl
    def approximate_diffusion(self,cell_one,cell_two,diffusion_coef,time_step):
        """Diffuese the species based on an approxiamte description, where rate is just proportional to the concentration difference"""
        for sp in cell_one.get_species():
            conc_change=time_step*diffusion_coef*(cell_two.species[sp]-cell_one.species[sp])
            cell_one.species[sp]+=conc_change
            cell_two.species[sp]-=conc_change            

# Set up the system with hte appropriate initial concentrations

oregonator=Reactions()
specify_system("oregonator.txt",oregonator)
oregonator.set_concentration("A",0.06)
oregonator.set_concentration("B",0.06)
oregonator.set_concentration("P",0.0)
oregonator.set_concentration("Q",0.0)
oregonator.set_concentration("X",10**(-9.8))
oregonator.set_concentration("Y",10**(-6.52))
oregonator.set_concentration("Z",10**(-7.32))

#Simulate the oregonator reaction for 90 seconds
X=[]
Y=[]
Z=[]
time=[]
delta_t=2*10**(-6)
t0=tm.time()
for i in range(0,45000000):
    oregonator.rate_of_change(delta_t)
    if i%130000==0:
        X.append(oregonator.get_species()["X"])
        Y.append(oregonator.get_species()["Y"])
        Z.append(oregonator.get_species()["Z"])
        time.append(i*delta_t)
t1=tm.time()
print("The calculation took "+ str(t1-t0)+" seconds.")
np.save("oregonator_conc_x_larger_t.npy",X)
np.save("oregonator_y_larger_t.npy",Y)
np.save("oregonator_z_larger_t.npy",Z)
np.save("oregonator_t_larger_t.npy",time)

plt.semilogy(time,X,label="X",marker="o")
plt.semilogy(time,Y,label="Y")
plt.semilogy(time,Z,label="Z",marker="x")
plt.legend()
plt.xlim((time[0],time[-1]))
plt.xlabel("time / s")
plt.ylabel("Concentration")
plt.savefig('oregonator_kinetics.pdf')
plt.show()
