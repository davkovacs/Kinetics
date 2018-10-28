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
        """Returns the list of elementary reactions involved in the system"""
        return self.rx
    def rate_of_change_equilibrator(self,cutoff,time_step):
        """Updates the concentrations and returns whether the equilibrium is reached (within the cutoff)"""
        i=0
        prod=[]
        verdict=True
        concentrs0={}
        delta_c={}
        for xx in self.rx:
            prod.append(xx.get_rate_const())
            for yy in xx.get_reactants():
                if yy not in concentrs0:
                    concentrs0[yy]=self.species[yy]
                prod[i]=prod[i]*self.species[yy]
            i+=1
        for j in range(0,len(prod)):
            for sp in self.rx[j].get_reactants():
                self.species[sp]-=time_step*prod[j]
            for sp in self.rx[j].get_products():
                self.species[sp]+=time_step*prod[j]
        for sp in self.species:
            delta_c[sp]=self.species[sp]-concentrs0[sp]
            if np.abs(delta_c[sp])>cutoff:
                verdict=False
        return verdict

      
def specify_system(file,rx_system_name):
    "Builds the reaction system from an appropriate file"
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

protein=Reactions()
specify_system("protein.txt",protein)
protein.set_concentration("D",0.0)
protein.set_concentration("N",1.0)
urea=[0.25*x for x in range(0,33)]
delta_t=2*10**(-6)
cutoff=0.00000005
max_step=5000000
j=0
D=[]
I=[]
N=[]
k=[]
for rx in protein.get_elem_rx():
    k.append(rx.get_rate_const())
t0=tm.time()
for c in urea:
    D.append(0)
    I.append(0)
    N.append(0)
    protein.get_elem_rx()[0].set_rate_const(k[0]*np.exp(-1.68*c))
    protein.get_elem_rx()[1].set_rate_const(k[1]*np.exp(0.95*c))
    protein.get_elem_rx()[2].set_rate_const(k[2]*np.exp(-1.72*c))
    protein.get_elem_rx()[3].set_rate_const(k[3]*np.exp(1.2*c))
    judgement=False
    steps=0
    while judgement==False and steps<max_step:
        judgement=protein.rate_of_change_equilibrator(cutoff,delta_t) 
        steps+=1
    print(steps)
    print("Step number "+str(j)+" converged: "+ str(judgement))
    D[j]=protein.get_species()["D"]
    I[j]=protein.get_species()["I"]        
    N[j]=protein.get_species()["N"]        
    j+=1
t1=tm.time()
print("The calculation took: "+ str(t1-t0)+ " seconds")
np.save("protein_conc_d.npy",D)
np.save("protein_conc_i.npy",I)
np.save("protein_conc_n.npy",N)
np.save("protein_conc_urea.npy",urea)
plt.scatter(urea,D,label="D",marker="1")
plt.scatter(urea,I,label="I",marker="x")
plt.scatter(urea,N,label="N",marker="o")
plt.legend()
plt.xlim((urea[0]-0.05,urea[-1]+0.05))
plt.ylim((-0.1,1.1))
plt.xlabel("[Urea] / M")
plt.ylabel("Fraction of species")
plt.savefig('protein_folding_kinetics.pdf')
plt.show()
