def s(string,length):
    string = str(string)
    while len(string)<length: string=' ' + string;
    #if len(string)>length: string=string[0:length];
    return string
class gmx_topology:
    def __init__(self):  
        #->Here, arrays are created. Each array element consists of a 
        #->class that stores different elements from the .top file (a bond, dihedral, angle, etc. )
        #->The possible classes are listed here, and will be defined below: 
        self.defaults_list       = []
        self.atomtypes_list      = []
        self.moleculetype_list   = []
        self.atoms_list          = []
        self.pairs_list          = []
        self.bonds_list          = []
        self.angles_list         = []
        self.dihedrals_list      = []
        self.impropers_list      = []
        self.system_list         = []
        self.molecules_list      = []
        self.virtual_sites2_list = []
        self.constraints_list    = []
        #-> these arrays store all the comment and descriptive lines,
        #-> e.g. '[ bonds ]' or ';this is a comment'
        self.defaults_meta      = []
        self.atomtypes_meta     = []
        self.moleculetype_meta  = []
        self.atoms_meta         = []
        self.pairs_meta         = []
        self.bonds_meta         = []
        self.angles_meta        = []
        self.dihedrals_meta     = []
        self.impropers_meta     = []
        self.system_meta        = []
        self.molecules_meta     = []
        self.virtual_sites2_meta= []
        self.constraints_meta   = []
    #There are two methods defined for this topology class, read_top and write_top
    #read_top uses a lot of .strip(), .remove(), .split() operations, which were
    #tested on several topologies and should do the job.
    def read_top(self,filename):
        block=[]
        blocklist=[]
        #Read Topology and separate it into blocks [ atoms ], [ bonds ], etc.
        for i in open(filename,'r'):
            if len(i.strip('\n'))==0:    # blank line indicates the end of a block [ atoms ], [ bonds ], etc.
                if len(block)>0: blocklist.append(block); 
                block = []
            elif len(i.strip('\n'))>0:    # read block
                block.append(i.strip('\n'))
        blocklist.append(block);

        for i in range(blocklist.count([])): blocklist.remove([]);
        #Read [ atomtypes ]
        for block in blocklist:
                if block[0].find('atomtypes') >= 0:
                    for i in range(len(block)): #lines
                        if block[i][0]==';' or block[i][0]=='[':
                            self.atomtypes_meta.append(block[i])           #save comment lines, [...], etc.
                        else:
                            top=block[i].split(' ')
                            for i in range(top.count('')): top.remove(''); #remove blanks from array
                            self.atomtypes_list.append(self.atomtype(top[0],top[1],top[2],top[3],top[4],top[5],top[6]))

        #Read [ moleculetype ]
        for block in blocklist:
            if block[0].find('moleculetype') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.moleculetype_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        self.moleculetype_list.append(self.moleculetype(top[0],top[1]))
        
        #Read [ atoms ]
        for block in blocklist:
            if block[0].find('atoms') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.atoms_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #nr,type,resnr,residue,atom,cgnr,charge,mass
                        self.atoms_list.append(self.atoms(top[0],top[1],top[2],top[3],top[4],top[5],top[6],top[7]))
        
        #Read [ bonds ]
        for block in blocklist:
            if block[0].find('bonds') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.bonds_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #ai,aj,funct,r,k
                        try:
                            self.bonds_list.append(self.bonds(top[0],top[1],top[2],top[3],top[4]))
                        except Exception:   
                            self.bonds_list.append(self.bonds(top[0],top[1],top[2]))

        #Read [ pairs ]
        for block in blocklist:
            if block[0].find('pairs') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.pairs_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #ai,aj,funct,r,k
                        self.pairs_list.append(self.pairs(top[0],top[1],top[2]))
    
        #Read [ angles ]
        for block in blocklist:
            if block[0].find('angles') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.angles_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #ai,aj,ak,funct,theta,cth
                        try:
                            self.angles_list.append(self.angles(top[0],top[1],top[2],top[3],top[4],top[5]))
                        except Exception:
                            self.angles_list.append(self.angles(top[0],top[1],top[2],top[3]))

        #Read [ dihedrals ]
        for block in blocklist:
            if block[0].find('dihedrals') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.dihedrals_meta.append(block[i])           #save comment lines, [...], etc.
                        self.impropers_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')                        #split line
                        for i in range(top.count('')): top.remove(''); #remove empty elements
                        #i,j,k,l,funct,C0,C1,C2,C3,C4,C5
                        if int(top[4])==3:
                            try:
                                self.dihedrals_list.append(self.dihedrals(top[0],top[1],top[2],top[3],top[4],top[5],top[6],top[7],top[8],top[9],top[10]))
                            except Exception:
                                self.dihedrals_list.append(self.dihedrals(top[0],top[1],top[2],top[3],top[4]))
                        if int(top[4])==1:
                            self.impropers_list.append(self.impropers(top[0],top[1],top[2],top[3],top[4]))
        
        #Read [ system ]
        for block in blocklist:
            if block[0].find('system') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.system_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #nr, name 
                        try: 
                            self.system_list.append(self.system(top[0],top[1]))
                        except Exception: 1;
        #Read [ molecules ]
        for block in blocklist:
            if block[0].find('molecules') >= 0:
                for i in range(len(block)): #lines
                    if block[i][0]==';' or block[i][0]=='[':
                        self.molecules_meta.append(block[i])           #save comment lines, [...], etc.
                    else:
                        top=block[i].split(' ')
                        for i in range(top.count('')): top.remove(''); #remove blanks from array
                        #nr, name 
                        try: 
                            self.molecules_list.append(self.molecules(top[0],top[1]))
                        except Exception: 1;
    #Each of the topology elements contains a method called cprint() which contains 
    #the proper printing format for it's type, this makes the write_top function 
    #very simple
    def write_top(self,filename):
        f=open(filename,'w')
        print >> f, ';topology written using timo\'s topology tool ttt.py - tgraen@gwdg.de\n'
        #write [ defaults ]
        for i in self.defaults_meta: print >> f,i; 
        for i in self.defaults_list: print >> f,i.cprint();
        print >> f, ''
        #write [ atomtypes ]
        for i in self.atomtypes_meta: print >> f,i;
        for i in self.atomtypes_list:print >> f,i.cprint();
        print >> f, ''
        #write [ moleculetype ]
        for i in self.moleculetype_meta: print >> f,i;
        for i in self.moleculetype_list:print >> f,i.cprint();
        print >> f, ''
        #write [ atoms ]
        for i in self.atoms_meta: print >> f,i;
        for i in self.atoms_list:print >> f,i.cprint();
        print >> f, ''
        #write [ bonds ]
        for i in self.bonds_meta: print >> f,i;
        for i in self.bonds_list:
            if i.funct>0:
                print >> f,i.cprint();
        print >> f, ''    
        #write [ pairs ]
        for i in self.pairs_meta: print >> f,i;
        for i in self.pairs_list:print >> f,i.cprint();
        print >> f, ''               
        #write [ angles ]
        for i in self.angles_meta: print >> f,i;
        for i in self.angles_list:print >> f,i.cprint();
        print >> f, ''               
        #write [ dihedrals ]
        for i in self.dihedrals_meta[:2]: print >> f,i;
        for i in self.dihedrals_list:print >> f,i.cprint();
        print >> f, ''               
        #write [ impropers ]
        for i in self.impropers_meta[2:]: print >> f,i;
        for i in self.impropers_list:print >> f,i.cprint();
        print >> f, ''               
        #write [ virtual_sites2 ]
        for i in self.virtual_sites2_meta: print >> f,i;
        for i in self.virtual_sites2_list:print >> f,i.cprint();
        print >> f, ''
        #write [ constraints ]
        for i in self.constraints_meta: print >> f,i;
        for i in self.constraints_list:print >> f,i.cprint();
        print >> f, ''
        #write [ system ]
        for i in self.system_meta: print >> f,i;
        for i in self.system_list:print >> f,i.cprint();
        print >> f, ''                
        #write [ molecules ]
        for i in self.molecules_meta: print >> f,i;
        for i in self.molecules_list:print >> f,i.cprint();
        print >> f, ''
        f.close()

    # classes for the data types inside the GMX Topology file:   
    class defaults:
        def __init__(self,nbfunc,comb_rule,gen_pairs,fudgeLJ,fudgeQQ):
            self.nbfunc    = (int)(nbfunc.strip(' '))
            self.comb_rule = (int)(comb_rule.strip(' '))
            self.gen_pairs = gen_pairs.strip(' ')
            self.fudgeLJ   = (float)(fudgeLJ.strip(' '))
            self.fudgeQQ   = (float)(fudgeQQ.strip(' '))
        def cprint(self):
            return str(self.nbfunc) + s(self.comb_rule,15) + s(self.gen_pairs,18) + s(self.fudgeLJ,16) + s(self.fudgeQQ,11)
                                                               
    class atomtype:
        def __init__(self,name,bond_type,mass,charge,ptype,sigma,epsilon):
            self.name    = name.strip(' ')
            self.bond_type= bond_type.strip(' ')
            self.mass    = (float)(mass.strip(' '))
            self.charge  = (float)(charge.strip(' '))
            self.ptype   = ptype.strip(' ')
            self.sigma   = (float)(sigma.strip(' '))
            self.epsilon = (float)(epsilon.strip(' '))
        def cprint(self):
            return s(str(self.name),2) + s(self.bond_type,14) + s(self.mass,12) + s(self.charge,8) + s(self.ptype,3) + s(self.sigma,14) + s(self.epsilon,13)
    
    class moleculetype:
        def __init__(self,name,nrexcl):
            self.name    = name.strip(' ')
            self.nrexcl  = (int)(nrexcl.strip(' '))
        def cprint (self):
            return str(self.name) + ' ' + s(self.nrexcl,13)
    
    class atoms:
        def __init__(self,nr,type,resnr,residue,atom,cgnr,charge,mass):
            self.nr      = (int)(nr.strip(' '))
            self.type    = type.strip(' ')
            self.resnr   = (int)(resnr.strip(' '))
            self.residue = residue.strip(' ')
            self.atom    = atom.strip(' ')
            self.cgnr    = (int)(cgnr.strip(' '))
            self.charge  = (float)(charge.strip(' '))
            self.mass    = (float)(mass.strip(' ').strip('\n'))
        def cprint(self):
            return s(self.nr,6) + s(self.type,11) + s(self.resnr,7) + s(self.residue,7) + s(self.atom,7) + s(self.cgnr,7) + s(self.charge,11) + s(self.mass,11)
    
    class bonds:
        def __init__(self,ai,aj,funct,r='',k=''):
            self.ai      = (int)(ai.strip(' '))
            self.aj      = (int)(aj.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.r       = (float)(r.strip(' '))
                self.k       = (float)(k.strip(' '))
            except Exception:
                self.r       = ''
                self.k       = ''
                
        def cprint(self):
            return s(self.ai,5) + s(self.aj,6) + s(self.funct,6) + s(self.r,12) + s(self.k,12)
            
    class pairs:
        def __init__(self,ai,aj,funct):
            self.ai      = (int)(ai.strip(' '))
            self.aj      = (int)(aj.strip(' '))
            self.funct   = (int)(funct.strip(' '))
        def cprint(self):
            return s(self.ai,6) + s(self.aj,7)+ s(self.funct,7)
    
    class angles:
        def __init__(self, ai,aj,ak,funct,theta='',cth=''):
            self.ai      = (int)(ai.strip(' '))
            self.aj      = (int)(aj.strip(' '))
            self.ak      = (int)(ak.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.theta   = (float)(theta.strip(' '))
                self.cth     = (float)(cth.strip(' '))
            except Exception:
                self.theta   = ''
                self.cth     = ''
                
        def cprint(self):
            return s(self.ai,5) + s(self.aj,6) + s(self.ak,6) + s(self.funct,6) + s(self.theta,12) + s(self.cth,12)
        
    class dihedrals:
        def __init__(self,i,j,k,l,funct,C0='',C1='',C2='',C3='',C4='',C5=''):
            self.i       = (int)(i.strip(' '))
            self.j       = (int)(j.strip(' '))
            self.k       = (int)(k.strip(' '))
            self.l       = (int)(l.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.C0      = (float)(C0.strip(' '))
                self.C1      = (float)(C1.strip(' '))
                self.C2      = (float)(C2.strip(' '))
                self.C3      = (float)(C3.strip(' '))        
                self.C4      = (float)(C4.strip(' '))
                self.C5      = (float)(C5.strip(' ').strip('\t;'))
            except Exception:
                self.C0      = ''
                self.C1      = ''
                self.C2      = ''
                self.C3      = ''        
                self.C4      = ''
                self.C5      = ''
    
        def cprint(self):
            return s(self.i,5) + s(self.j,5) + s(self.k,5) + s(self.l,6)+ s(self.funct,6) + s(self.C0,12) + s(self.C1,12) + s(self.C2,12) + s(self.C3,12) + s(self.C4,12) + s(self.C5,12)
    
    class impropers:
        def __init__(self,i,j,k,l,funct,C0='',C1='',C2='',C3='',C4='',C5=''):
            self.i       = (int)(i.strip(' '))
            self.j       = (int)(j.strip(' '))
            self.k       = (int)(k.strip(' '))
            self.l       = (int)(l.strip(' '))
            self.funct   = (int)(funct.strip(' '))
            try:
                self.C0      = (float)(C0.strip(' '))
                self.C1      = (float)(C1.strip(' '))
                self.C2      = (float)(C2.strip(' '))
                self.C3      = (float)(C3.strip(' '))        
                self.C4      = (float)(C4.strip(' '))
                self.C5      = (float)(C5.strip(' ').strip('\t;'))
            except Exception:
                self.C0      = ''
                self.C1      = ''
                self.C2      = ''
                self.C3      = ''        
                self.C4      = ''
                self.C5      = ''
    
        def cprint(self):
            return s(self.i,5) + s(self.j,5) + s(self.k,5) + s(self.l,6)+ s(self.funct,6) + s(self.C0,12) + s(self.C1,12) + s(self.C2,12) + s(self.C3,12) + s(self.C4,12) + s(self.C5,12)
            
            
    class system:
        def __init__(self,nr,name):
            self.nr      = (int)(nr.strip(' '))
            self.name    = name.strip(' ')
        def cprint(self):
            return str(self.nr) + ' ' + str(self.name)
            
    class molecules: 
        def __init__(self,compound,nmols):
            self.compound= compound.strip(' ')
            self.nmols   = (int)(nmols.strip(' ')) 
        def cprint(self):
            return str(self.compound) + s(self.nmols,13)
    
    class virtual_sites2:
        def __init__(self,site,qmatom,mmatom,type,a):
            self.site   = site
            self.qmatom = qmatom
            self.mmatom = mmatom
            self.type   = type
            self.a      = a
        def cprint(self):
            return str(self.site) +' '+ str(self.qmatom) +' '+ str(self.mmatom) +' '+ str(self.type) +' '+ str(self.a) 

    class constraints:
        def __init__(self,qmatom,mmatom,type,a):
            self.qmatom = qmatom
            self.mmatom = mmatom
            self.type   = type
            self.a      = a
        def cprint(self):
            return str(self.qmatom) +' '+ str(self.mmatom) +' '+ str(self.type) +' '+ str(self.a) 
         