import sys,os
from gmx_top import *

class gmx_pdb:
    def __init__(self,filename='NONE'):
        self.box     = [0.0,0.0,0.0]
        if filename != 'NONE':
            self.atoms   = self.read_pdb(filename)
            self.residues= self.assign_residues()
        else:
            self.atoms      = []
            self.residues   = []

    class pdb_atom: #creates blank atom 
        def __init__(self,nr=0,name='',res='',seq=0,x=0.0,y=0.0,z=0.0,occup=0.0,tempfac=0.0,element=''):
            self.nr      = ''
            self.name    = ''
            self.res     = ''
            self.seq     = ''
            self.x       = ''
            self.y       = ''
            self.z       = ''
            self.occup   = '' 
            self.tempfac = ''
            self.element = ''
            self.qm_group= -1
            self.bonds   = {}
            
        def s(self,string,length): # prepares string 'string' with length 'length'
            string = str(string)
            while len(string)<length: string=' ' + string;
            if len(string)>length: string=string[0:length];
            return string
        def cprint(self): #creates string in correct PDB format version 3.2 
            return self.s('ATOM  ',6) + self.s(self.nr,5) + ' ' + self.s(self.name,4) + ' ' + self.s(self.res,3) + '  ' + self.s(self.seq,4) + '    ' + self.s(self.x,8) + self.s(self.y,8) + self.s(self.z,8) + self.s(self.occup,6) + self.s(self.tempfac,6) + self.s('',10) + self.s(self.element,2)  
        def mprint(self): #creates pdbobject lines
            return ['ATOM',str(self.nr),' '+str(self.name)+' ',str(self.res),str(self.seq),str(self.x),str(self.y),str(self.z),str(self.occup),str(self.tempfac),str(self.res)]
         
    def read_pdb(self,filename):
        print '\n..Reading PDB file',filename
        pdb_atoms=[]        
        for line in open(filename,'r').readlines():
            if line.find('CRYST1')>=0:#try to find a box definition is pdb file
                try:
                    self.box = [ float(line[9:15].strip(' ')), float(line[18:24].strip(' ')),float(line[27:33].strip(' '))]
                except Exception:0;
                
            if line[:6].find('ATOM')>=0:
                atom = self.pdb_atom()
                try:atom.nr      = (int)    (line[6:11].strip(' ' ));
                except Exception:  print 'Error reading PDB file',filename,'atom nr. missing\n-->',line; sys.exit(1); 
                try:
                    atom.name    =          (line[12:16].strip(' '));
                    #here, the water bonds are added, since these do not appear in the topol.top file
                    if atom.name == 'OW':
                        atom.bonds[atom.nr]=[atom.nr+1,atom.nr+2]
                    elif atom.name == 'HW1':
                        atom.bonds[atom.nr]=[atom.nr-1]
                    elif atom.name == 'HW2':
                        atom.bonds[atom.nr]=[atom.nr-2]
                except Exception:  print 'Error reading PDB file',filename,'atom name. missing\n-->',line; sys.exit(1);
                try:atom.res     =          (line[17:20].strip(' '));
                except Exception:  print 'Error reading PDB file',filename,'atom res. missing\n-->',line; sys.exit(1);
                try:atom.seq     = (int)    (line[22:26].strip(' '));
                except Exception: 0 ;
                try:
                    atom.x       = (float)  (line[30:38].strip(' '));
                    atom.y       = (float)  (line[38:46].strip(' '));
                    atom.z       = (float)  (line[46:54].strip(' '));
                except Exception:  print 'Error reading PDB file',filename,'coordinates missing\n-->',line; sys.exit(1);
                try:atom.occup   = (float)  (line[54:60].strip(' '));
                except Exception: 0;
                try:atom.tempfac = (float)  (line[60:66].strip(' '));
                except Exception: 0;
                try:atom.element =          (line[76:78].strip(' '));
                except Exception: 0;
                pdb_atoms.append(atom)
        print '..Identified',len(pdb_atoms),'atoms from PDB file\n'
        return pdb_atoms
    
    def assign_residues(self):
        #creates a residues[] list which contains residue[] lists
        #each residue list contains pdb_atom() objects
        residues = []
        residue  = []
        first_run= 1
        curr_atom= ''
        last_atom= self.atoms[0].name
        curr_res = ''
        last_res = self.atoms[0].res
        for atom in self.atoms:
            curr_res = atom.res
            curr_atom= atom.name
            if curr_res!=last_res or (curr_atom==last_atom and first_run==0):
                residues.append(residue)
                residue=[]
                last_res=curr_res
                last_atom=curr_atom
            if curr_res==last_res :
                residue.append(atom)
            if first_run==1:
                first_run=0
        residues.append(residue)
        #for i in residues:print i[0].res;
        return residues
    
    def write_pdb(self,filename):
        #since each pdb_atom carries a cprint() method, the PDB file can be written very easily
        f=open(filename,'w')
        for atom in self.atoms:
            print >> f, atom.cprint()