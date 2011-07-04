from gmx_pdb import gmx_pdb
from gmx_top import gmx_topology
from math import *
from copy import *

class gentop ():
    def __init__(self,pdb='mol.pdb',top='topol.top'):
        self.pdb        = gmx_pdb(pdb)
        self.top        = gmx_topology()
        self.top.read_top('topol.top')
        self.zoneatoms  = []
        self.read_zoneatoms('qm_atomlist.dat')
        self.assign_bonds()
        self.qm_connect = self.assign_qmgroups()  # list of [ [groupmembers],[partners] ] with [partners]=[[groupmember],[outsideatom]] so [ [1,2,3], [[2,5],[3,7]]
        self.azone,self.zone_id_list,self.azone_center_id=self.find_azone()
        self.create_zones_top_working()

    def read_zoneatoms(self,filename):
        z=[]
        for line in open(filename,'r'):
            z.append(int(line))
        self.zoneatoms=z
        print z
    def find_azone(self):
        #creates and completes a zone of radius r around a center atom
        zone_id_list={}
        azone_center_id=-1
        azone={}
        for i in self.zoneatoms:    #cuts QM_zone
            azone[i] = ''

        #-----writes_pdb_for_illustrational_purposes-----#
        f = open('azone_pre.pdb','w')
        tmp=[]
        for i in self.pdb.atoms:
            if i.nr in azone.keys():
                tmp.append(i)
        for i in tmp:
            print >>f, i.cprint()
        f.close()
        #-----END-----#

        #here, the QM_zone is completed (or T-Zone in the adaptive case)
        for key in azone.keys():        #loops over QM_zone atom numbers
            atm = self.pdb.atoms[key-1] #call pdb atom with matching atom number
            if atm.qm_group!= -1 :      #exclude all backbone atoms (also all PRO & GLY residues)
                zone_id_list[atm.qm_group]='' #creates a list of qm_group ids which are in the A&T zone
                k = self.qm_connect[atm.qm_group][0]
                for atom in k:
                        #if atom==center_atom_number:
                        #    azone_center_id=atm.qm_group    #gets the qm_group id of the A-zone qm_group
                        azone[atom]=''
        
        #-----writes_pdb_for_illustrational_purposes-----#
        f = open('azone_post.pdb','w')
        tmp=[]
        
        for atm in azone.keys():
            tmp.append(self.pdb.atoms[atm-1])
        for i in tmp:
            print >>f, i.cprint()
        f.close()
        #-----END-----#
        return azone,zone_id_list,azone_center_id

    def assign_qmgroups(self):
        resgroup=[]
        for res in xrange(len(self.pdb.residues)):
            tmpgroup=[]
            for atm in self.pdb.residues[res]:
                if atm.name not in ['N','H','C','O','CA','HA','LA'] and atm.qm_group==-1:#'OW','HW1','HW2'
                    tmpgroup.append(atm.nr)
                #    if atm.name not in ['OW','HW1','HW2']:print atm.nr,atm.name #<-- here, everything still works ;-)
            for i in tmpgroup:
                self.pdb.atoms[i-1].qm_group=len(resgroup)
            if tmpgroup:
                resgroup.append(tmpgroup)
            
        qmconnect=[]
        for group in resgroup:
            partners = []
            partners_bonds =[]
            for atm in group:
                try:
                    for partner in self.pdb.atoms[atm-1].bonds[atm]:
                        if partner not in group:
                            partners.append(partner) #####################
                except Exception:a=0;
            qmconnect.append([group,partners])
        return qmconnect
    
        
    def assign_bonds(self):
        for at in self.top.bonds_list:
            try:
                i =(at.ai)
                ip=(at.ai-1)
                j =(at.aj)
                self.pdb.atoms[ip].bonds[i].append(j)
            except Exception:
                self.pdb.atoms[ip].bonds[i] = [j]
            try:
                j =(at.aj)
                jp=(at.aj-1)
                i =(at.ai)
                self.pdb.atoms[jp].bonds[j].append(i)
            except Exception:
                self.pdb.atoms[jp].bonds[j] = [i]
    def create_zones_top_working(self):
        for zone in ['TEST']: #get all atoms in all zones
            top = deepcopy(self.top)                        ######## THIS IS EXTREMELY SLOW... here one should us an 'undo' function 
            pdb = deepcopy(self.pdb)                        ######## THIS IS EXTREMELY SLOW... here one should us an 'undo' function
            qm_zone_atoms = self.zoneatoms
            links=[[987,988],[961,959],[3208,3206],[3485,3482]]
            #==TOPology==# bonds section fix
            for bond in top.bonds_list: #set all the bonds in the qm zone to function type 5
                if bond.ai in qm_zone_atoms: 
                    bond.funct = 5
                if bond.aj in qm_zone_atoms:
                    bond.funct = 5
                for link in links: # remove explicit bonds between qm and mm zone
                    if (bond.ai==link[0] and bond.aj==link[1]) or (bond.ai==link[1] and bond.aj==link[0]):
                        bond.funct = -1
            #==TOPology==# atoms section fix & vsites section fix & constraints section fix
            top.virtual_sites2_meta.append('[ virtualsites2 ]')             
            top.constraints_meta.append('[ constraints ]')                  
            for i in links:
            #pdb_file
                tmp4 = pdb.pdb_atom()
                tmp4.nr=str(len(top.atoms_list)+1)
                tmp4.name='LA'
                tmp4.res='XXX'
                tmp4.seq=str(top.atoms_list[-1].resnr+1)
                tmp4.x='0.0'
                tmp4.y='0.0'
                tmp4.z='0.0'
                tmp4.occup='1.0'
                tmp4.tempfac='0.0'
                pdb.atoms.insert(len(top.atoms_list),tmp4)
            #linkatoms
                tmp = top.atoms(\
                               nr       = str(len(top.atoms_list)+1) ,\
                               type     = 'LA',\
                               resnr    = str(top.atoms_list[-1].resnr+1),\
                               residue  = 'XXX',\
                               atom     = 'LA',\
                               cgnr     = str(top.atoms_list[-1].cgnr+1),\
                               charge   = '0.0',\
                               mass     = '0.0' )
                top.atoms_list.append(tmp)
            #vsites
                tmp2 = top.virtual_sites2(\
                                site    = str(len(top.atoms_list)),\
                                qmatom  = i[0],\
                                mmatom  = i[1],\
                                type    = 1,\
                                a       = 0.65)
                top.virtual_sites2_list.append(tmp2)
            #constraints
                tmp3 =top.constraints(\
                                qmatom  = i[0],\
                                mmatom  = i[1],\
                                type    = 2,\
                                a       = 0.153)
                top.constraints_list.append(tmp3)
            print 'writing:','link_'+str(zone)+'.top/.pdb'
            top.write_top('link_'+str(zone)+'.top')
            pdb.write_pdb('link_'+str(zone)+'.pdb')

gentop(pdb='mol.pdb',top='topol.top')
