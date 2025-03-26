from rdkit import Chem
from features import parameters
from functional_groups import groups



Param=parameters()
features=Param.descriptor(file,groups())

class coordinates:
    
    def coord(self,file):
        co_ord={}
        
        sdf_to_mol=Chem.SDMolSupplier(file,removeHs=False)
        for mol in sdf_to_mol:
            if mol is None:
                continue
            
            for i in range(mol.GetNumAtoms()):
                conf=mol.GetConformer()
                r=conf.GetAtomPosition(i)
                co_ord[i]=(r.x,r.y,r.z)
        return co_ord

    def fea_coord(self,co_ord,features):    
        atm_coord={}
        feature_coord={} 
        for key,value in features.items():
            atm_coord[key]=[co_ord[i] for i in value if i in co_ord]
                  
        for key,value in atm_coord.items():
            summ= tuple(sum(x) for x in zip(*value))
            feature_coord[key]=tuple(x/len(value) for x in summ)    
        return feature_coord
    




