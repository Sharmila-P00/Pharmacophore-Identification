from functional_groups import groups
from rdkit import Chem
from collections import defaultdict


class parameters:
    
    def descriptor(self,file,descriptor_group):
        results=defaultdict(list)
        again_dict={}
        sdf_to_mol=Chem.SDMolSupplier(file,removeHs=False)
        for mol in sdf_to_mol:
            if mol is None:
                continue
            for name,smarts in descriptor_group.items():
                for groups in smarts:
                    pattern=Chem.MolFromSmarts(groups)
                    if mol.HasSubstructMatch(pattern):
                            match=mol.GetSubstructMatches(pattern)
                            for i in match:
                                results[name].append(i)
                                
            results=dict(results)

            for key,values in results.items():
                for j,val in enumerate(values,start=1):
                    new_ke=f'{key} {j if len(values) > 1 else ""}'
                    again_dict[new_ke]=val
                
            return again_dict
        




