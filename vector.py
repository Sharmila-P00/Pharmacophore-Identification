import numpy as np
from rdkit import Chem
from features import parameters
from functional_groups import groups    



Param=parameters()
features=Param.descriptor(file,groups())

class vector:    
    class HBD_vector:    
        def HBD(self,file,features):    
            neighbor_h={}
            sdf_to_mol=Chem.SDMolSupplier(file,removeHs=False)
            for mol in sdf_to_mol:
                if mol is None:
                    continue
            for key,value in features.items():
                feature_li=value[:]
                if key[0:3] == "HBD":
                    for z in value:
                        feature_li.append(z)
                        a=mol.GetAtomWithIdx(z)
                        for c in a.GetNeighbors():
                            if c.GetSymbol() == "H":
                                feature_li.append(c.GetIdx())
                neighbor_h[key]=feature_li
        
            neighbor_coord={}
            for key,value in neighbor_h.items():
                conf=mol.GetConformer()    
                neighbor_coord[key]=[conf.GetAtomPosition(i) for i in value]  
            
            for i in neighbor_coord.values():
                if len(i)==3:
                    vec_HBD2=[]
                    vec1={}
                    vec2={}
                    
                    base,pt1,pt2=i[0],i[1],i[2]
                    
                    v1=np.subtract(np.array(pt1),np.array(base))
                    v2=np.subtract(np.array(pt2),np.array(base))
                   
                    v1_sq=np.square(v1)
                    v2_sq=np.square(v2)
                    
                    add1=sum(v1_sq)
                    add2=sum(v2_sq)
                   
                    v1_mag=np.sqrt(add1)
                    v2_mag=np.sqrt(add2)
                    
                    x=np.divide(v1[0],v1_mag)
                    vec1['x']=x
                    y=np.divide(v1[1],v1_mag)
                    vec1['y']=y
                    z=np.divide(v1[2],v1_mag)
                    vec1['z']=z
                    x=np.divide(v2[0],v2_mag)
                    vec2['x']=x
                    y=np.divide(v2[1],v2_mag)
                    vec2['y']=y
                    z=np.divide(v2[2],v2_mag)
                    vec2['z']=z
                    vec_HBD2.append(vec1)
                    vec_HBD2.append(vec2)
                    
                    yield vec_HBD2
                    
                elif len(i)==2:
                    vec_HBD1=[]
                    vec01={}
                    base,pt1=i[0],i[1]

                    v1=np.subtract(np.array(pt1),np.array(base))
                    v1_sq=np.square(v1)
                    add1=sum(v1_sq)
                    v1_mag=np.sqrt(add1)
                    
                    x=np.divide(v1[0],v1_mag)
                    vec01['x']=x
                    y=np.divide(v1[1],v1_mag)
                    vec01['y']=y
                    z=np.divide(v1[2],v1_mag)
                    vec01['z']=z
                    vec_HBD1.append(vec01)
                    
                    yield vec_HBD1

    class HBA_vector:
        def HBA(self,file,features):
            neighbor_h={}
            supplier=Chem.SDMolSupplier(file,removeHs=False)
            for mol in supplier:
                if mol is None:
                    continue
    
            for key,value in features.items():
                feature_li=value[:]
                lis=[]
                if key[0:3] == "HBA":
                    for z in value:
                        feature_li.append(z)
                        a=mol.GetAtomWithIdx(z)
                        for c in a.GetNeighbors():
                            if c.GetSymbol() == "H":
                                feature_li.append(c.GetIdx())
                neighbor_h[key]=feature_li



            neighbor_coord={}
            for key,value in neighbor_h.items():
                conf=mol.GetConformer()    
                neighbor_coord[key]=[conf.GetAtomPosition(i) for i in value]
             
            for i in neighbor_coord.values():
                
                if len(i)==3:
                    vec1_2=[]
                    vec1={}
                    vec2={}
                    base,pt1,pt2=i[0],i[1],i[2] 
                      
                    v1=np.subtract(np.array(pt1),np.array(base))
                    v2=np.subtract(np.array(pt2),np.array(base))
            
                    v1_sq=np.square(v1)
                    v2_sq=np.square(v2)
                    
                    add1=sum(v1_sq)
                    add2=sum(v2_sq)
                    
                    v1_mag=np.sqrt(add1)
                    v2_mag=np.sqrt(add2)
                   
                    x=np.divide(v1[0],v1_mag)
                    vec1['x']=x
                    y=np.divide(v1[1],v1_mag)
                    vec1['y']=y
                    z=np.divide(v1[2],v1_mag)
                    vec1['z']=z
                    
                    x==np.divide(v2[0],v2_mag)
                    vec2['x']=x
                    y=np.divide(v2[1],v2_mag)
                    vec2['y']=y
                    z=np.divide(v2[2],v2_mag)
                    vec2['z']=z
                   
                    vec1_2.append(vec1)
                    vec1_2.append(vec2)
                
                    yield vec1_2
                    
                elif len(i)==2:
                    vec0=[]
                    vec={}
                    base,pt1=i[0],i[1]
                    
                    v1=np.subtract(np.array(pt1),np.array(base))
                    v1_sq=np.square(v1)
                    add1=sum(v1_sq)
                    v1_mag=np.sqrt(add1)
                    
                    x=np.divide(v1[0],v1_mag)
                    vec['x']=x
                    y=np.divide(v1[1],v1_mag)
                    vec['y']=y
                    z=np.divide(v1[2],v1_mag)
                    vec['z']=z
                    
                    vec0.append(vec)
                    
                    yield vec0
           
                if len(feature_li) == 1:
                    neighbor={}
                    for k in feature_li:
                        a = mol.GetAtomWithIdx(k)
                        for g in a.GetNeighbors():
                            lis.append(g.GetIdx())
                    neighbor[key] = lis

                    neighbor_coord={}
                    for key,value in neighbor.items():
                        conf=mol.GetConformer()    
                        neighbor_coord[key]=[conf.GetAtomPosition(i) for i in value]

                        
                    for i in neighbor_coord.values():
                        if len(i)==2:
                            vec0_1=[]
                            vec={}
                            base,pt1=i[0],i[1]
                                
                            v1=np.subtract(np.array(pt1),np.array(base))
                            v1_sq=np.square(v1)
                            add1=sum(v1_sq)
                            v1_mag=np.sqrt(add1)
                          
                            x=(np.divide(v1[0],v1_mag))*(-1)
                            vec['x']=x
                            y=(np.divide(v1[1],v1_mag))*(-1)
                            vec['y']=y
                            z=(np.divide(v1[2],v1_mag))*(-1)
                            vec['z']=z
                           
                            vec0_1.append(vec)
                            
                            yield vec0_1
    
    class Aromatic_vector:
        def aromatic(self,file,features):
            Aro_vec={}   
            sdf_to_mol=Chem.SDMolSupplier(file,removeHs=False)
            for mol in sdf_to_mol:
                if mol is None:
                    continue
            for key,value in features.items():
                new_li=value[:]
                if i[0:3] == "Aro":    
                    for z in value:  
                        new_li.append(z)
                Aro_vec[key]=new_li
                      
                                
            coord_aro={}
            for key,value in Aro_vec.items():
               conf=mol.GetConformer()
               coord_aro[key]=[conf.GetAtomPosition(i) for i in value]
                               
            for i in coord_aro.values():
                #print(len(i))
                if len(i)==6:
                    vec01_2=[]
                    vec1={}
                    vec2={}
                    base,pt1,pt2=i[0],i[2],i[4] 
                                  
                    v1=np.subtract(np.array(pt1),np.array(base))
                    v2=np.subtract(np.array(pt2),np.array(base))
                    
                    cross_pdt=np.cross(v1,v2)
                   
                    v_sq=np.square(cross_pdt)
                 
                    add1=sum(v_sq)
                   
                    v_mag=np.sqrt(add1)
                   
                    x=np.divide(cross_pdt[0],v_mag)
                    vec1['x']=x
                    y=np.divide(cross_pdt[1],v_mag)
                    vec1['y']=y
                    z=np.divide(cross_pdt[2],v_mag)
                    vec1['z']=z
                  
                    x=(np.divide(cross_pdt[0],v_mag))*(-1)
                    vec2['x']=x
                    y=(np.divide(cross_pdt[1],v_mag))*(-1)
                    vec2['y']=y
                    z=(np.divide(cross_pdt[2],v_mag))*(-1)
                    vec2['z']=z
                   
                    vec01_2.append(vec1)
                    vec01_2.append(vec2)
                    yield vec01_2
                    
                    
                elif len(i)==5:
                    vec11_2=[]
                    vec1_1={}
                    vec1_2={}
                   
                    base,pt1,pt2=i[0],i[2],i[3]
                                   
                    v1=np.subtract(np.array(pt1),np.array(base))
                    v2=np.subtract(np.array(pt2),np.array(base))
                          
                    cross_pdt=np.cross(v1,v2)
                   
                    v_sq=np.square(cross_pdt)
                   
                    add1=sum(v_sq)
                   
                    v_mag=np.sqrt(add1)
                   
                    x=np.divide(cross_pdt[0],v_mag)
                    vec1_1['x']=x
                    y=np.divide(cross_pdt[1],v_mag)
                    vec1_1['y']=y
                    z=np.divide(cross_pdt[2],v_mag)
                    vec1_1['z']=z
                    
                    x=(np.divide(cross_pdt[0],v_mag))*(-1)
                    vec1_2['x']=x
                    y=(np.divide(cross_pdt[1],v_mag))*(-1)
                    vec1_2['y']=y
                    z=(np.divide(cross_pdt[2],v_mag))*(-1)
                    vec1_2['z']=z
                    
                    vec11_2.append(vec1_1)
                    vec11_2.append(vec1_2)
                    yield vec11_2

