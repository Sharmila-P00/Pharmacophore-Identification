from features import parameters
from functional_groups import groups


Param=parameters()
features=Param.descriptor(file,groups())

class radius:
    def radii(self,features):
        radii={}
        radii_map = {
            "Aro": 1.1,
            "HBD": 0.5,
            "HBA": 0.5,
            "Negative": 0.75,
            "Positive": 0.75,
            "Hyd": 1.0
        }
        for i in features.keys():
            j=i.strip()
            key = j[:3] if j[:3] in radii_map else j
            if key in radii_map:
                radii[j] = radii_map[key]   
        return radii
    


