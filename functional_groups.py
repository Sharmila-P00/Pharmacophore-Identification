

def groups():
    functional_groups={
        "Aromatic":["a1aaaaa1","a1aaaa1"],
                    "HBD":["[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]",
                            "[#8!H0&!$([OH][C,S,P]=O)]",
                            "[#16!H0]"],
                    "HBA":["[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]",
                            "[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]"],
                    "Positive":["[+,+2,+3,+4]",
                                "[$(CC)](=N)N",
                                "[$(C(N)(N)=N)]",
                                "[$(n1cc[nH]c1)]"],
                                
                    "Negative":["[-,-2,-3,-4]",
                                "C(=O)[O-,OH,OX1]",
                                "[$([S,P](=O)[O-,OH,OX1])]",
                                "c1[nH1]nnn1",
                                "c1nn[nH1]n1",
                                "C(=O)N[OH1,O-,OX1]",
                                " C(=O)N[OH1,O-]",
                                "CO(=N[OH1,O-])",
                                "[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]"],
                                
                    "Hydrophobic":["a1aaaaa1",
                                    "a1aaaa1",
                                    "[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]",
                                    "[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
                                    "*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
                                    "[C&r3]1~[C&r3]~[C&r3]1",
                                    "[C&r4]1~[C&r4]~[C&r4]~[C&r4]1",
                                    "[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1",
                                    "[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1",
                                    "[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1",
                                    "[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1",
                                    "[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]",
                                    "[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
                                    "[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]",
                                    "[$([S]~[#6])&!$(S~[!#6])]"] }
    

    return(functional_groups)