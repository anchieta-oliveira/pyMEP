#Modulos


#############################



class ORCA:
    def __init__(self,
                path:str = "orca ",
                name:str = ""
                ) -> None:
        self.name = name
        self.outputname:str
        self.conf:str
        self.inp_file:str
        self.dir:str
        self.args_config:str 
        self.args_config_run:str 
        self.path:str = path


    def config(self, 
               coordinates:str,
               method:str = "B3LYP",
               basis:str = "6-31G*",
               basisLow:str = "",
               charge:int = 0,
               mult:str = "1.00",
               units:str = "Angs",
               points_charge:list = [],
               opt:bool = False,
               oniom:bool = False,
               high_oniom:str = "B3LYP",
               low_oniom:str = "PM3",
               AutoFF_QM2_Method:str = "XTB",
               Charge_Method:str = "Hirshfeld",
               Scale_QM2Charges_MMAtom:str = "1",
               convergence: str ="TightSCF",
               AutoTRAHIter:int = 50,
               MaxIter: int = 2000,
               DIISMaxEq:int = 15,
               AutoTRAHNInter: int = 10,
               convergence_strategies:str = "NormalConv",
               output:str = "None", 
               p_overlap:bool = False,
               p_mos:bool = False,
               p_basis: bool = False,
               p_mulliken:bool = False,
               p_atcharges_m:bool = False,
               atoms_high:list = [],
               atoms_low:list = [],
               qm_xtb:bool = False,
               xtb_keywords: str = "--iterations 500",
               tddft:bool = False,
               mdci:bool = False,
               nroots:int = 30,
               donto:bool = False,
               freq:bool = False,
               freq_temp:list = [77, 298, 330, 450],
               solvent:str = "None",
               cpcm:bool = False,
               cpcm_custom:str = "EPSILON 80.4\n REFRAC 1.33",
               smd:bool = False,
               neb_ts:bool = False,
               keepdens:bool = False,
               previous_name:str = "",
               cores: int = 1,
               maxRam: int = 4000,
               directresetfreq:int = 10,
               extra_keys:str = ""
               ) -> None:


        local_vars = locals()
        args = {arg: local_vars[arg] for arg in self.config.__annotations__.keys() if arg in local_vars}
        self.args_config = args
        return args


    def make_input(self) -> str:
        conf = self.args_config
        #conf_run = self.args_config_run
        inp = ""
        # Cálculos ONIMO internos do ORCA
        if conf['oniom'] or conf['qm_xtb']:
            if conf['qm_xtb']:
                inp += "!QM/XTB "
                inp += f"{conf['high_oniom']} "
                inp += f"{conf['basis']} \n"
                inp += f"%xtb\n\tXTBINPUTSTRING \"{conf['xtb_keywords']}\"\n"
                inp += "\nend\n"
                inp += "\n%QMMM\n"
            # Adicionando coisas do QMMMM ONIOM
            else:
                inp += "!QM/QM2 "
                inp += f"{conf['high_oniom']} "
                inp += f"{conf['basis']} "
                inp += "\n%QMMM\n"
                inp += f"\tQM2CUSTOMMETHOD \"{conf['low_oniom']}\"\n"
                inp += f"\tQM2CUSTOMBASIS \"{conf['basisLow']}\"\n"
                inp += f"\tAutoFF_QM2_Method {conf['AutoFF_QM2_Method']}\n"
                inp += f"\tCharge_Method {conf['Charge_Method']}\n"
                inp += f"\tScale_QM2Charges_MMAtom {conf['Scale_QM2Charges_MMAtom']}\n"

            inp += "\tQMATOMS {" + str(conf['atoms_high']).replace(",", "")[1:-1] + "} end"
            if len(conf['atoms_low']) > 0 and not len(conf['atoms_high']) > 0:
                inp += "\tQM2ATOMS {" + str(conf['atoms_low']).replace(",", "") + "} end"
            
            inp += "\nEND\n"
        else:
            #Add. mtd de basis se não for ONIOM
            inp += f"!{conf['method']} {conf['basis']}\n"
        
        inp += "!" + conf['convergence'] + "\n"

        if conf['points_charge']:
            pcs_file = self.make_points_charge_file(points_charge=conf['points_charge'])
            inp += f"%pointcharges \"./{pcs_file}\""

        if conf['opt']: 
            inp += "!OPT\n"

        if conf['output'] != "None" or bool(conf['p_mos']) or bool(conf['p_overlap']) or bool(conf['p_basis']) or bool(conf['p_mulliken']) or bool(conf['p_atcharges_m']):
            if conf['output'] != "None":
                inp += f"\n%output\n{conf['output']}\n"
            else:
                inp += f"%output\n"

            if bool(conf['p_mos']):
               inp += "print[p_mos] 1\n"

            if bool(conf['p_overlap']):
                inp += "print[p_overlap] 5\n"
            
            if bool(conf['p_atcharges_m']):
                inp += "Print[p_atcharges_m] 1\n"
            
            if bool(conf['p_basis']):
                inp += "print[p_basis] 2\n"

            if bool(conf['p_mulliken']):
                inp += "print[p_mulliken] 1\n"
            
            inp += "end\n"


        inp += f"\n\n%scf\nAutoTRAHIter {conf['AutoTRAHIter']}\nMaxIter {conf['MaxIter']}\nAutoTRAHNInter {conf['AutoTRAHNInter']}\nDIISMaxEq {conf['DIISMaxEq']}\ndirectresetfreq {conf['directresetfreq']}\nend\n"
        inp += "\n#Keywords run\n"
        inp += f"%MAXCORE {conf['maxRam']}\n"
        inp += f"%PAL NPROCS {conf['cores']} END\n"
        inp += f"!{conf['convergence_strategies']}\n"
        #Coisas extras 
        inp += "\n\n#Keywords extras\n"
        inp += conf['extra_keys']
    
        # Corrds vem por último para facilitar a leitura 
        inp += "\n#Coords JOB\n"
        inp += f"%coords\n\tCTyp xyz\n\tCharge {conf['charge']}\n\tMult {conf['mult']}\n\tUnits {conf['units']}\n\tcoords\n"
        inp += conf['coordinates']
        inp += "\nEND\n"
        

        inp += "\nEND" # End finallll

        # FIMMMMMM
        self.inp_file = inp
        return inp


    def config_run(self,
                   cores:int,
                   maxRam:int,
                   
                   ) -> None:
        local_vars = locals()
        args = {arg: local_vars[arg] for arg in self.config_run.__annotations__.keys() if arg in local_vars}
        self.args_config_run = args

        return args

    def make_points_charge_file(self, points_charge:list) -> str:
        pcs = f"{len(points_charge)}\n"
        for point in points_charge:
            pcs += f"{point.charge}\t{point.coordinates.x}\t{point.coordinates.y}\t{point.coordinates.z}\n"
        
        with open('pointscharge', "w") as file_pdb:
            file_pdb.write(pcs)

        return "pointscharge"


    def write_inp(self, name:str) -> str:
        with open(name + '.inp', "w") as file_pdb:
            file_pdb.write(self.inp_file)
        return name + ".inp"




