#Modulos


#############################



class MOPAC:
    def __init__(self,
                path:str = "mopac ",
                name:str = "",
                ) -> None:
        self.name = name
        self.outputname:str
        self.conf:str
        self.mop_file:str
        self.dir:str
        self.args_config:str 
        self.args_config_run:str 
        self.path:str = path

    def config(self, 
               coordinates:str,
               method:str = "PM7",
               charge:int = 0,
               scf:str = "1SCF",
               qmmm:bool = False,
               aux:bool = False,
               geo_ok:bool = False,
               esp: bool = False,
               graphf: bool= False,
               pdbout:bool = False,
               opt:bool = False, 
               mozyme:bool = False,
               ) -> None:


        local_vars = locals()
        args = {arg: local_vars[arg] for arg in self.config.__annotations__.keys() if arg in local_vars}
        self.args_config = args
        return args


    def make_input(self) -> str:
        conf = self.args_config
        #conf_run = self.args_config_run
        mop = ""

        # Key  Method 
        mop += f"{conf['method']} "
        # key QMMMM
        if conf['qmmm']:
            mop += "QMMM "

        # key Charge
        if conf['charge']:
            mop += f"CHARGE={conf['charge']} "

        # Key SCF    
        if conf['scf']:
            mop += f"{conf['scf']} "

        # Key MOZYME
        if conf['mozyme']:
            mop += "MOZYME "
        # Key ESP
        if conf['esp']:
            mop += "ESP "
        # Key GEO-OK
        if conf['geo_ok']:
            mop += "GEO-OK "
        # Key PDBOUT 
        if conf['pdbout']:
            mop += "PDBOUT "
        # Key graphf
        if conf['graphf']:
            mop += "GRAPHF "
        # Key AUX
        if conf['aux']:
            mop += "AUX "

        mop += "\n\n\n"    
        mop += conf['coordinates']
            

        # FIMMMMMM
        self.mop_file = mop
        return mop


    def config_run(self,
                   cores:int,
                   maxRam:int,
                   
                   ) -> None:
        local_vars = locals()
        args = {arg: local_vars[arg] for arg in self.config_run.__annotations__.keys() if arg in local_vars}
        self.args_config_run = args

        return args

    def make_points_charge_file(self, points_charge:list, atoms_qm:list) -> str:
        # http://openmopac.net/manual/QMMM.html
        pcs = f"\n{len(atoms_qm)} 0\n"

        def calculate_potential(atom) -> float:
            p =  332 * sum(map(lambda pc:pc.charge / pc.coordinates.measure_distance(atom.coordinates.x, atom.coordinates.y, atom.coordinates.z), points_charge))
            return p 

        pots = list(map(calculate_potential, atoms_qm))
        for n, i in enumerate(atoms_qm):
            pcs += f"{i.element.symbol}\t{i.coordinates.x}\t{i.coordinates.y}\t{i.coordinates.z}\t{pots[n]}\n"

        with open('mol.in', "w") as file_pdb:
            file_pdb.write(pcs)

        return "mol.in"


    def write_mop(self, name:str) -> str:
        with open(name + '.mop', "w") as file:
            file.write(self.mop_file)
        return name + ".mop"