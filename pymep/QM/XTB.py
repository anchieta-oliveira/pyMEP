#Modulos


#############################



class XTB:
    def __init__(self,
                path:str = "mopac ",
                name:str = ""
                ) -> None:
        self.name = name
        self.outputname:str
        self.conf:str
        self.file_conf:str
        self.dir:str
        self.args_config:str 
        self.path:str = path