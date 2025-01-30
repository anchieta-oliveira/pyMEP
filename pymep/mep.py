import os
import numba
import numpy as np
from tqdm import tqdm
from pymep.VOL.dx import DX
from pymep.QM.aux import AUX
from pymep.MOL.PDB import PDB
from pymep.VOL.cube import Cube
from pymep.QM.orca_out import OrcaOut
from pymep.FF.forceField import ForceField

class MEP:
    def __init__(self) -> None:
        pass

    
    def calculate(self, pdb:PDB, aux:AUX=AUX(), orca_out:OrcaOut=OrcaOut(), res:float=0.5, gpu:bool=False, charges:np.array=np.array([]), FF:str="", form:str="cube", margim:float=0.4, cutoff: float = 15):
        if gpu:
            import cupy as cp
            
        #matrix_dist = pdb.get_distance_matrix()
        #dist_max = matrix_dist.max()
        coords_atoms = pdb.coordinates

        x_min, y_min, z_min = coords_atoms.min(axis=0)
        x_max, y_max, z_max = coords_atoms.max(axis=0)

        dist_max = np.linalg.norm([x_max - x_min, y_max - y_min, z_max - z_min])

        xorg, yorg, zorg = coords_atoms[:, 0].min() - (dist_max * margim), coords_atoms[:, 1].min() - (dist_max * margim), coords_atoms[:, 2].min()-(dist_max * margim)
        xmax, ymax, zmax = coords_atoms[:, 0].max() + (dist_max * margim), coords_atoms[:, 1].max() + (dist_max * margim), coords_atoms[:, 2].max()+ (dist_max * margim)

        xn = int((xmax - xorg) / res) + 1
        yn = int((ymax - yorg) / res) + 1
        zn = int((zmax - zorg) / res) + 1
        d = DX()
        d.make(xdel=res , ydel=res , zdel=res,
               xn=xn, yn=yn, zn=zn,
               xorg=xorg, yorg=yorg, zorg=zorg,
               )
        x = d.coordinates[:, 0]
        y = d.coordinates[:, 1]
        z = d.coordinates[:, 2]
            
        if len(aux.atom_charges) > 0:
            zatoms = np.array(aux.atom_charges)

        elif len(orca_out.mulliken_atomic_charges) > 0:
            zatoms = orca_out.mulliken_atomic_charges
        elif FF == "charmm":
            zatoms = np.zeros([len(pdb.atoms)],dtype=np.float32)
            dir_app = os.path.dirname(os.path.realpath(__file__))
            ff = ForceField(path_ff=f"{dir_app}/FF/top_all36_prot.rtf")
            
            for i, at in enumerate(pdb.atoms):
                zatoms[i] = ff.get_atom_charge(at.resname, at.name)
            print(zatoms)
            zatoms = np.nan_to_num(zatoms, nan=.0)
        elif charges.size > 0:
            zatoms = charges

        if gpu:
            import cupy as cp
            x = cp.asarray(x, dtype=cp.float32)
            y = cp.asarray(y, dtype=cp.float32)
            z = cp.asarray(z, dtype=cp.float32)
            coords_atoms = cp.asarray(coords_atoms, dtype=cp.float32)
            zatoms = cp.asarray(zatoms, dtype=cp.float32)

        if gpu:
            fun = self.comput_mep_gpu
        else:
            fun = self.comput_mep
        
        mask = zatoms != 0  
        coords_atoms = coords_atoms[mask]
        zatoms = zatoms[mask]

        gmep = fun(catoms=coords_atoms, zatoms=zatoms, x=x, y=y, z=z, xn=d.xn, yn=d.yn, zn=d.zn)

        if gpu:
            gmep = cp.asnumpy(gmep)
        d.values = gmep.flatten()

        if form == "cube":
            a_to_bohr = 1.8897259886
            c = Cube()
            c.make(xdel=res*a_to_bohr , ydel=res*a_to_bohr , zdel=res*a_to_bohr,
                xn=xn, yn=yn, zn=zn,
                xorg=xorg*a_to_bohr, yorg=yorg*a_to_bohr, zorg=zorg*a_to_bohr,
                )
            c.natoms = len(pdb.atoms)
            c.header = "MEP\n"
            for at in pdb.atoms:
                at.atomic_number
                at.coordinates.convert_to(unit="bohr")
            c.atoms = pdb.atoms
            c.values = d.values
            return c
        return d

    @staticmethod
    @numba.njit(parallel=True, cache=True, fastmath=True)
    def comput_mep(catoms: np.array, zatoms: np.array, x: np.array, y: np.array, z: np.array, xn: int, yn: int, zn: int) -> np.array:
        grid = np.stack((x, y, z), axis=-1)
        gmep = np.zeros((xn, yn, zn),  dtype=np.float32)
        
        for i in numba.prange(catoms.shape[0]):
            r = np.sqrt(np.sum((grid - catoms[i]) ** 2, axis=-1))
            r[r == 0] = np.inf
            contribution = zatoms[i] / r
            gmep += contribution.reshape((xn, yn, zn))
        return gmep

    @staticmethod    
    def comput_mep_gpu(catoms, zatoms, x, y, z, xn:int, yn:int, zn:int):
        import cupy as cp
        grid = cp.stack((x, y, z), axis=-1) 
        gmep = cp.zeros((xn, yn, zn), dtype=cp.float32)  

        for i in tqdm(range(catoms.shape[0]), desc="Processando átomos", unit="átomo"):
            r = cp.sqrt(cp.sum((grid - catoms[i]) ** 2, axis=-1))
            r = cp.where(r == 0, cp.inf, r)
            contribution = zatoms[i]/r
            gmep += contribution.reshape((xn, yn, zn))
        return gmep
