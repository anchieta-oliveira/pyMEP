import sys
import time
import numpy as np
import seaborn as sns
sys.path.append('../')
from pymep.mep import MEP
from pymep.MOL.PSF import PSF
from pymep.MOL.PDB import PDB
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def read_charges_file(charge_file) -> np.array:
    with open(charge_file, "r") as file:
        c = file.readlines()
        charges = np.array(list(map(lambda x: float(x),c)))
        return charges

def run_trp_cage():
    pdb = PDB()
    mep = MEP()

    print("Lendo PDB and PSF...")
    pdb.read("2luf_ph74_autopsf.pdb")
    psf = PSF("2luf_ph74_autopsf.psf")

    ts_cpu = []
    gps_cpu = []
    gs = 1

    while gs >= 0.1:
        print(f"MEP res: {gs}")
        start_time = time.time() 
        dx = mep.calculate(pdb=pdb, charges=psf.charges, res=gs, gpu=False, margim=0.3, cutoff=30, gpus_id=[0], form="dx")
        end_time = time.time()
        t = end_time - start_time
        gp = dx.xn * dx.yn * dx.zn

        ts_cpu.append(t)
        gps_cpu.append(gp)
        print(f"Grid spacing {t} s - {gp} grid points")
        gs -= 0.1

    ts_gpu = []
    gps_gpu = []
    gs = 1
    while gs >= 0.1:
        print(f"MEP res: {gs}")
        start_time = time.time() 
        dx = mep.calculate(pdb=pdb, charges=psf.charges, res=gs, gpu=True, margim=0.3, cutoff=30, gpus_id=[0], form="dx")
        end_time = time.time()
        t = end_time - start_time
        gp = dx.xn * dx.yn * dx.zn

        ts_gpu.append(t)
        gps_gpu.append(gp)
        print(f"Grid spacing {t} s - {gp} grid points")
        gs -= 0.1

    
    plt.figure(figsize=(10,6))
    sns.set(style="whitegrid")

    sns.lineplot(x=[gp / 1e6 for gp in gps_cpu], y=ts_cpu, marker='o', color='b', label='CPU - Execution Time')

    sns.lineplot(x=[gp / 1e6 for gp in gps_gpu], y=ts_gpu, marker='o', color='r', label='GPU - Execution Time')

    plt.title('Benchmark: Execution Time vs. Number of Grid Points (CPU vs GPU)')
    plt.xlabel('Number of Grid Points (millions)')
    plt.ylabel('Execution Time (seconds)')

    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:.1f}M'))

    plt.legend()
    plt.savefig("1o0h_cpu_gpu.png", dpi=250)
    plt.show()


def run_1o0h():
    pdb = PDB()
    mep = MEP()

    print("Lendo PDB and Charges...")
    pdb.read("1o0h_protein.pdb")
    charges = read_charges_file("1o0h_protein_xtb.charges")

    ts_cpu = []
    gps_cpu = []
    gs = 1

    while gs >= 0.3:
        print(f"MEP gs: {gs}")
        start_time = time.time() 
        dx = mep.calculate(pdb=pdb, charges=charges, res=gs, gpu=False, margim=0.3, cutoff=30, gpus_id=[0], form="dx")
        end_time = time.time()
        t = end_time - start_time
        gp = dx.xn * dx.yn * dx.zn

        ts_cpu.append(t)
        gps_cpu.append(gp)
        print(f"Grid spacing {t} s - {gp} grid points")
        gs -= 0.1

    ts_gpu = []
    gps_gpu = []
    gs = 1
    while gs >= 0.3:
        print(f"MEP gs: {gs}")
        start_time = time.time() 
        dx = mep.calculate(pdb=pdb, charges=charges, res=gs, gpu=True, margim=0.3, cutoff=30, gpus_id=[0], form="dx")
        end_time = time.time()
        t = end_time - start_time
        gp = dx.xn * dx.yn * dx.zn

        ts_gpu.append(t)
        gps_gpu.append(gp)
        print(f"Grid spacing {t} s - {gp} grid points")
        gs -= 0.1

    dx.write(f"1o0h_gs{gs}_gpu_mep.dx")
    
    plt.figure(figsize=(10,6))
    sns.set(style="whitegrid")

    sns.lineplot(x=[gp / 1e6 for gp in gps_cpu], y=ts_cpu, marker='o', color='b', label='CPU - Execution Time')

    sns.lineplot(x=[gp / 1e6 for gp in gps_gpu], y=ts_gpu, marker='o', color='r', label='GPU - Execution Time')

    plt.title('Benchmark: Execution Time vs. Number of Grid Points (CPU vs GPU)')
    plt.xlabel('Number of Grid Points (millions)')
    plt.ylabel('Execution Time (seconds)')

    plt.gca().xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f'{x:.1f}M'))

    plt.legend()
    plt.savefig("1o0h_cpu_gpu.png", dpi=250)
    plt.show()


# Run...
run_trp_cage()
run_1o0h()