import os, sys
import argparse
from pathlib import Path
import pandas as pd
from .wrappers import MakeFlexReceptor, PrepareReceptor, PredictActiveSites, PrepareLigand, Vina, SetSeed

HERE = Path("/".join(os.path.realpath(__file__).split('/')[:-1]))
CMD = "moldo"
NAME = "simple molecular docking with Autodock Vina & P2Rank"

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\n%s: error: %s\n' % (self.prog, message))

def fprint(x):
    return print(f"{x}".replace("  ", ""))

_ver = None
def _get_ver():
    global _ver
    if _ver is None:
        with open(HERE.joinpath("version.txt")) as f:
            _ver = f.read()
    return _ver

_head = None
def _get_header():
    global _head
    if _head is None:
        _head = f"""\
        {NAME}
        v{_get_ver()}
        """
    return _head

def print_header():
    fprint("#"*os.get_terminal_size().columns)
    fprint(_get_header)

def main():
    parser = ArgumentParser(prog=f'{CMD}')
    parser.add_argument("--receptor", "-r", metavar="PATH", required=True,  help=".pdb file of receptor")
    parser.add_argument("--ligand", "-l",   metavar="PATH", required=True,  help=".sdf file of ligand")
    parser.add_argument("--out", "-o",      metavar="PATH", default=os.getcwd(),   help="output folder")
    for param, short, default, help in [
        ("active-sites", "a", 3, "number of active sites to try"),
        ("flexes", "f", 6, "number of residues near the active site to simulate dynamically"),
        ("attempts", "n", 100, "number of docking attempts to try"),
        ("box-len", "b", 20, "simulation box's side length in angstroms"),
        ("exhaustiveness", "x", 16, "how hard autodock should try for each attempt"),
        ("threads", "t", 8, "number of cpu threads to utilize"),
        ("seed", "s", 42, "seed for random number generator"),
    ]:
        parser.add_argument(f"--{param}", f"-{short}", metavar="INT", default=default, help=help)

    # parse paths
    def check_path(p, ext_wl: list[str]=list(), is_file=True):
        if is_file and not os.path.exists(p):
            print(f"path [{p}] does not exist"); return False
        if not is_file and os.path.exists(p):
            print(f"output folder [{p}] already exists"); return False

        p = Path(p)
        if is_file:
            if not p.is_file():
                print(f"expected file, path [{p}] is not a file"); return False
            if len(ext_wl)>0 and not any(str(p).endswith(ext) for ext in ext_wl):
                print(f"read file [{p}] doesn't end in one of {ext_wl}"); return False
        else:
            if p.is_file():
                print(f"expected folder, path [{p}] is a file"); return False
            if not p.exists(): os.makedirs(p, exist_ok=True)
        return p
    
    args = parser.parse_args(sys.argv[1:])
    _paths: list[Path] = []
    for p, ext, is_f in [
        (args.receptor, ".pdb", True),
        (args.ligand, ".sdf", True),
        (args.out, None, False),
    ]:
        if not check_path(p, [ext], is_f): return
        _paths.append(Path(p))
    receptor, ligand, out = _paths

    # parse nums
    _nums: list[int] = []
    for k in "active_sites, attempts, exhaustiveness, flexes, box_len, threads, seed".split(", "):
        n = args.__dict__[k]
        try:
            _nums.append(int(n))
        except:
            print(f"expected int for {k}, got [{n}]")
            return
    n_active_sites, n_poses, n_ex, n_flexes, box_len, n_threads, rng_seed = _nums
    SetSeed(rng_seed)

    # dock
    work = out.joinpath("docking_work")
    receptor_pdbqt, receptor_pdb = PrepareReceptor(receptor, work)
    active_sites = PredictActiveSites(receptor_pdb, work, top_k=n_active_sites)
    positioned_ligands = PrepareLigand(ligand, active_sites, work)
    results = []
    for i, (asite, lig) in enumerate(zip(active_sites, positioned_ligands)):
        acache = work.joinpath(f"flexible_receptor_for_active_site_{i+1:02}")

        if n_flexes <=0:
            v = Vina(receptor_pdbqt, active_site=asite, cpus=n_threads, box_len=box_len)
        else:
            flex, rigid = MakeFlexReceptor(receptor_pdb, receptor_pdbqt, acache, asite, num_residues=n_flexes)
            v = Vina(rigid, flex_path=flex, active_site=asite, cpus=n_threads, box_len=box_len)

        docking_out = work.joinpath(f"results_for_active_site_{i+1:02}")
        _rec, _ligs, energies = v.Dock(ligand_path=lig, exhaustiveness=n_ex, n_poses=n_poses, cache=docking_out)
        
        # print(energies, [type(e) for e in energies])
        x, y, z = asite
        for docked_ligand, e in zip(_ligs, energies):
            results.append((e, x, y, z, f"{docked_ligand}"))

    df_results = pd.DataFrame(results, columns="kcal_per_mol, x, y, z, ligand_path".split(", "))
    df_results.to_csv(out.joinpath("results.csv"), index=False)
