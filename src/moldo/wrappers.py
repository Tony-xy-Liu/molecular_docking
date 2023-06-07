import os
from pathlib import Path
from vina import Vina as _vina
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from rdkit.Geometry import Point3D
import pandas as pd
import numpy as np
import re

Vector3 = tuple[float, float, float]

_RAND_SEED = 42

def SetSeed(n: int):
    global _RAND_SEED
    _RAND_SEED = n

def _get_file_name(file_path: Path, no_periods=False):
    if no_periods: # only return first segment between periods "a.b.c" -> "a"
        return str(file_path.name).split('.')[0]
    else: # remove file extension only "a.b.c" -> "a.b"
        return ".".join(str(file_path.name).split('.')[:-1])

def _ensure_folder_exists(folder: Path):
    if not os.path.isdir(folder):
        os.mkdir(folder)

def PrepareReceptor(receptor_path: Path, cache: Path):
    assert receptor_path.name.endswith('.pdb'), "pdb format required"
    _ensure_folder_exists(cache)

    fname = _get_file_name(receptor_path)
    prepped_path = cache.joinpath(f"{fname}.raw.pdbqt")
    os.system(f'prepare_receptor -r {receptor_path} -o {prepped_path} -A hydrogens')

    out_pdbqt = cache.joinpath(f"{fname}.pdbqt")
    with open(prepped_path) as pdbqt:
        with open(out_pdbqt, 'w') as out:
            for l in pdbqt:
                l_len = len(l)
                if l_len == 81:
                    l = "".join((l[:21], l[22:]))
                out.write(l)

    out_pdb = cache.joinpath(f"{fname}.prepped.pdb")
    os.system(f"obabel {out_pdbqt} -opdb > {out_pdb} && rm {prepped_path}")
    # os.system(f'singularity run -B {CACHE}:/ws {_LIB}/obabel obabel /ws/{out_path.split("/")[-1]} -opdb > {readable}')
    return out_pdbqt, out_pdb

def PrepareLigand(file_path: Path, active_sites: list[Vector3], cache: Path):
    assert file_path.name.endswith('.sdf'), "sdf format required"
    _ensure_folder_exists(cache)

    lig = Chem.MolFromMolFile(f"{file_path}")
    lig = Chem.AddHs(lig)
    AllChem.EmbedMolecule(lig,randomSeed=_RAND_SEED)  # 3d projection
    lig_conf = lig.GetConformer()

    lig_atom_pos = []
    for i in range(lig_conf.GetNumAtoms()):
        pos = lig_conf.GetAtomPosition(i)
        lig_atom_pos.append((pos.x, pos.y, pos.z))

    to_tup = lambda c: (c.x, c.y, c.z)
    lig_center = to_tup(ComputeCentroid(lig_conf))
    original_pos = [(pos.x, pos.y, pos.z) for pos in [lig_conf.GetAtomPosition(i) for i in range(lig_conf.GetNumAtoms())]]

    out_paths: list[Path] = []
    for j, active_site in enumerate(active_sites):
        for i in range(lig_conf.GetNumAtoms()):
            pos = original_pos[i]
            x, y, z = [v-c+d for v, c, d in zip(pos, lig_center, active_site)] # center atom, then to active site
            lig_conf.SetAtomPosition(i,Point3D(x,y,z))

        lig_name = _get_file_name(file_path)
        lig_th_path = cache.joinpath(f"{lig_name}_AS{j+1:02}.hy.sdf")
        w = Chem.SDWriter(f"{lig_th_path}")
        w.write(lig)
        w.flush()
        w.close()

        out_path = Path(f'{lig_th_path}'.replace(".sdf", ".pdbqt"))
        # # os.system(f'{_ADFR}/pythonsh {_LIB}/mk_prepare_ligand.py -i {lig_th_path} -o {out_path}')
        # os.system(f'singularity run -B {CACHE}:/ws {_LIB}/obabel obabel /ws/{lig_th_path.split("/")[-1]} -opdbqt > {out_path}')
        os.system(f'''\
            obabel {lig_th_path} -opdbqt > {out_path}
            rm {lig_th_path}
        ''')
        out_paths.append(out_path)
    return out_paths

def MakeFlexReceptor(receptor_pdb: Path, receptor_pdbqt: Path, cache: Path, active_site: Vector3, num_residues: int=10):
    _ensure_folder_exists(cache)

    residues = []
    with open(receptor_pdb) as recf:
        for l in recf:
            if not l.startswith('ATOM'): continue
            coords = [float(t) for t in l[31:54].split(' ') if t != '']
            atom = l[22:26].strip()
            aa = l[17:20]
            residues.append((aa, atom, coords))
    coords = np.array([c for _, _, c in residues])
    dists = np.abs(np.linalg.norm(active_site - coords, axis=1))
    dists = [(i, d) for i, d in enumerate(dists)]
    dists = sorted(dists, key=lambda t: t[1])

    to_make_flexible = set()
    for i, d in dists:
        if len(to_make_flexible) >= num_residues: break
        res = "".join(residues[i][:2])
        if res in to_make_flexible: continue
        to_make_flexible.add(res)
    
    with open(cache.joinpath("residues_made_flexible.txt"), "w") as f:
        for res in to_make_flexible:
            f.write(f"{res}\n") 

    flex, rigid = cache.joinpath('flex.pdbqt'), cache.joinpath('rigid.pdbqt')
    cmd = f"""\
        pythonsh /opt/adfr/prepare_flexreceptor.py \
        -r {receptor_pdbqt} -s {'_'.join(to_make_flexible)} \
        -x {flex} -g {rigid}
    """
    # print(cmd)
    os.system(cmd)
    return flex, rigid

def PredictActiveSites(receptor_path: Path, cache: Path, top_k=1) -> list[Vector3]:
    assert receptor_path.name.endswith('.pdb'), "pdb format required"

    fname = _get_file_name(receptor_path, no_periods=True)
    csv_dir = cache.joinpath(f'active_site_predictions__{fname}')
    os.system(f'''\
        prank predict -f {receptor_path} -o {csv_dir}
        cd {csv_dir}
        mv *predictions.csv predictions.csv
        mv *residues.csv residues.csv
    ''')

    pock_df = pd.read_csv(f'{csv_dir}/predictions.csv')
    active_sites = []
    count = 0
    for i, row in pock_df.iterrows():
        loc = tuple(float(str(x)) for x in (row['   center_x'], row['   center_y'], row['   center_z']))
        active_sites.append(loc)
        count += 1
        if count >= top_k: break
    return active_sites

class Vina:
    def __init__(self, receptor_path: Path, active_site:Vector3, flex_path=None, cpus:int=0, box_len:float=20) -> None:
        v = _vina(sf_name='vina', cpu=cpus, seed=_RAND_SEED, verbosity=1)
        self._vina = v
        self._receptor_name = _get_file_name(receptor_path, no_periods=True)
        self._receptor_path = receptor_path
        self._flex_path = flex_path

        if flex_path is None:
            v.set_receptor(f"{receptor_path}")
        else:
            v.set_receptor(flex_pdbqt_filename=f"{flex_path}", rigid_pdbqt_filename=f"{receptor_path}")

        box = [box_len]*3
        print(f'computing vina grid for {self._receptor_name}')
        v.compute_vina_maps(center=active_site, box_size=box)

    def Dock(self, ligand_path: Path, cache: Path, exhaustiveness: int=16, n_poses=3):
        lig_name = _get_file_name(ligand_path, no_periods=True)
        # out_name = f'{self._receptor_name}__{lig_name}'
        out_dir = cache
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        v = self._vina
        v.set_ligand_from_file(f"{ligand_path}")
        print(f'docking {lig_name}')
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        energies = v.energies(n_poses=n_poses)[:, 0]
        vina_out = out_dir.joinpath(f'vina_out.pdbqt')
        v.write_poses(f"{vina_out}", overwrite=True, n_poses=n_poses, energy_range=9999999)
        # _results = []
        # with open(vina_out) as f:
        #     for l in f:
        #         if "VINA RESULT" not in l: continue
        #         _results.append([float(v) for v in re.findall(r'-?\d+\.\d+', l)])

        # results = np.array(_results, dtype=float)
        # energies = results[:, 0]

        models = []
        current_model = []
        def submit():
            nonlocal current_model
            models.append(current_model)
            current_model = []

        with open(vina_out) as pdbqt:
            for l in pdbqt:
                if l.startswith('MODEL'): submit()
                current_model.append(l)
            submit()

        def run_if_safe(fn):
            try: fn()
            except FileExistsError: pass
        run_if_safe(lambda: os.mkdir(out_dir))
        run_if_safe(lambda: os.mkdir(f'{out_dir}/pdbqt'))
        run_if_safe(lambda: os.mkdir(f'{out_dir}/sdf'))

        print(f'separating {len(models)-1} poses')
        poses:list[str] = []
        for i, model in enumerate(models[1:]): # first will be blank
            i += 1
            with open(f'{out_dir}/pdbqt/{i}.pdbqt', 'w') as out:
                out.writelines(model[:31])
            pose_path = f'{out_dir}/sdf/{lig_name}_{i:02}.sdf'
            os.system(f'obabel {out_dir}/pdbqt/{i}.pdbqt -osdf > {pose_path}')
            poses.append(pose_path)
        os.system(f"rm -r {out_dir}/pdbqt")

        return vina_out, poses, energies
