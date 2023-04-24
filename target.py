import itertools
import os
import glob

from metrics import GDT
from metrics import INF
from metrics import Clashscores
from metrics import TMScore
from metrics import LDDT

from utils import preprocess_rec
from utils import is_valid_pdb

def consolidate_target(target_path, metric="all"):
    """
    Consolidate all metrics for a target.

    Parameters:
    target_path (str): Path to directory containing PDB files.

    Returns:
    None
    """
    
    target_name = os.path.basename(target_path)

    gdt = GDT()
    inf = INF()
    clashscore = Clashscores()
    tm_score = TMScore()
    lddt = LDDT()

    if metric == "all":
        gdt.consolidate(target_name)
        clashscore.consolidate(target_name)
        lddt.consolidate(target_name)
        inf.consol_bulk(target_name)
        tm_score.consol_bulk(target_name)
        
        df = gdt.consol_gdt_percent(target_name)
        gdt.save_consol_gdt_percent_fig(target_name, df)
    elif metric == "gdt":
        gdt.consolidate(target_name)
        
        df = gdt.consol_gdt_percent(target_name)
        gdt.save_consol_gdt_percent_fig(target_name, df)
    elif metric == "clashscore" or metric == "clash":
        clashscore.consolidate(target_name)
    elif metric == "lddt":
        lddt.consolidate(target_name)
    elif metric == "tm_score":
        tm_score.consol_bulk(target_name)
    elif metric == "inf":
        inf.consol_bulk(target_name)

def calculate_target(target_path, metric="all", preprocess=True, force=False, run_summary=False):
    """
    Calculate all metrics for a target.

    Parameters:
    target_path (str): Path to directory containing PDB files.
    metric (str): Which metric to calculate. Options are "all", "gdt", "clashscore", "lddt", "tm_score", "inf".
    preprocess (bool): Whether to preprocess the PDB files.
    force (bool): Whether to force recalculation of metrics.
    run_summary (bool): Whether to only consolidate metrics and skip metrics calculation

    Returns:
    None
    """
    assert os.path.isdir(target_path), f"{target_path} does not exist."
    if preprocess: 
        preprocess_rec(target_path)

    gdt = GDT()
    inf = INF()
    clashscore = Clashscores()
    tm_score = TMScore()
    lddt = LDDT()

    # print(target_path)
    r_glob = f"{target_path}/references/processed/*.pdb"
    m_glob = f"{target_path}/models/processed/*.pdb"
    
    # Loop through each file that matches t_glob or m_glob
    
    glob_reference = glob.glob(r_glob)
    glob_models = glob.glob(m_glob)

    target_name = os.path.basename(target_path)

    for reference in glob_reference:
        if metric == "inf":
            print("Running (bulk) inf.")
            inf.calc_bulk(reference, glob_models, force=force)
            
        if metric == "tm_score":
            print("Running (bulk) tm_score.")
            tm_score.calc_bulk(reference, glob_models, force=force)

    for reference, model in itertools.product(glob_reference, glob_models):
        print(f"reference: {reference}")
        print(f"model: {model}")
        assert is_valid_pdb(reference) and is_valid_pdb(model), "reference and model PDB must be valid."

        print(f"Running {metric} on reference {reference} and model {model}")

        if metric == "all":
            print("Running all methods.")
            gdt.calculate(reference, model, force=force)
            clashscore.calculate(reference, model, force=force)
            lddt.calculate(reference, model, force=force)

        elif metric == "gdt":
            print("Running gdt.")
            gdt.calculate(reference, model, force=force)
        elif metric == "clashscore" or metric == "clash":
            print("Running clashscore.")
            clashscore.calculate(reference, model, force=force)
        elif metric == "lddt":
            print("Running lddt.")
            lddt.calculate(reference, model, force=force)
        elif metric == "inf" or metric == "tm_score":
            return
        else:
            raise ValueError(f"{metric} is not a valid metric.")
            
    if (run_summary): consolidate_target(target_path, metric=metric)
    
