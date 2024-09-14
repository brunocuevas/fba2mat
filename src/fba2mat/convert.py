from typing import List, Dict
import pandas as pd
import numpy as np


def index_species(species):
    m = pd.DataFrame.from_records(species).set_index("id")
    return m


def read_reactions(reactions: List[Dict], index: pd.DataFrame):
    # Part 1: Creating new index, to remove unused molecules
    new_index = []
    for i, reaction in enumerate(reactions):
        for key in reaction["metabolites"].keys():
            species = index.loc[key].to_dict()
            species["id"] = key
            new_index.append(species)
    new_index = (
        pd.DataFrame.from_records(new_index)
        .drop_duplicates(subset=["id"])
        .set_index("id")
    )
    new_index["position"] = np.arange(len(new_index))
    # Part 2: Build stoichiometric matrix
    if len(reactions) * len(index) > 10000:
        print("Caution: the resulting stoichiometric matrix is going to be very big")
    M = np.zeros((len(reactions), len(new_index)))
    for i, reaction in enumerate(reactions):
        for key, value in reaction["metabolites"].items():
            j = new_index.loc[key]["position"]
            M[i, j] = value
    return M, new_index.reset_index()


def classify_reaction(reaction: Dict, index: pd.DataFrame):
    stoichiometry = reaction["metabolites"]
    k_substrates = [k for k, s in stoichiometry.items() if s < 0]
    k_products = [k for k, s in stoichiometry.items() if s > 0]
    if len(k_substrates) == 0:
        return "export"
    elif len(k_products) == 0:
        return "import"
    else:
        names_in = sorted([index.loc[k]["name"] for k in k_substrates])
        names_out = sorted([index.loc[k]["name"] for k in k_products])
        check_in = all([i in names_in for i in names_out])
        check_out = all([i in names_out for i in names_in])
        if check_in and check_out:
            return "transport"
        else:
            return "metabolic"


def generate_string(metabolites):
    substrates = " + ".join(
        [f"{int(-s):2d}:{k:s}" for k, s in metabolites.items() if s < 0]
    )
    products = " + ".join([f"{int(s):2d}:{k:s}" for k, s in metabolites.items() if s > 0])
    return substrates + "->" + products
