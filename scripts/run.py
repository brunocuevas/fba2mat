import click
import json
from fba2mat.convert import *
import numpy as np
import pandas as pd


@click.command()
@click.option("--filter-export", default=True)
@click.option("--filter-import", default=True)
@click.option("--filter-transport", default=True)
@click.argument("FILEIN", type=click.File("r"))
@click.argument("MATRIX_OUT", type=click.File("w"))
@click.argument("SPECIES_INDEX_OUT", type=click.File("w"))
@click.argument("REACTION_INDEX_OUT", type=click.File("w"))
def process(
    filein,
    matrix_out,
    species_index_out,
    reaction_index_out,
    filter_export,
    filter_import,
    filter_transport,
):
    """
    Converts a Flux Balance Analysis model into a matrix
    for stoichiometric / flux-balance analysis.
    """
    u = json.load(filein)

    species_index = index_species(u["metabolites"])
    reactions = u["reactions"]
    for r in reactions:
        r["type"] = classify_reaction(r, species_index)
    if filter_export:
        reactions = [r for r in reactions if r["type"] != "export"]
    if filter_import:
        reactions = [r for r in reactions if r["type"] != "import"]
    if filter_transport:
        reactions = [r for r in reactions if r["type"] != "transport"]
    print(f"-- reactions: {len(reactions)}")
    reactions_index = pd.DataFrame.from_records(reactions)

    reactions_index["reaction_string"] = reactions_index["metabolites"].apply(
        generate_string
    )
    M, new_species_index = read_reactions(reactions, species_index)

    new_species_index = new_species_index[["position", "id", "name", "compartment", "charge", "formula"]]
    reactions_index = reactions_index[["id", "name", "reaction_string"]]

    np.savetxt(matrix_out, M, fmt="%4.2f")
    new_species_index.to_csv(species_index_out, sep=";")
    reactions_index.to_csv(reaction_index_out, sep=";")


if __name__ == "__main__":
    process()
