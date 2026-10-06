"""Protein-chain inspection and extraction helpers."""

from Bio.PDB import PDBIO, PDBParser, Select


STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
}


class ProteinChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def _protein_residue_count(chain):
    return sum(
        residue.resname.strip() in STANDARD_AMINO_ACIDS
        for residue in chain.get_residues()
    )


def available_protein_chains(structure):
    """Return protein-containing chains in the first model and their residue counts."""
    model = next(structure.get_models(), None)
    if model is None:
        return {}
    return {
        chain.id: _protein_residue_count(chain)
        for chain in model
        if _protein_residue_count(chain) > 0
    }


def select_chain(structure, requested_chain_id="A"):
    """Use the requested protein chain, or the largest available protein chain."""
    chains = available_protein_chains(structure)
    if not chains:
        raise ValueError("The structure contains no standard amino-acid protein chains.")

    if requested_chain_id in chains:
        selected = requested_chain_id
        print(f"Using requested chain {selected} ({chains[selected]} protein residues).")
    else:
        selected = max(chains, key=chains.get)
        available = ", ".join(f"{key} ({value})" for key, value in chains.items())
        print(
            f"Requested chain {requested_chain_id!r} is unavailable. "
            f"Using largest protein chain {selected!r} ({chains[selected]} residues). "
            f"Available chains: {available}"
        )
    return selected


def extract_chain(input_pdb, output_pdb, chain_id="A"):
    """Extract a verified protein chain and return its path and selected chain ID."""
    print("\n[7] Loading protein structure...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    print("Protein structure loaded successfully.")

    selected_chain_id = select_chain(structure, chain_id)
    print(f"[8] Extracting chain {selected_chain_id}...")
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), ProteinChainSelect(selected_chain_id))

    print(f"Chain {selected_chain_id} extracted successfully.")
    print("Saved as:", output_pdb)
    return output_pdb, selected_chain_id
