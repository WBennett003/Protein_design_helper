from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

# Function to filter residues based on the given list
def filter_residues(structure, residue_list):
    selected_residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Generate residue identifier string, e.g., "A13"
                res_id = f"{chain.id}{residue.id[1]}"
                if res_id in residue_list:
                    selected_residues.append(residue)
    return selected_residues

def get_protein_target(structure, chains=[], name='target'):
    new_structure = Structure.Structure(name)
    new_model = Model.Model(0) # create a new model at index 0

        # Extract the specified chain
    for model in structure:
        for chain in model:
            if chain.id in chains:
                # Detach the chain from its parent and add it to the new model
                new_chain = chain.copy()  # Make a copy to prevent modifying the original
                new_chain.detach_parent()
                new_model.add(new_chain)
                break  # Stop after finding the desired chain
        break  # Stop after processing the first model
    
    # Add the model to the new structure
    new_structure.add(new_model)
    return new_structure

def scrape_residues(pdb, target_chain=["A"], residues=[], output_file=''):
    parser = PDBParser()
    complex = parser.get_structure("Complex", pdb)

    residues_objs = filter_residues(complex, residues)

    target = get_protein_target(complex, target_chain, 'target')
    
    new_chain = Chain.Chain("B")
    for res in residues_objs:
        new_chain.add(res)

    target[0].add(new_chain)

    io = PDBIO()
    io.set_structure(target)
    io.save(output_file)


if __name__ == '__main__':
    scrape_residues('test/gC1q_l78_s920008.pdb', ['A'], ["B1","B2","B7","B35","B42","B43","B46","B49","B51","B52","B53","B54","B55","B58","B59","B61","B62","B63","B65","B66","B69","B72","B73","B76"], 
                    'test.pdb')
