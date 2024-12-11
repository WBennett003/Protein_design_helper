from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain

class RfDiffusion_parser:
    def __init__(self, structure, target_chain=['A'], binder_chain='', new_pdb='tmp.pdb', diffusion_t=50, partial_diffusion=1,model_weights='Default'):
        self.model_weights = model_weights
        self.diffusion_t = diffusion_t
        self.partial_diffusion = partial_diffusion

        self.interchain_seperater = '/' # some uses :
        self.chain_seperater = ' '
        self.chain_terminator = ' /0'

        self.structure = structure
        self.hotspots = set()

        self.inpainting = False
        self.inpainting_configs = {
            "start_res" : [],
            "end_res" : [],
            "n_res" : []
        }

        self.target_shaved = False
        self.target_segments = {
            "start_res_idx" : [],
            "end_res_idx" :[],
            "chain" : []
        }

        self.target_chain = target_chain
        self.binder_chain = binder_chain
        self.new_pdb = new_pdb

    def add_hotspot(self, residues):
        """
        residues format will be an array of ChainIdx i.e ["A11","A14","B5"]
        """
        for res in residues:
            self.hotspots.add(res)

    def remove_hotspot(self, residues):
        """
        residues format will be an array of ChainIdx i.e ["A11","A14","B5"]
        """
        for res in residues:
            self.hotspots.remove(res)

    def reset_hotspot(self):
        self.hotspots = set()

    def add_inpainting(self, starting_residue, ending_residue, n_residues):
        self.inpainting = True
        self.inpainting_config.append({
            "start_res" : starting_residue,
            "end_res"  :ending_residue,
            "n_res" : n_residues
        })

    def reset_inpainting(self):
        self.inpainting = False
        self.inpainting_config = []

    def shave_target(self, segments):
        """
        segments of the target protein to create a binder too.
        """
        self.target_shaved = True
        
        for seg in segments:
            chain = seg[0]
            seg_idx = [int(x[1:]) for x in seg] # "A11"
            seg_idx = seg_idx.sort()
            assert seg_idx == list(range(seg_idx[0], seg_idx[-1]))
            
            self.target_segments["start_res_idx"].append(seg_idx[0])
            self.target_segments["end_res_idx"].append(seg_idx[-1])
            self.target_segments["chain"].append(chain)
            
    def reset_shave(self):
        self.target_shaved = False
        self.target_segments = []

    def create_rfdiffusion_script(self):

        if self.target_shaved:
            if self.in_painting:
                before_motif = self.inpainting_config['']
        


class sequence_designer:
    def __init__(self):
        pass

    