import argparse
import pathlib 

def merge_binder_target(fasta_seq, target_seq, binder_chain):
    chains = ["A", "B", "C"]
    binder_idx = chains.index(binder_chain)
    binder_seq = fasta_seq.split(':')[binder_idx]

    target_seq = target_seq.split(':')
    fasta_line = target_seq[:binder_idx] + ':' + binder_seq + ":" + target_seq[binder_seq:]
    return fasta_line

def merge_batch(fasta_files, target_seq, binder_chain, output_file):
    for i,fasta_file in enumerate(fasta_files):
        with open(fasta_file, 'r') as f:
            fasta = f.read().split('\n')
        
        new_txt = fasta[0] + '\n' + fasta[1] + '\n'
        for i,line in enumerate(fasta[2:]):
            if line[0] != '>':
                fasta += merge_binder_target(line, target_seq, binder_chain) + '\n'
            else:
                line[0] == '>':
                fasta += f'{line}_design_{i}\n'

    with open(output_file, 'w+') as f:
        f.write(fasta)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_path')
    parser.add_argument('--target_seq')
    parser.add_argument('--output_file')
    parser.add_argument('--chain', default='B')
    args = parser.parse_args()
    
    path = pathlib.Path(args.fasta_path)
    fasta_files = list(path.rglob("*.pdb"))

    merge_batch(fasta_files, args.target_seq, args.chain, args.output_file)