import pandas as pd
import json
from pathlib import Path
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_gtf',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/WBcel235/Caenorhabditis_elegans.WBcel235.112.gtf')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/reference_genomes/WBcel235')
    
    args = parser.parse_args()
    
    df = pd.read_csv(Path(args.input_gtf),
                     sep='\t',
                     comment='#',
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    genes = df[df['feature']=='gene']
    protein_coding_genes: list[str] = []
    for x in genes['attribute']:
        if 'protein_coding' in x:
            protein_coding_genes.append(x.split(';')[0].replace('gene_id "', '').replace('"', '').strip())

    with open(Path(args.output_folder) / 'protein_coding_genes.json', 'w') as file:
        json.dump(sorted(protein_coding_genes), file)


