import urllib.request
import urllib.parse
import sys
import subprocess
import os
import argparse

try:
    import yaml
except ImportError:
    print("Installing PyYAML for config loading...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml"])
    import yaml

try:
    import pandas as pd
except ImportError:
    print("Installing pandas and openpyxl...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "openpyxl"])
    import pandas as pd

def get_biomart_rRNA(dataset):
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "gene_biotype" />
    </Dataset>
</Query>"""
    url = "http://www.ensembl.org/biomart/martservice"
    data = urllib.parse.urlencode({'query': xml}).encode('utf-8')
    req = urllib.request.Request(url, data)
    try:
        response = urllib.request.urlopen(req)
        result = response.read().decode('utf-8')
        ribo_genes = set()
        for line in result.strip().split('\n'):
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    gene_id = parts[0]
                    biotype = parts[1]
                    if 'rRNA' in biotype or biotype == 'rRNA' or biotype == 'Mt_rRNA':
                        ribo_genes.add(gene_id)
        return ribo_genes
    except Exception as e:
        print(f"Failed to query Biomart for {dataset}: {e}")
        return set()

def calc_prop(file_path, ribo_set):
    with open(file_path, 'r') as f:
        total_counts = None
        ribo_counts = None
        samples = []
        for line in f:
            if line.startswith('#'):
                continue
            if line.startswith('Geneid'):
                samples = line.strip().split('\t')[6:]
                total_counts = [0] * len(samples)
                ribo_counts = [0] * len(samples)
                continue
            parts = line.strip().split('\t')
            gene_id = parts[0]
            counts = [int(x) for x in parts[6:]]
            for i, c in enumerate(counts):
                total_counts[i] += c
                if gene_id in ribo_set:
                    ribo_counts[i] += c
    
    data = []
    for s, t, r in zip(samples, total_counts, ribo_counts):
        s_name = os.path.basename(s.strip())
        prop = (r / t * 100) if t > 0 else 0
        data.append({
            "Sample": s_name,
            "Total_Mapped_Reads": t,
            "Ribosomal_Reads": r,
            "Percentage_Ribosomal": round(prop, 2)
        })
    return pd.DataFrame(data)

def main():
    parser = argparse.ArgumentParser(description="Calculate rRNA proportions from featureCounts matrices using a YAML config.")
    parser.add_argument('config', help='Path to the YAML configuration file')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    out_file = config.get('output_file', 'rRNA_proportions_output.xlsx')
    runs = config.get('runs', [])

    if not runs:
        print("No runs defined in the config. Exiting.")
        sys.exit(1)

    results = {}
    for run in runs:
        species_name = run.get('species', 'unknown')
        dataset = run.get('biomart_dataset')
        matrix_file = run.get('matrix_file')

        if not dataset or not matrix_file:
            print(f"Skipping run for {species_name}: missing 'biomart_dataset' or 'matrix_file'")
            continue
        
        if not os.path.exists(matrix_file):
            print(f"Skipping run for {species_name}: file {matrix_file} does not exist.")
            continue

        print(f"Querying Ensembl Biomart for {species_name} rRNA features ({dataset})...")
        rRNA_genes = get_biomart_rRNA(dataset)
        print(f"Got {len(rRNA_genes)} {species_name} rRNA genes.")

        print(f"Processing count matrix for {species_name}: {matrix_file}")
        df = calc_prop(matrix_file, rRNA_genes)
        results[species_name] = df

    if not results:
        print("No valid results to save.")
        sys.exit(1)

    print(f"Saving to {out_file}...")
    with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
        for species_name, df in results.items():
            # Excel sheet names have a max length of 31 characters
            sheet_name = f"{species_name}_rRNA"[:31]
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print("Done!")

if __name__ == "__main__":
    main()
