import urllib.request
import urllib.parse
import sys
import os
import argparse

try:
    import yaml
except ImportError:
    import subprocess
    print("Installing PyYAML for config loading...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyyaml"])
    import yaml

try:
    import pandas as pd
except ImportError:
    import subprocess
    print("Installing pandas and openpyxl...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "openpyxl"])
    import pandas as pd

def get_biomart_annotations(dataset):
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "external_gene_name" />
        <Attribute name = "description" />
        <Attribute name = "gene_biotype" />
    </Dataset>
</Query>"""
    url = "http://www.ensembl.org/biomart/martservice"
    data = urllib.parse.urlencode({'query': xml}).encode('utf-8')
    req = urllib.request.Request(url, data)
    
    print(f"Querying Biomart for {dataset} annotations...")
    try:
        response = urllib.request.urlopen(req)
        result = response.read().decode('utf-8')
        
        from io import StringIO
        df_annot = pd.read_csv(StringIO(result), sep='\t')
        
        df_annot.columns = ['Geneid', 'geneSymbol', 'Description', 'biotype']
        df_annot = df_annot.drop_duplicates(subset=['Geneid'])
        
        return df_annot
    except Exception as e:
        print(f"Failed to query Biomart for {dataset}: {e}")
        return pd.DataFrame(columns=['Geneid', 'geneSymbol', 'Description', 'biotype'])

def annotate_matrix(matrix_path, annot_df):
    print(f"Loading count matrix: {matrix_path}")
    
    df_matrix = pd.read_csv(matrix_path, sep='\t', comment='#')
    
    print("Merging annotations...")
    df_merged = pd.merge(df_matrix, annot_df, on='Geneid', how='left')
    
    cols = df_merged.columns.tolist()
    
    cols.remove('geneSymbol')
    cols.remove('Description')
    cols.remove('biotype')
    
    cols.insert(1, 'geneSymbol')
    cols.insert(2, 'Description')
    cols.insert(3, 'biotype')
    
    df_final = df_merged[cols]
    
    dir_name = os.path.dirname(matrix_path)
    base_name = os.path.basename(matrix_path)
    name_no_ext = os.path.splitext(base_name)[0]
    out_file = os.path.join(dir_name, f"annotated_{name_no_ext}.xlsx")
    
    print(f"Saving to {out_file}...")
    df_final.to_excel(out_file, index=False, engine='openpyxl')
    print(f"Saved successfully!")

def main():
    parser = argparse.ArgumentParser(description="Annotate featureCounts matrices with BioMart Gene Symbols, Descriptions, and Biotypes using a YAML config.")
    parser.add_argument('config', help='Path to the YAML configuration file')
    args = parser.parse_args()

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    runs = config.get('runs', [])

    if not runs:
        print("No runs defined in the config. Exiting.")
        sys.exit(1)

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
            
        annot_df = get_biomart_annotations(dataset)
        
        if not annot_df.empty:
            annotate_matrix(matrix_file, annot_df)
        else:
            print(f"Warning: Annotation dataframe was empty for {dataset}.")

if __name__ == "__main__":
    main()
