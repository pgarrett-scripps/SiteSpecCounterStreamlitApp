from collections import Counter

import requests
import streamlit as st
import pandas as pd
from peptacular.sequence import strip_modifications, get_modifications
from peptacular.protein import find_peptide_indexes

# make wide
st.set_page_config(page_title="Site Spec Counter")

st.title('Site Spec Counter :syringe:')
st.write('This app counts the number of spectra for each modification at each site in a protein sequence')

def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        return None
    return ''.join(response.text.split('\n')[1:])  # Remove the header line

files = st.file_uploader("Upload a file", type=["csv"], accept_multiple_files=True)

protein_id = st.text_input("Enter a protein ID (e.g. P02769)")

fetched_protein = ''
if protein_id:
    try:
        fetched_protein = fetch_sequence_from_uniprot(protein_id)
    except:
        st.error("Failed to fetch protein sequence from Uniprot")
        st.stop()

protein = st.text_area("Enter a protein sequence", value=fetched_protein)

HEADER = ['unique', 'sequence', 'spec count', 'confidence (%)', 'scan', 'charge', 'evaluation', 'fileName', 'primary score', 'DeltCN', 'M+H+(calculated)', 'M+H+(measured)', 'm/z(calculated)', 'm/z(measured)', 'ppm', 'RetTime']


if len(files) != 0 and protein is not None:

    protein = protein.replace('\n', '').replace('\r', '').replace(' ', '')

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep=',', header=None, names=HEADER)
        dfs.append(df)

    df = pd.concat(dfs)

    df['Peptide'] = [sequence[2:-2] for sequence in df['sequence']]
    df['Stripped.Peptide'] = df['Peptide'].apply(strip_modifications)
    df['Modifications'] = df['Peptide'].apply(get_modifications)
    df['Protein.Index'] = df.apply(lambda row: find_peptide_indexes(protein, row['Stripped.Peptide']), axis=1)

    # For rows which have multiple protein indexes, we need to create a new row for each protein index
    df = df.explode('Protein.Index')

    # add Protein.Index to modifications
    df['Modifications'] = df.apply(lambda row: {row['Protein.Index'] + i: mod for i, mod in row['Modifications'].items()}, axis=1)

    # make modifications a list of tuples
    df['Modifications'] = df['Modifications'].apply(lambda mods: [(site, mod) for site, mod in mods.items()])

    # explode the modifications
    df = df.explode('Modifications')

    # remove rows with no modifications
    df = df[df['Modifications'].notna()]

    # Define a Site column
    df['Site'] = df['Modifications'].apply(lambda x: x[0])

    # Define a Modification column
    df['Modification'] = df['Modifications'].apply(lambda x: x[1])

    # Define a spec count column
    st.subheader('Raw data')
    st.dataframe(df, use_container_width=True)

    # sum the spec count for each modification per site
    df = df.groupby(['fileName', 'Modification', 'Site']).agg({'spec count': 'sum'}).reset_index()

    st.subheader('Spec Counts')
    st.dataframe(df, use_container_width=True)
