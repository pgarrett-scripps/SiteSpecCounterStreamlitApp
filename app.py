from collections import Counter

import streamlit as st
import pandas as pd
from peptacular.sequence import strip_modifications, get_modifications
from peptacular.protein import find_peptide_indexes

# make wide
st.set_page_config(page_title="Site Spec Counter")

st.title('Site Spec Counter :syringe:')
st.write('This app counts the number of spectra for each modification at each site in a protein sequence')


f = st.file_uploader("Upload a file")
protein = st.text_area("Enter a protein sequence")

if f is not None and protein is not None:

    protein = protein.replace('\n', '').replace('\r', '').replace(' ', '')

    if f.name.endswith('.xlsx'):
        df = pd.read_excel(f)
    else:
        df = pd.read_csv(f)

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
