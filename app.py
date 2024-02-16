from collections import Counter

import streamlit as st
import pandas as pd
from peptacular.sequence import strip_modifications, get_modifications
from peptacular.protein import find_peptide_indexes

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

    # now for each file, I want to accumulate all modification
    mod_per_file = {}
    for file in df['fileName']:
        mod_counter = Counter()
        file_df = df[df['fileName'] == file]
        spec_count = file_df.shape[0]
        for mods, spec_count in file_df[['Modifications', 'spec count']].values:
            for i, mod in mods.items():
                mod_counter[mod] += spec_count

        mod_per_file[file] = mod_counter

    mod_data = []
    for file, mod_counter in mod_per_file.items():
        for mod, count in mod_counter.items():
            mod_data.append({'file': file, 'mod': mod, 'count': count})

    mod_df = pd.DataFrame(mod_data)
    st.dataframe(mod_df)



