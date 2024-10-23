import pandas as pd
from Bio import SeqIO
import sys
import os

script_directory = os.path.dirname(os.path.abspath(__file__))


mirs_fasta = sys.argv[1]  # First argument
line_name = sys.argv[2]  # Second argumen

def parse_fasta_to_dataframe(fasta_file):
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        records.append((seq_id, sequence))
    
    df = pd.DataFrame(records, columns=["Sequence_ID", "Sequence"])
    return df


fasta_file = mirs_fasta
df = parse_fasta_to_dataframe(fasta_file)

df_precursor = df[~df['Sequence_ID'].str.contains('mature', case=False, na=False)]
df_precursor = df_precursor[~df_precursor['Sequence_ID'].str.contains('star', case=False, na=False)]
df_precursor.rename(columns={'Sequence_ID': 'Cluster', 'Sequence': 'Precursor_sequence'}, inplace=True)


df_star = df[df['Sequence_ID'].str.contains('star', case=False, na=False)]
df_star.rename(columns={'Sequence_ID': 'Cluster', 'Sequence': 'Star_sequence'}, inplace=True)

df_mature = df[df['Sequence_ID'].str.contains('mature', case=False, na=False)]
df_mature.rename(columns={'Sequence_ID': 'Cluster', 'Sequence': 'Mature_sequence'}, inplace=True)


f"{line_name}_results.txt"

df_precursor.to_csv(os.path.join(script_directory, f"{line_name}_precursors.csv"), index=False)
df_star.to_csv(os.path.join(script_directory, f"{line_name}_2gg_star.csv"), index=False)
df_mature.to_csv(os.path.join(script_directory, f"{line_name}_2gg_mature.csv"), index=False)

df_full = pd.read_csv(os.path.join(script_directory,'raw_counts_named_sequences.csv'))
df_prec = pd.read_csv(os.path.join(script_directory,f"{line_name}_precursors.csv"))
df_star2 = pd.read_csv(os.path.join(script_directory,f"{line_name}_2gg_star.csv"))
df_mature2 = pd.read_csv(os.path.join(script_directory,f"{line_name}_2gg_mature.csv"))

df_combined = pd.concat([df_full, df_mature2], axis=1)
df_combined = pd.concat([df_combined, df_star2], axis=1)
df_combined = pd.concat([df_combined, df_prec], axis=1)

df_combined.to_csv(f"{line_name}_combined.csv", index=False)
df_combined_short = df_combined[['known_miRNAs', 'MajorRNA.x', 'Locus', "Strand", 'Mature_sequence', 'Star_sequence', 'Precursor_sequence']]
df_combined_short.to_csv(os.path.join(script_directory,f"{line_name}_combined_short.csv"), index=False)
df_combined_short['Combined'] = df_combined_short.apply(lambda x: f"{x['MajorRNA.x']}_{x['Locus']}_{x['Strand']}", axis=1)
df_combined_short.to_csv(os.path.join(script_directory,f"{line_name}_combined_short_combined.csv"), index=False)
