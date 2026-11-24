import pandas as pd
from Bio.Seq import Seq

df = pd.read_csv("/home/wkkongaa/data/notebooks/250718_TRC_like_shRNAone_lib/trc_public.05Apr11.processed.txt", sep="\t", dtype=str)
df = df.drop(columns=["plateName", "plateBatch", "row", "col", "vectorId"]).drop_duplicates(subset=['cloneId'], keep='first').fillna('NaN').astype(str)

def concat_cloneId(df):
    # Group by 'targetSeq' and concatenate 'cloneId' values
    df_concat = df.copy()
    df_concat['cloneId'] = df_concat.groupby('targetSeq')['cloneId'].transform(lambda x: '|'.join(x))
    # Drop duplicates based on 'targetSeq', keeping the first occurrence
    return df_concat.drop_duplicates(subset=['targetSeq'], keep='first')

df_hsa = concat_cloneId(df[df["targetTaxonId"] == "9606"])
df_mmu = concat_cloneId(df[df["targetTaxonId"] == "10090"])
df_cntrl = concat_cloneId(df[df["targetTaxonId"] == "CONTROL"])
df_empty = concat_cloneId(df[df["targetTaxonId"] == "EMPTY"])

dfs = [df_hsa, df_mmu, df_cntrl]
dfs_output = []

shRNAone_5 = "aaagctcccatagtggctgtggctggat".upper() ##first nucleotide was deleted to customise the first nucleotide in the later steps
shRNAone_Loop = "cctgtttaagcttcagtggagg".upper()
shRNAone_3 = "ttctggaggcagctgatggtccgccgcc".upper() ##first nucleotide was deleted to customise the first nucleotide in the later steps


for df_type in dfs:
    df_type = df_type[df_type["targetSeq"].str.len() == 21]
    guides = [str(Seq(guide[:-4]).reverse_complement()) + "t" for guide in df_type["targetSeq"].to_list()]
    passengers = []
    shRNAone = []
    primerF = []
    primerR = []
    for guide in guides:
        
        passenger = list(str(Seq(guide).reverse_complement()))
        passenger[3] = str(Seq(passenger[3]).reverse_complement())
        passenger[7] = str(Seq(passenger[7]).reverse_complement())
        passenger[8] = str(Seq(passenger[8]).reverse_complement())
        passenger[17] = str(Seq(passenger[17]).reverse_complement())
        
        passenger = "".join(passenger)
        passengers.append(passenger)

        shRNAone.append(shRNAone_5 + guide + shRNAone_Loop + passenger + shRNAone_3)
    
    for shRNA in shRNAone:
        primerF.append(shRNA[0:71].lower())
        primerR.append(str(Seq(shRNA[-71:]).reverse_complement()).lower())
    
    df_type["shRNAone"] = shRNAone
    df_type["ForwardPrimer"] = primerF
    df_type["ReversePrimer"] = primerR
    dfs_output.append(df_type)
df_empty["shRNAone"] = ["acactcgagcac".upper()] 
df_empty["ForwardPrimer"] = ["acactcgagcac".lower()]
df_empty["ReversePrimer"] = [str(Seq("acactcgagcac").reverse_complement())]

dfs_output[0].to_csv("TRC_like_shRNAone_lib_hsa.txt", index=False, sep="\t")
dfs_output[1].to_csv("TRC_like_shRNAone_lib_mmu.txt", index=False, sep="\t")
dfs_output[2].to_csv("TRC_like_shRNAone_lib_cntrl.txt", index=False, sep="\t")
df_empty.to_csv("TRC_like_shRNAone_lib_empty.txt", index=False, sep="\t")