#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import re
from tqdm import tqdm
import pandas as pd
from pandas import IndexSlice as IDX
from BCBio import GFF


GTF_REGEX = re.compile(r'(gene_id|transcript_id|exon_id|exon_assignment) "([^;]+)"')
COLNAMES = [
    "chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "gene_id", "transcript_id", "exon_id"
]
COLS_TO_USE = [
    "gene_id", "transcript_id", "exon_id", "start", "end", "chr", "strand"
]
COLS_TO_USE_RENAMED = [
    "gene id", "transcript id", "exon id", "exon start", "exon end", "chr", "strand"
]
COLS_TO_USE_MAPPING = dict(zip(COLS_TO_USE, COLS_TO_USE_RENAMED))


def ParseGFF(gff):
    data = []
    f = open(gff)
    for rec in GFF.parse(f):
        for gene_feature in rec.features:
            for transcript_feature in gene_feature.sub_features:
                e_number = 1
                for feature in transcript_feature.sub_features:
                    if feature.type != "exon":
                        continue
                    try:
                        data.append(
                            {
                                "gene id": gene_feature.id, 
                                "transcript id": transcript_feature.id,
                                "exon id": feature.id, 
                                "exon start": feature.location.start, 
                                "exon end": feature.location.end, 
                                "chr": rec.id, 
                                # "strand": "+" if feature.strand >= 0 else "-",
                            }
                        )
                    except:
                        continue
    f.close()

    return pd.DataFrame(data).set_index(
        ["gene id", "transcript id", "exon id"], drop=False
    )
    
def get_gene_transcript_exon_id(attrs):
    name, value = zip(*re.findall(GTF_REGEX, attrs))
    attrs_parsed = pd.Series(value, name)
    
    return attrs_parsed
    

def ParseGTF(gtf):
    data = []
    for chunk in tqdm(pd.read_csv(
        gtf, sep="\t", comment="#", header=None, chunksize=1000
    )):
        chunk.columns = COLNAMES[:-3]
        chunk = chunk.loc[chunk["feature"] == "exon"]
        if len(chunk) == 0:
            continue
        parsed_attributes = chunk["attribute"].str.findall(GTF_REGEX).apply(dict).apply(pd.Series)
        chunk = pd.DataFrame(chunk.join(parsed_attributes), columns = COLNAMES)

        # filling empty values
        chunk.loc[chunk["gene_id"].isna(), "gene_id"] = chunk.loc[chunk["gene_id"].isna(), "transcript_id"]
        chunk.loc[chunk["gene_id"].isna(), "gene_id"] = chunk.loc[chunk["gene_id"].isna(), "chr"]
        chunk.loc[chunk["transcript_id"].isna(), "transcript_id"] = chunk.loc[chunk["transcript_id"].isna(), "gene_id"]
        chunk.loc[chunk["transcript_id"].isna(), "transcript_id"] = chunk.loc[chunk["transcript_id"].isna(), "chr"]
        chunk.loc[chunk["exon_id"].isna(), "exon_id"] = chunk.loc[chunk["exon_id"].isna(), "transcript_id"]
        
        data.append(chunk.rename(columns=COLS_TO_USE_MAPPING)[COLS_TO_USE_RENAMED])
    
    data = pd.concat(data, axis=0)
    # data.sort_values(by="exon start", inplace=True)
    # data["exon number"] = data.groupby("transcript id").apply(lambda x: x.index.to_series().rank(method="first").astype("int"))
    # data.loc[data["exon id"] == data["transcript id"], "exon id"] = data["transcript id"] + "_" + data["exon number"].astype("str")

    return data.set_index(
        ["gene id", "transcript id", "exon id"], drop=False
    )
