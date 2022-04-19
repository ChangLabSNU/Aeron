#1       havana  exon    11869   12227   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "1"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002234944"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    12613   12721   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "2"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00003582793"; exon_version "1"; tag "basic"; transcript_support_level "1";
#1       havana  exon    13221   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; transcript_id "ENST00000456328"; transcript_version "2"; exon_number "3"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-202"; transcript_source "havana"; transcript_biotype "processed_transcript"; exon_id "ENSE00002312635"; exon_version "1"; tag "basic"; transcript_support_level "1";

import os
import re
from tqdm import tqdm
from optparse import OptionParser
import numpy as np
import pandas as pd
from Bio import SeqIO
from ParseGTF import ParseGTF, ParseGFF
        
class ParseOptions():
	def getoptions(self):
		parser = OptionParser()
		wd = os.getcwd()
		parser.add_option('-e', '--genome',dest='ef',help='File containing genome sequences divided into chromosomes', action="store")
		parser.add_option('-g', '--gtf',dest='gf',help='GTF file for the genome', action="store")
		parser.add_option('-o', '--out',dest='out',help='Output Folder', action="store", default=wd)
		(options, args) = parser.parse_args()
		return options

class SequenceAnalyser():
	def FParser(self, seq):
		Sequences = {}
		seqparse = SeqIO.parse(seq, 'fasta')
		
		for seq in seqparse:
			seq_id, sequence = seq.id, seq.seq
			
			Sequences[seq_id] = sequence
		return Sequences


def GraphBuild(exons: pd.DataFrame, Sequences):
	chr = exons["chr"][0]
	splice_sites = exons[["transcript id", "exon start", "exon end"]].melt("transcript id").sort_values("value")
	nodes = getNodePositions(splice_sites)
	# nodes = getNodeIDs(nodes, exons)
	nodes_out, connections_out = getNodeConnections(nodes, Sequences[chr])
	
	return nodes_out, connections_out


# def getNodeStartEnd(two_splice_sites: pd.DataFrame):
# 	node = two_splice_sites["value"]
# 	node.index = ["start", "end"]
# 	node.rename(two_splice_sites["transcript id"].iloc[1], inplace=True)
# 	if node["start"] == node["end"]:
# 		node.loc[["start", "end"]] = np.nan
# 	elif two_splice_sites["variable"].iloc[0] and two_splice_sites["variable"].iloc[1]: # start start
# 		node["end"] -= 1
# 	elif not (two_splice_sites["variable"].iloc[0] or two_splice_sites["variable"].iloc[1]): # end end
# 		node["start"] += 1
# 	elif two_splice_sites["variable"].iloc[0] and not two_splice_sites["variable"].iloc[1]: # start end
# 		pass
# 	else:
# 		node.loc[["start", "end"]] = np.nan

# 	return node


def getNodePositions(splice_sites: pd.DataFrame):
	nodes = []
	splice_sites["variable"] = splice_sites["variable"] == "exon start"
	for two_splice_sites in splice_sites.rolling(2):
		if len(two_splice_sites) < 2:
			continue
		node = two_splice_sites["value"].copy()
		node.index = ["start", "end"]
		node.rename(
			f'{two_splice_sites["transcript id"].iloc[1]}-{node["start"]}', inplace=True
		)
		if node["start"] == node["end"]:
			continue
		if two_splice_sites["variable"].iloc[0] and two_splice_sites["variable"].iloc[1]: # start start
			node["end"] -= 1
		elif not (two_splice_sites["variable"].iloc[0] or two_splice_sites["variable"].iloc[1]): # end end
			node["start"] += 1
		elif two_splice_sites["variable"].iloc[0] and not two_splice_sites["variable"].iloc[1]: # start end
			pass
		else:
			continue
		nodes.append(node)

	return pd.DataFrame(nodes).dropna().astype("int")


def getNodeConnections(nodes: pd.DataFrame, sequence: str):
	nodes_out = [
		f'S	{idx}	{sequence[node["start"] - 1 : node["end"]]}' 
		for idx, node in nodes.iterrows()
	]
	connections_out = []
	nids = nodes.index.to_list()
	for idx, nid_src in enumerate(nids):
		for nid_dest in nids[idx + 1:]:
			connections_out.append(f"L	{nid_src}	+	{nid_dest}	+	0M")

	return nodes_out, connections_out
	

def sanityCheck(nodes, connections):
	if len(nodes)==0:
		print("Warning: No sequence information added")
	if len(connections)==0:
		print("Warning: No connection information added")
		

if __name__ == "__main__":
	po  = ParseOptions().getoptions()
	sq  = SequenceAnalyser() 
	fn  = po.gf
	seq = po.ef
	out = po.out

	f = open(out, "w")
	nodes = []
	connections = []	
			
	print("Reading Sequences")
	fasta = sq.FParser(seq)
	print("Done reading sequences")

	print("Reading and processing gtf")
	pg = ParseGTF(fn) if fn.endswith(".gtf") else ParseGFF(fn)
	print("Done reading gtf")

	print("Building graph")
	genes = pg.index.unique("gene id")
	for gene in tqdm(genes):
		nodes_of_gene, connections_of_gene = GraphBuild(pg.loc[gene], fasta)
		nodes += nodes_of_gene
		connections += connections_of_gene

	f.write("\n".join(nodes))
	f.write("\n")
	f.write("\n".join(connections))

	sanityCheck(nodes, connections)

	f.close()
