import os

configfile: "config.yaml"
VGPATH = config["vgpath"]
SCRIPTPATH=config["scripts"]
ALIGNERBINPATH = config["binaries"]
READFILE_FULLNAME = config["reads"]
TRANSCRIPTFILE_FULLNAME = config["transcripts"]
GRAPHFILE_FULLNAME = config["graph"]
SEEDSIZE = config["seedsize"]
MAXSEEDHITS = config["maxseeds"]
ALIGNERBANDWIDTH = config["aligner_bandwidth"]
GTFPATH = config["gtffile"]
ALIGNMENTSELECTION = config["alignment_selection"]
ECUTOFF = config["alignment_E_cutoff"]
OUTDIR = config["outdir"]

if isinstance(READFILE_FULLNAME, str): READFILE_FULLNAME = [READFILE_FULLNAME]

READFILE_NAME = [os.path.basename(path).split('.', maxsplit=1)[0] for path in READFILE_FULLNAME]
TRANSCRIPTFILE_NAME = os.path.basename(TRANSCRIPTFILE_FULLNAME).split('.', maxsplit=1)[0]
GRAPHFILE_NAME = os.path.basename(GRAPHFILE_FULLNAME).split('.', maxsplit=1)[0]

def readfile_withending(wildcards):
	for r in READFILE_FULLNAME:
		if os.path.basename(r).split('.', maxsplit=1)[0] == wildcards.reads:
			return r
	if wildcards.reads == os.path.basename(TRANSCRIPTFILE_NAME):
		return TRANSCRIPTFILE_FULLNAME
	assert False

def graph_path(wildcards):
	if wildcards.graph == os.path.basename(GRAPHFILE_FULLNAME).split(".", maxsplit=1)[0]:
		return GRAPHFILE_FULLNAME
	assert False

rule all:
	input:
		expand(OUTDIR + "/aln.{reads}.{graph}.all.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/aln.{reads}.{graph}.selected.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/aln.{reads}.{graph}.full_length.gam", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		#expand(OUTDIR + "/alignmentstats.{reads}.{graph}.txt", reads=READFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/aln.{transcript}.{graph}.full_length.gam", transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		#expand(OUTDIR + "/alignmentstats_{transcript}_{graph}.txt", transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/matrix.{reads}.{transcript}.{graph}.all.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		#expand(OUTDIR + "/matrix.{reads}.{transcript}.{graph}.bestmatch.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/matrixstats.{reads}.{transcript}.{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME),
		expand(OUTDIR + "/CountMatrix.{reads}.{transcript}.{graph}.txt", reads=READFILE_NAME, transcript=TRANSCRIPTFILE_NAME, graph=GRAPHFILE_NAME)

rule align:
	input:
		graph = graph_path,
		reads = readfile_withending
	output:
		OUTDIR + "/aln.{reads}.{graph}.all.gam"
	benchmark:
		OUTDIR + "/benchmark/aln.{reads}.{graph}.all.txt"
	threads: 15
	run:
		shell("mkdir -p tmp")
		shell("{ALIGNERBINPATH}/GraphAligner -g {input.graph} -f {input.reads} --try-all-seeds --seeds-mxm-length {SEEDSIZE} --seeds-mem-count {MAXSEEDHITS} --seeds-mxm-cache-prefix tmp/seedcache -a {output} -t {threads} -b {ALIGNERBANDWIDTH} {ALIGNMENTSELECTION} 1> tmp/aligner_stdout.txt 2> tmp/aligner_stderr.txt")

rule postprocess:
	input:
		all_alns = OUTDIR + "/aln.{reads}.{graph}.all.gam",
		reads = readfile_withending
	output:
		selected_alns = OUTDIR + "/aln.{reads}.{graph}.selected.gam",
		full_len_alns = OUTDIR + "/aln.{reads}.{graph}.full_length.gam",
		summary = "tmp/run.{reads}.{graph}.summary.txt",
	benchmark:
		OUTDIR + "/benchmark/postprocess.{reads}.{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/Postprocess {input.all_alns} {input.reads} {output.selected_alns} {output.full_len_alns} {output.summary}"

rule pick_longest:
	input:
		OUTDIR + "/aln.{reads}.{graph}.selected.gam"
	output:
		OUTDIR + "/aln.{reads}.{graph}.longest.gam"
	benchmark:
		OUTDIR + "/benchmark/pick_longest.{reads}.{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/SelectLongestAlignment {input} {output}"

rule assign_reads_to_transcripts:
	input:
		readalns = OUTDIR + "/aln.{reads}.{graph}.selected.gam",
		transcripts = OUTDIR + "/aln.{transcripts}.{graph}.full_length.gam",
		readfa = readfile_withending
	output:
		OUTDIR + "/matrix.{reads}.{transcripts}.{graph}.all.txt"
	benchmark:
		OUTDIR + "/benchmark/assign_reads_to_transcript.{reads}.{transcripts}.{graph}.txt"
	shell:
		"{ALIGNERBINPATH}/AlignmentSubsequenceIdentity {input.transcripts} {input.readalns} {input.readfa} > {output}"

#rule find_best_assignments:
#	input:
#		OUTDIR + "/matrix_{runid}_all.txt"
#	output:
#		OUTDIR + "/matrix_{runid}_bestmatch.txt"
#	benchmark:
#		OUTDIR + "/benchmark/bestmatch_{runid}.txt"
#	shell:
#		"{SCRIPTPATH}/find_matrix_bestmatch.py {input} 0.4 0.95 > {output}"

rule output_assignment_statistics:
	input:
		OUTDIR + "/matrix.{reads}.{transcripts}.{graph}.all.txt",
		#OUTDIR + "/matrix.{reads}.{transcripts}.{graph}_bestmatch.txt",
		#OUTDIR + "/matrix.{reads}.{transcripts}.{graph}_unambiguous.txt",
		#read_stdout = "tmp/aligner_stdout.{reads}.{graph}.txt",
		#transcript_stdout = "tmp/aligner_stdout.txt"
	output:
		OUTDIR + "/matrixstats.{reads}.{transcripts}.{graph}.txt"
	params:
		allmatrix = OUTDIR + "/matrix.{reads}.{transcripts}.{graph}.all.txt",
		#bestmatchmatrix = OUTDIR + "/matrix.{reads}.{transcripts}.{graph}_bestmatch.txt",
		#unambiguousmatrix = OUTDIR + "/matrix.{reads}.{transcripts}.{graph}_unambiguous.txt"
	run:
		shell("echo 'Number of reads considered:' >> {output}")
		#shell("grep 'Reads with an alignment:' < {input.read_stdout} | cut -d ' ' -f 5 >> {output}")
		shell("echo 'Number of transcripts considered:' >> {output}")
		#shell("grep 'Output end-to-end alignments:' < {input.transcript_stdout} | cut -d ' ' -f 3 >> {output}")
		shell("echo 'Number of reads which overlap a transcript:' >> {output}"),
		shell("cut -f 1 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of transcripts which overlap a read:' >> {output}"),
		shell("cut -f 2 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of total read-transcript overlaps:' >> {output}"),
		shell("wc -l < {params.allmatrix} >> {output}"),
		#shell("echo 'Number of bestmatch read-transcript overlaps:' >> {output}"),
		#shell("wc -l < {params.bestmatchmatrix} >> {output}"),
		#shell("echo 'Number of unambiguous overlaps:' >> {output}")
		#shell("wc -l < {params.unambiguousmatrix} >> {output}")

rule generateCountMatrix:
	input:
		matrix = OUTDIR + "/matrix.{reads}.{transcripts}.{graph}.all.txt"
	output:
		OUTDIR + "/CountMatrix.{reads}.{transcripts}.{graph}.txt"
	benchmark:
		OUTDIR + "/benchmark/generateCount.{reads}.{transcripts}.{graph}.txt"
	shell:
		"python {SCRIPTPATH}/ThreePrime.py -m {input.matrix} >> {output}"