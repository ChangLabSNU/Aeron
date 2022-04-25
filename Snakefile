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

def readfile_withending(wildcards):
	for r in READFILE_FULLNAME:
		if os.path.basename(r).split('.', maxsplit=1)[0] == wildcards.reads:
			return r
	if wildcards.reads == "refs":
		return TRANSCRIPTFILE_FULLNAME
	assert False


rule all:
	input:
		expand(OUTDIR + "/{reads}/aln.all.gam", reads=READFILE_NAME + ["refs"]),
		expand(OUTDIR + "/{reads}/aln.selected.gam", reads=READFILE_NAME + ["refs"]),
		expand(OUTDIR + "/{reads}/aln.full_length.gam", reads=READFILE_NAME + ['refs']),
		expand(OUTDIR + "/{reads}/matrix.all.txt", reads=READFILE_NAME),
		expand(OUTDIR + "/{reads}/matrixstats.txt", reads=READFILE_NAME),
		expand(OUTDIR + "/{reads}/CountMatrix.txt", reads=READFILE_NAME)

rule align_reads:
	input:
		graph = GRAPHFILE_FULLNAME,
		reads = readfile_withending
	output:
		OUTDIR + "/{reads}/aln.all.gam"
	benchmark:
		OUTDIR + "/benchmark/{reads}/aln.all.txt"
	threads: 15
	run:
		shell("mkdir -p tmp")
		shell("{ALIGNERBINPATH}/GraphAligner -g {input.graph} -f {input.reads} --try-all-seeds --seeds-mxm-length {SEEDSIZE} --seeds-mem-count {MAXSEEDHITS} --seeds-mxm-cache-prefix tmp/seedcache -a {output} -t {threads} -b {ALIGNERBANDWIDTH} {ALIGNMENTSELECTION} 1> tmp/aligner_stdout.txt 2> tmp/aligner_stderr.txt")


rule postprocess:
	input:
		all_alns = OUTDIR + "/{reads}/aln.all.gam",
		reads = readfile_withending
	output:
		selected_alns = OUTDIR + "/{reads}/aln.selected.gam",
		full_len_alns = OUTDIR + "/{reads}/aln.full_length.gam",
		summary = "tmp/{reads}/run.summary.txt",
	benchmark:
		OUTDIR + "/benchmark/{reads}.postprocess.txt"
	shell:
		"{ALIGNERBINPATH}/Postprocess {input.all_alns} {input.reads} {output.selected_alns} {output.full_len_alns} {output.summary}"

rule pick_longest:
	input:
		OUTDIR + "/{reads}/aln.selected.gam"
	output:
		OUTDIR + "/{reads}/aln.longest.gam"
	benchmark:
		OUTDIR + "/benchmark/{reads}/pick_longest.txt"
	shell:
		"{ALIGNERBINPATH}/SelectLongestAlignment {input} {output}"

rule assign_reads_to_transcripts:
	input:
		readalns = OUTDIR + "/{reads}/aln.selected.gam",
		transcripts = OUTDIR + "/refs/aln.full_length.gam",
		readfa = readfile_withending
	output:
		OUTDIR + "/{reads}/matrix.all.txt"
	benchmark:
		OUTDIR + "/benchmark/{reads}/assign_reads_to_transcript.txt"
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
		OUTDIR + "/{reads}/matrix.all.txt",
	output:
		OUTDIR + "/{reads}/matrixstats.txt"
	params:
		allmatrix = OUTDIR + "/{reads}/matrix.all.txt",
	run:
		shell("echo 'Number of reads considered:' >> {output}")
		shell("echo 'Number of transcripts considered:' >> {output}")
		shell("echo 'Number of reads which overlap a transcript:' >> {output}"),
		shell("cut -f 1 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of transcripts which overlap a read:' >> {output}"),
		shell("cut -f 2 {params.allmatrix} | sort | uniq | wc -l >> {output}"),
		shell("echo 'Number of total read-transcript overlaps:' >> {output}"),
		shell("wc -l < {params.allmatrix} >> {output}"),

rule generateCountMatrix:
	input:
		matrix = OUTDIR + "/{reads}/matrix.all.txt"
	output:
		OUTDIR + "/{reads}/CountMatrix.txt"
	benchmark:
		OUTDIR + "/benchmark/{reads}/generateCount.txt"
	shell:
		"python {SCRIPTPATH}/ThreePrime.py -m {input.matrix} >> {output}"