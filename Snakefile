import os
import glob

configfile: "config.yaml"

DATA=config['data'] if 'data' in config else "data"
PIPELINE=os.path.dirname(workflow.snakefile)

def get_config(name, default=None):
    return config[name] if name in config else default

def flatten(alist):
    if alist == []:
        return []
    elif type(alist) is not list:
        return [alist]
    else:
        return flatten(alist[0]) + flatten(alist[1:])
    
def get_input_files(name):
    if name not in config:
        return []
    l = config[name]
    if not isinstance(l, list):
        l = [l]
    l = [os.path.expanduser(p) for p in l]
    l = [glob.glob(p) for p in l]
    l = flatten(l)
    if l == []:
        raise FileNotFoundError('input files for \''+name+'\' not found.')
    return l

# ------------------------------------------------------------------------
# Run Unicycler
# ------------------------------------------------------------------------

rule run_unicycler:
    input:
        short_r1=get_input_files('trimmed_R1_fq'),
        short_r2=get_input_files('trimmed_R2_fq'),
        long_reads=get_input_files('filtered_long_fq')
    output: DATA+"/unicycler/assembly.fasta"
    threads: 9999
    conda: "envs/unicycler.yaml"
    shell:
        """
        unicycler -t {threads} \
          -1 {input.short_r1} \
          -2 {input.short_r2} \
          -l {input.long_reads} \
          -o $(dirname {output})
        """

# ------------------------------------------------------------------------
# Compare input genome and Unicycler results with `dnadiff`
# ------------------------------------------------------------------------

rule run_dnadiff:
    input:
        raw=get_input_files('assembly_fa'),
        unic=DATA+"/unicycler/assembly.fasta"
    output: DATA+"/dnadiff/out.report"
    conda: "envs/mummer4.yaml"
    shell:
        """
        dir=$(dirname {output})
        cp {input.raw} $dir/raw.fasta
        cp {input.unic} $dir/unicycler.fasta
        cd $dir
        dnadiff raw.fasta unicycler.fasta
        """

# ------------------------------------------------------------------------
# "Normalize" the genome.
# ------------------------------------------------------------------------

rule normalize_genome:
    input: get_input_files('assembly_fa')
    output: DATA+"/normalized/normalized.fasta"
    params:
        strain=get_config('strain'),
        version=get_config('version',''),
        linear=get_config('linear_contigs','')
    conda: "envs/normalize.yaml"
    shell:
        """
        dir=$(dirname {output})

        cat {input} \
            | {PIPELINE}/scripts/dephix \
                  > $dir/unnormalized.fasta

        {PIPELINE}/scripts/normalize-assembly \
            -d $dir/tmp \
            -f {PIPELINE}/inputs/starts.faa \
            -l "{params.linear}" \
            $dir/unnormalized.fasta {params.strain}{params.version}_ \
            > {output}
        """

# ------------------------------------------------------------------------
# Run PGAP
# ------------------------------------------------------------------------

if get_config('pgap_dir') != None:

    rule run_pgap:
        input:
            genome=DATA+"/normalized/normalized.fasta"
        output:
            faa=DATA+"/pgap/annot.faa",
            fna=DATA+"/pgap/annot.fna",
            gbk=DATA+"/pgap/annot.gbk",
            gff=DATA+"/pgap/annot.gff",
        params:
            pgap_dir=get_input_files('pgap_dir'),
            strain=get_config('strain'),
            version=get_config('version',''),
            genus=get_config('genus','FIXME'),
            species=get_config('species','FIXME'),
            pgap_args=' '.join(get_config('pgap_args',[]))
        threads: 9999
        shell:
            """
            dir=$(dirname {output.gbk})
            rm -rf $dir/tmp
            {params.pgap_dir}/pgap.py \
                --genome {input.genome} \
                --organism "{params.genus} {params.species}" \
                --output $dir/tmp \
                --taxcheck --report-usage-false --quiet --no-internet \
                --docker apptainer --cpus {threads} --no-self-update
            mv $dir/tmp/* $dir/.
            rmdir $dir/tmp
            """

# ------------------------------------------------------------------------
# Run Prokka
# ------------------------------------------------------------------------

rule run_prokka:
    input:
        genome=DATA+"/normalized/normalized.fasta"
    output:
        faa=DATA+"/prokka/output.faa",
        fna=DATA+"/prokka/output.fna",
        gbk=DATA+"/prokka/output.gbk",
        gff=DATA+"/prokka/output.gff",
    params:
        strain=get_config('strain'),
        version=get_config('version',''),
        gram=get_config('gram'),
        genus=get_config('genus','FIXME'),
        species=get_config('species','FIXME'),
    threads: 9999
    conda: "envs/prokka.yaml"
    shell:
        """
        dir=$(dirname {output.gbk})
        rm -rf $dir/tmp
        prokka --cpus {threads} --quiet \
                 --outdir $dir/tmp \
                 --prefix output \
                 --genus {params.genus} \
                 --species {params.species} \
                 --strain {params.strain} \
                 --locustag {params.strain}{params.version}_prokka \
                 --rfam --addgenes \
                 {input.genome}
        mv $dir/tmp/* $dir/.
        rmdir $dir/tmp
        """

# ------------------------------------------------------------------------
# Run Bakta
# ------------------------------------------------------------------------

if get_config('bakta_db') != None:
    rule run_bakta:
        input:
            genome=DATA+"/normalized/normalized.fasta"
        output:
            faa=DATA+"/bakta/output.faa",
            fna=DATA+"/bakta/output.fna",
            gbk=DATA+"/bakta/output.gbff",
            gff=DATA+"/bakta/output.gff3",
        params:
            db=get_config('bakta_db'),
            strain=get_config('strain'),
            version=get_config('version',''),
            gram=get_config('gram'),
            genus=get_config('genus','FIXME'),
            species=get_config('species','FIXME'),
        threads: 9999
        conda: "envs/bakta.yaml"
        shell:
            """
            dir=$(dirname {output.gbk})
            bakta --db {params.db} \
                  --prefix output \
                  --output $dir --force \
                  --genus {params.genus} \
                  --species {params.species} \
                  --strain {params.strain} \
            	  --complete \
            	  --threads {threads} \
	 	  --keep-contig-headers \
	 	  --compliant \
                  {input.genome}
            """
    
# ------------------------------------------------------------------------
# The final genome
# ------------------------------------------------------------------------

if get_config('pgap_dir') != None:
    rule make_final_using_pgap:
        input:
            faa=DATA+"/pgap/annot.faa",
            fna=DATA+"/pgap/annot.fna",
            gbk=DATA+"/pgap/annot.gbk",
            gff=DATA+"/pgap/annot.gff",
        output:
            faa=DATA+"/final.faa",
            fna=DATA+"/final.fna",
            gbk=DATA+"/final.gbk",
            gff=DATA+"/final.gff",
        shell:
            """
            cp {input.faa} {output.faa}
            cp {input.fna} {output.fna}
            cp {input.gbk} {output.gbk}
            cp {input.gff} {output.gff}
            """
elif get_config('bakta_db') != None:
    rule make_final_using_bakta:
        input:
            faa=DATA+"/bakta/output.faa",
            fna=DATA+"/bakta/output.fna",
            gbk=DATA+"/bakta/output.gbff",
            gff=DATA+"/bakta/output.gff3",
        output:
            faa=DATA+"/final.faa",
            fna=DATA+"/final.fna",
            gbk=DATA+"/final.gbk",
            gff=DATA+"/final.gff",
        shell:
            """
            cp {input.faa} {output.faa}
            cp {input.fna} {output.fna}
            cp {input.gbk} {output.gbk}
            cp {input.gff} {output.gff}
            """
else:
    rule make_final_using_prokka:
        input:
            faa=DATA+"/prokka/output.faa",
            fna=DATA+"/prokka/output.fna",
            gbk=DATA+"/prokka/output.gbk",
            gff=DATA+"/prokka/output.gff",
        output:
            faa=DATA+"/final.faa",
            fna=DATA+"/final.fna",
            gbk=DATA+"/final.gbk",
            gff=DATA+"/final.gff",
        shell:
            """
            cp {input.faa} {output.faa}
            cp {input.fna} {output.fna}
            cp {input.gbk} {output.gbk}
            cp {input.gff} {output.gff}
            """

# ------------------------------------------------------------------------
# Run BUSCO
# ------------------------------------------------------------------------

rule run_busco:
    input: DATA+"/final.faa"
    output: DATA+"/busco/done.txt"
    params:
        lineage_arg = ('--lineage_dataset '+config['busco_lineage']) if 'busco_lineage' in config else '--auto-lineage-prok'
    threads: 9999
    conda: "envs/busco.yaml"
    shell:
        """
        dir=$(dirname {output})
        busco \
            -q \
            -i {input} \
            -o $dir/output \
            -m proteins \
            {params.lineage_arg} \
            -c {threads} \
            --download_path $dir/downloads
        touch {output}
        """

# ------------------------------------------------------------------------
# Generate BUSCO summary
# ------------------------------------------------------------------------

rule generate_busco_summary:
    input: DATA+"/busco/done.txt"
    output: DATA+"/busco/report.txt"
    params:
        strain=get_config('strain'),
    shell:
        """
        dir=$(dirname {output})
        (
            cd $dir

            echo -e "Name\tdb\tC\tS\tD\tF\tM\tn"
            echo -n {params.strain}
            egrep '^'$'\t''C:' /dev/null output/short_summary.specific.*.txt \
                | sed -r \
                      -e 's/^output//' \
                      -e 's|/short_summary.specific.(.+)_odb10.output.txt:|\t\\1|' \
                      -e 's/[ \t]+$//' \
                      -e 's/C:([0-9.%]+)\\[S:([0-9.%]+),D:([0-9.%]+)\\],F:([0-9.%]+),M:([0-9.%]+),n:([0-9.%]+)/\\1\t\\2\t\\3\t\\4\t\\5\t\\6/'
        ) > {output}
        """

# ------------------------------------------------------------------------
# Compute stats
# ------------------------------------------------------------------------

rule make_stats:
    input:
        raw_long=get_input_files('raw_long_fq'),
        filtered_long=get_input_files('filtered_long_fq'),
        raw_r1=get_input_files('raw_R1_fq') if 'raw_R1_fq' in config else [],
        raw_r2=get_input_files('raw_R2_fq') if 'raw_R2_fq' in config else [],
        trimmed_r1=get_input_files('trimmed_R1_fq') if 'trimmed_R1_fq' in config else [],
        trimmed_r2=get_input_files('trimmed_R2_fq') if 'trimmed_R2_fq' in config else [],
        final_fna=DATA+"/final.fna",
        final_gff=DATA+"/final.gff",
    output: DATA+"/stats.tsv"
    params:
        args="-q -s",
        strain=get_config('strain'),
        version=get_config('version',''),
    threads: 9999
    conda: "envs/make_stats.yaml"
    shell:
        """
        {PIPELINE}/scripts/compute-assembly-stats \
            -t {threads} \
            {params.args} -S {params.strain}{params.version} \
            {input.raw_long} {input.filtered_long} \
            {input.raw_r1} {input.raw_r2} \
            {input.trimmed_r1} {input.trimmed_r2} \
            {input.final_fna} \
            {input.final_gff} \
            | tee {output}
        """

# ------------------------------------------------------------------------
# Generate the summary
# ------------------------------------------------------------------------

rule make_summary:
    input:
        dnadiff_report=(DATA+"/dnadiff/out.report" if 'skip_unicycler' not in config and 'trimmed_R1_fq' in config else []),
        stats_tsv=DATA+"/stats.tsv",
        busco_txt=DATA+"/busco/report.txt",
        _pgap_gbk=(DATA+"/pgap/annot.gbk" if 'pgap_dir' in config else []),
        _prokka_gbk=DATA+"/prokka/output.gbk",
        _bakta_gbk=(DATA+"/bakta/output.gbff" if 'bakta_db' in config else []),
        _final_gbk=DATA+"/final.gbk",
    output: DATA+"/summary-post-assembly.log"
    conda: "envs/datamash.yaml"
    shell:
        """
        (
            if [ "{input.dnadiff_report}" ] ; then
                echo 
                echo === autocycler vs. unicycler ===
                head -n13 {input.dnadiff_report}  | tail -n+4
            fi
            echo
            echo === stats ===
            cat {input.stats_tsv} | datamash transpose | ( printf '%s\n' .TS 'box;LL.' ; cat ; printf '%s\n' .TE '.pl 0') | tbl | nroff
            echo
            echo === busco ===
            cat {input.busco_txt} | datamash transpose | ( printf '%s\n' .TS 'box;LL.' ; cat ; printf '%s\n' .TE '.pl 0') | tbl | nroff
            echo
        ) | tee {output}
        """


# ------------------------------------------------------------------------
# Check Git status
# ------------------------------------------------------------------------

rule run_git:
    input: DATA+"/summary-post-assembly.log"
    output: DATA+"/git-post-assembly.log"
    shell:
        """
	(
	    cd {PIPELINE}
	    echo
	    ( set -x ; git status )
	    echo
	    ( set -x ; git log -n1 )
	) 2>&1 | tee {output}
        """ 

# ------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------

rule all:
    input:
        DATA+"/summary-post-assembly.log",
        DATA+"/git-post-assembly.log"
    default_target: True

