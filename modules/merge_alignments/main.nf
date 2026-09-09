process MERGE_ALIGNMENTS {
    tag "merge"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/7920d351ee6307611f471013dcddff6d8a5c7ec7bf723d854e99a545376b63d8/data' :
        'community.wave.seqera.io/library/python_numpy_faiss-cpu_mkl_libblas:078dd4f35c795d6e' }"

    input:
    path alignments, stageAs: 'aln_tsv/*'
    path texts,      stageAs: 'aln_txt/*'
    path evd
    val split_for_plots

    output:
    path "alignments.tsv", emit: alignments
    path "alignments.txt", emit: text
    path "by_query/*.tsv", emit: plot_alignments, optional: true
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def plot_script = split_for_plots && !(alignments instanceof Collection && alignments.isEmpty()) ? """
    mkdir -p by_query
    split_alignments.py \\
        --input aln_tsv/*.tsv \\
        --outdir by_query
    """ : ''
    def tsv_arg = alignments instanceof Collection && alignments.isEmpty() ? '' : '--alignments aln_tsv/*'
    def txt_arg = texts instanceof Collection && texts.isEmpty() ? '' : '--alignment-text aln_txt/*'
    """
    export PYTHONPATH="${projectDir}/bin"
    merge_alignments.py \\
        ${tsv_arg} \\
        ${txt_arg} \\
        --evd ${evd} \\
        --output alignments.tsv \\
        --alignment-text-output alignments.txt \\
        ${args}

    ${plot_script}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """

    stub:
    def plot_stub = split_for_plots ? """
    mkdir -p by_query
    echo -e "cluster_id\\tquery_id\\ttarget_id" > by_query/query.tsv
    """ : ''
    """
    echo -e "cluster_id\\tquery_id\\ttarget_id\\tscore\\ttotal_score\\tmax_score\\tbit_score\\tevalue\\tevalue_pair\\tlog_evalue\\talignment_count\\tcluster_ids\\thsp_scores\\thsp_spans\\tquery_start\\tquery_end\\ttarget_start\\ttarget_end\\tquery_length\\ttarget_length\\tmatch_count\\taligned_columns\\tungapped_columns\\tbase_matches\\tstructure_matches\\tseed_count\\tmax_seed_score\\tquery_sequence\\tquery_structure\\ttarget_sequence\\ttarget_structure\\tquery_aligned\\ttarget_aligned" > alignments.tsv
    touch alignments.txt

    ${plot_stub}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import sys; print('.'.join(map(str, sys.version_info[:3])))")
    END_VERSIONS
    """
}
