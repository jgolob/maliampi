nextflow.enable.dsl=2

params.help = false
params.dedup_jplace = null
params.sv_h5ad = null
params.placement_identity = null
params.dataset_id = null
params.phylotype_thresholds = [0.1, 0.3, 1.0]
params.phylotype_kr_thresholds = []
params.phylotype_add_thresholds = []
params.phylotype_distance = 'legacy'
params.phylotype_lwr_overlap = 0.1
params.output = '.'

include { Phylotypes; ImportPhylotypeMapping; AggregatePhylotype; BuildPhylotypeSet } from '../../modules/phylotypes'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / phylotypes
    ─────────────────────────────────────
    Cluster placed sequence variants into phylotypes at one or more
    distance thresholds.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/phylotypes.nf [options]

    Required:
      --dedup_jplace          Deduplicated jplace file from placement
      --sv_h5ad               H5AD file with sequence variants
      --placement_identity    Placement identity JSON from placement step
      --dataset_id            Dataset identifier

    Options:
      --phylotype_thresholds      Distance thresholds (default: [0.1, 0.3, 1.0])
      --phylotype_add_thresholds  Additional thresholds to merge in (default: [])
      --phylotype_distance        Distance metric: legacy (default) or kr
      --phylotype_kr_thresholds   KR-calibrated thresholds (required when distance=kr)
      --phylotype_lwr_overlap     LWR overlap cutoff (default: 0.1)
      --output                    Output directory (default: .)
      --help                      Show this help message
    """.stripIndent()
}

def effective_phylotype_thresholds(base, additions, distance, kr_base) {
    if (!(distance in ['legacy', 'kr'])) {
        error "Unsupported --phylotype_distance ${distance}; use legacy or kr"
    }
    def selected = distance == 'kr' ? kr_base : base
    if (distance == 'kr' && !selected) {
        error '--phylotype_distance kr requires an explicit calibrated --phylotype_kr_thresholds list'
    }
    def values = (selected + additions).collect { value ->
        try {
            def number = value as BigDecimal
            if (number <= 0 || !Double.isFinite(number.doubleValue())) {
                error "Phylotype threshold must be positive and finite: ${value}"
            }
            number
        } catch (NumberFormatException _ignored) {
            error "Invalid phylotype threshold: ${value}"
        }
    }
    return values.unique().sort()
}

workflow phylotypes {
    take:
    dedup_jplace
    sv_h5ad
    placement_identity
    dataset_id
    thresholds
    distance
    lwr_overlap

    main:
    dedup_jplace.combine(channel.fromList(thresholds))
        .map { jplace, threshold -> tuple(threshold, threshold.toString().replace('.', 'p'), jplace) }
        .set { phylotype_inputs }
    Phylotypes(phylotype_inputs, distance, lwr_overlap)
    ImportPhylotypeMapping(Phylotypes.out.external_mapping)
    ImportPhylotypeMapping.out.mapping.combine(sv_h5ad).combine(placement_identity).set { phylotype_table_inputs }
    AggregatePhylotype(phylotype_table_inputs)
    ImportPhylotypeMapping.out.mapping.map { _threshold, _token, mapping -> mapping }.collect().set { all_mappings }
    BuildPhylotypeSet(dedup_jplace, placement_identity, all_mappings, dataset_id, thresholds, distance, lwr_overlap)

    emit:
    mappings = ImportPhylotypeMapping.out.mapping
    abundance_h5ad = AggregatePhylotype.out.abundance_h5ad
    set_bundle = BuildPhylotypeSet.out.set_bundle
}

workflow {
    def required = [
        '--dedup_jplace': params.dedup_jplace,
        '--sv_h5ad': params.sv_h5ad,
        '--placement_identity': params.placement_identity,
        '--dataset_id': params.dataset_id,
    ]
    def missing = required.findAll { _k, v -> v == null }.collect { k, _v -> k }

    if (params.help || missing) {
        helpMessage()
        if (!params.help) {
            error "Missing required parameter(s): ${missing.join(', ')}"
        }
        return
    }
    def thresholds = effective_phylotype_thresholds(
        params.phylotype_thresholds,
        params.phylotype_add_thresholds,
        params.phylotype_distance,
        params.phylotype_kr_thresholds,
    )
    phylotypes(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(params.sv_h5ad, checkIfExists: true),
        channel.fromPath(params.placement_identity, checkIfExists: true),
        params.dataset_id,
        thresholds,
        params.phylotype_distance,
        params.phylotype_lwr_overlap,
    )
}
