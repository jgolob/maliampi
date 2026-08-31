params.help = false
params.refpkg = null
params.output = '.'

nextflow.enable.dsl=2

include { ValidateRefpkg } from '../../modules/refpkg_validate'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / validate_refpkg
    ─────────────────────────────────────
    Validate a reference package archive.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/validate_refpkg.nf [options]

    Required:
      --refpkg        Reference package (.tar.gz)

    Options:
      --output        Output directory (default: .)
      --help          Show this help message
    """.stripIndent()
}

workflow validate_refpkg {
    take:
    refpkg

    main:
    ValidateRefpkg(refpkg)

    emit:
    refpkg = ValidateRefpkg.out.refpkg
    identity = ValidateRefpkg.out.identity
    report = ValidateRefpkg.out.report
}

workflow {
    if (params.help || params.refpkg == null) {
        helpMessage()
        if (!params.help) {
            error "Missing required parameter(s): --refpkg"
        }
        return
    }
    validate_refpkg(channel.fromPath(params.refpkg, checkIfExists: true))
}
