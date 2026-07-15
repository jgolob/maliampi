nextflow.enable.dsl=2

include { ValidateRefpkg } from '../../modules/refpkg_validate'

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

workflow validate_refpkg_entry {
    if (params.refpkg == null) {
        error 'Missing required --refpkg'
    }
    validate_refpkg(channel.fromPath(params.refpkg, checkIfExists: true))
}
