/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap       } from 'plugin/nf-schema'
include { ANALYSE                } from '../subworkflows/local/analyse'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialmibi_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALMIBI {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Construct channel for only ANALYSE subworkflow
    //
    // TODO: can we somehow preserve key-value pairs in the samplesheet channel?
    //       Then we could manage these args in a less cumbersome way
    ch_samplesheet.filter {
        it[1] == ['analyse']
    }.map {
        sample,
        function,
        an_expression_file,
        an_hierarchy_file,
        an_sample_id,
        an_marker_column,
        an_markers,
        an_are_markers_split,
        an_cell_types,
        an_parent_types,
        an_metadata_cols,
        an_plot_heatmaps,
        an_plot_props,
        an_plot_umap,
        an_plot_clusters,
        an_plot_spatial,
        an_save_rdata,
        seg_mibi_tiff,
        seg_nuclear_channel,
        seg_membrane_channels,
        seg_combine_method,
        seg_level,
        seg_interior_threshold,
        seg_maxima_smooth,
        seg_min_nuclei_area,
        seg_remove_border_cells,
        seg_include_measurements,
        seg_pixel_expansion,
        seg_padding -> [
            sample,
            an_expression_file,
            an_hierarchy_file,
            an_sample_id,
            an_marker_column,
            an_markers,
            an_are_markers_split,
            an_cell_types,
            an_parent_types,
            an_metadata_cols,
            an_plot_heatmaps,
            an_plot_props,
            an_plot_umap,
            an_plot_clusters,
            an_plot_spatial,
            an_save_rdata
        ]
    }.set { ch_analyse_samplesheet }

    //
    // Run the main ANALYSE subworkflow
    //
    ANALYSE(
        ch_analyse_samplesheet
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  ''  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
