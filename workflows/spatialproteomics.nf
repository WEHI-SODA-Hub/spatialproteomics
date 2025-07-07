/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap       } from 'plugin/nf-schema'
include { ANALYSE                } from '../subworkflows/local/analyse'
include { BACKGROUNDSUBTRACT     } from '../subworkflows/local/backgroundsubtract'
include { BACKSUBMESMER          } from '../subworkflows/local/backsubmesmer'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialproteomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALPROTEOMICS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Construct channel for background subtraction/segmentation workflow
    //
    // TODO: can we somehow preserve key-value pairs in the samplesheet channel?
    //       Then we could manage these args in a less cumbersome way
    ch_samplesheet.map {
        sample,
        run_backsub,
        run_mesmer,
        run_analyse,
        an_expression_file,
        an_hierarchy_file,
        an_sample_id,
        an_marker_column,
        an_markers,
        an_are_markers_split,
        an_cell_types,
        an_parent_types,
        an_metadata_cols,
        an_plot_metas,
        an_plot_heatmaps,
        an_plot_props,
        an_plot_umap,
        an_plot_clusters,
        an_plot_spatial,
        an_save_rdata,
        seg_tiff,
        seg_nuclear_channel,
        seg_membrane_channels,
        seg_combine_method,
        seg_level,
        seg_maxima_threshold,
        seg_interior_threshold,
        seg_maxima_smooth,
        seg_min_nuclei_area,
        seg_remove_border_cells,
        seg_include_measurements,
        seg_pixel_expansion,
        seg_padding,
        seg_skip_measurements -> [
            sample,
            run_backsub,
            run_mesmer,
            seg_tiff,
            seg_nuclear_channel,
            seg_membrane_channels,
            seg_combine_method,
            seg_level,
            seg_maxima_threshold,
            seg_interior_threshold,
            seg_maxima_smooth,
            seg_min_nuclei_area,
            seg_remove_border_cells,
            seg_pixel_expansion,
            seg_padding,
            seg_skip_measurements
        ]
    }.branch { it ->
        backsub_only: it[1].contains(true) && !it[2].contains(true)
        backsub_mesmer: it[1].contains(true) && it[2].contains(true)
        mesmer_only: !it[1].contains(true) && it[2].contains(true)
    }.set { ch_segmentation_samplesheet }

    // TODO: setup test data compatible with background subtraction
    //
    // Run the BACKGROUNDSUBTRACT subworkflow for samples that ONLY require
    // background subtraction (no segmentation)
    //
    BACKGROUNDSUBTRACT(
        ch_segmentation_samplesheet.backsub_only.map {
            sample,
            run_backsub,
            run_mesmer,
            seg_tiff,
            seg_nuclear_channel,
            seg_membrane_channels,
            seg_combine_method,
            seg_level,
            seg_maxima_threshold,
            seg_interior_threshold,
            seg_maxima_smooth,
            seg_min_nuclei_area,
            seg_remove_border_cells,
            seg_pixel_expansion,
            seg_padding,
            seg_skip_measurements -> [
                sample,
                seg_tiff,
            ]
        }
    )

    //
    // Run the BACKSUBMESMER subworkflow for samples that require
    // background subtraction and mesmer segmentation
    //
    BACKSUBMESMER(
        ch_segmentation_samplesheet.backsub_mesmer
    )

    //
    // Construct channel for only ANALYSE subworkflow
    //
    ch_samplesheet.filter {
        it[3].contains(true) // run_analyse true for sample
    }.map {
        sample,
        run_backsub,
        run_mesmer,
        run_analyse,
        an_expression_file,
        an_hierarchy_file,
        an_sample_id,
        an_marker_column,
        an_markers,
        an_are_markers_split,
        an_cell_types,
        an_parent_types,
        an_metadata_cols,
        an_plot_metas,
        an_plot_heatmaps,
        an_plot_props,
        an_plot_umap,
        an_plot_clusters,
        an_plot_spatial,
        an_save_rdata,
        seg_tiff,
        seg_nuclear_channel,
        seg_membrane_channels,
        seg_combine_method,
        seg_level,
        seg_maxima_threshold,
        seg_interior_threshold,
        seg_maxima_smooth,
        seg_min_nuclei_area,
        seg_remove_border_cells,
        seg_include_measurements,
        seg_pixel_expansion,
        seg_padding,
        seg_skip_measurements -> [
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
            an_plot_metas,
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
    // TODO: fix sample_id column issue with blank values (will pass a blank
    // list as sample_id if it is null to quarto notebook, which will crash)
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
