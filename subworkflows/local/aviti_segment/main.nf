include { AVITIDISCOVERTILES     } from '../../../modules/local/aviti/discovertiles/main.nf'
include { AVITIMERGETILECHANNELS } from '../../../modules/local/aviti/mergetilechannels/main.nf'
include { AVITIWHOLECELLSEGMENT  } from '../../../modules/local/aviti/wholecellsegment/main.nf'
include { AVITINUCLEARSEGMENT    } from '../../../modules/local/aviti/nuclearsegment/main.nf'
include { AVITISTITCHWELL        } from '../../../modules/local/aviti/stitchwell/main.nf'
include { CELLMEASUREMENT         } from '../../../modules/local/cellmeasurement/main.nf'
include { SEGMENTATIONREPORT      } from '../../../modules/local/segmentationreport/main.nf'

workflow AVITI_SEGMENT {

    take:
    ch_aviti_samplesheet // channel: [ meta, run_dir ] -- meta.id is the sample name; meta.wells (optional) restricts --wells
    ch_nuclear_model     // channel: value, path to the staged Cellpose 3.x custom nuclear model

    main:

    ch_versions = channel.empty()

    //
    // Discover wells/tiles for each AVITI sample and write a flat per-tile
    // manifest (one row per (well, tile)), including the stage coordinates
    // stitching needs later.
    //
    AVITIDISCOVERTILES(
        ch_aviti_samplesheet
    )
    ch_versions = ch_versions.mix(AVITIDISCOVERTILES.out.versions.first())

    //
    // Fan the manifest out into one channel item per tile. This is the point
    // where per-sample AVITI processing becomes per-tile parallel: every tile
    // in every well of every sample becomes an independent task from here
    // until the stitching step groups them back by well.
    //
    // actin_tif is emitted as a literal "NO_FILE" placeholder when empty
    // (2-channel mode), since a Nextflow tuple element cannot be "no path" --
    // downstream modules detect this by name and omit --actin-tif.
    //
    AVITIDISCOVERTILES.out.manifest
        .splitCsv(elem: 1, header: true)
        .combine(ch_aviti_samplesheet, by: 0)
        .map { sample_meta, row, run_dir ->
            def base = file(run_dir).parent
            def meta_tile = [
                id           : "${sample_meta.id}__Well${row.well}__${row.tile}",
                sample       : sample_meta.id,
                well         : row.well,
                tile         : row.tile,
                x_mm         : row.x_mm,
                y_mm         : row.y_mm,
                channel_mode : row.channel_mode,
            ]
            [
                meta_tile,
                file("${base}/${row.nucleus_tif}"),
                file("${base}/${row.membrane_tif}"),
                row.actin_tif ? file("${base}/${row.actin_tif}") : file('NO_FILE'),
            ]
        }
        .set { ch_tiles }

    //
    // Merge the raw per-channel tiles into one named multi-channel OME-TIFF.
    // This carries un-normalised pixel values (unlike the Cellpose input
    // stacks built independently by the two segmentation modules below) --
    // it is the per-tile intensity image later stitched for
    // CELLMEASUREMENT/KRONOS2EMBEDDINGS.
    //
    AVITIMERGETILECHANNELS(
        ch_tiles
    )
    ch_versions = ch_versions.mix(AVITIMERGETILECHANNELS.out.versions.first())

    //
    // Whole-cell/membrane segmentation (Cellpose v4 SAM). Runs on the full
    // tile (no internal sub-tiling): AVITI tiles are already HPC-friendly in
    // size, unlike the whole-slide COMET images the SOPA path patchifies.
    //
    AVITIWHOLECELLSEGMENT(
        ch_tiles
    )
    ch_versions = ch_versions.mix(AVITIWHOLECELLSEGMENT.out.versions.first())

    //
    // Nuclear segmentation with the Cellpose 3.x custom model, in its own
    // environment/container -- Cellpose 3.x and 4.x cannot coexist in one.
    //
    AVITINUCLEARSEGMENT(
        ch_tiles.map { meta_tile, nucleus_tif, _membrane_tif, _actin_tif -> [ meta_tile, nucleus_tif ] },
        ch_nuclear_model
    )
    ch_versions = ch_versions.mix(AVITINUCLEARSEGMENT.out.versions.first())

    //
    // Group the three per-tile outputs back by well. Joining on meta_tile
    // (by: 0) is safe here because every branch above originates from the
    // same ch_tiles, so the same meta_tile map is used as the join key by
    // all three.
    //
    AVITIWHOLECELLSEGMENT.out.cell_mask
        .join(AVITINUCLEARSEGMENT.out.nuclear_mask, by: 0)
        .join(AVITIMERGETILECHANNELS.out.image, by: 0)
        .map { meta_tile, cell_mask, nuclear_mask, image ->
            def well_meta = [
                id           : "${meta_tile.sample}__Well${meta_tile.well}",
                sample       : meta_tile.sample,
                well         : meta_tile.well,
                channel_mode : meta_tile.channel_mode,
            ]
            def row = [
                tile        : meta_tile.tile,
                x_mm        : meta_tile.x_mm,
                y_mm        : meta_tile.y_mm,
                cell_mask   : cell_mask.name,
                nuclear_mask: nuclear_mask.name,
                image_tif   : image.name,
            ]
            [ well_meta, row, cell_mask, nuclear_mask, image ]
        }
        // Collects each well's tile rows and matching files into parallel
        // lists -- exactly the shape AVITISTITCHWELL expects.
        .groupTuple(by: 0)
        .set { ch_well_stitch_input }

    //
    // Stitch per-tile masks/images into per-well outputs, using the stage
    // coordinates carried in each tile row.
    //
    AVITISTITCHWELL(
        ch_well_stitch_input
    )
    ch_versions = ch_versions.mix(AVITISTITCHWELL.out.versions.first())

    //
    // Feed stitched per-well outputs into the existing, unmodified
    // CELLMEASUREMENT module. meta.id is well-scoped
    // (`<sample>__Well<well>`), keeping AVITI outputs distinct from any
    // COMET/MIBI sample sharing the sample name.
    //
    AVITISTITCHWELL.out.image
        .join(AVITISTITCHWELL.out.nuclear_mask, by: 0)
        .join(AVITISTITCHWELL.out.cell_mask, by: 0)
        .set { ch_cellmeasurement }

    CELLMEASUREMENT(
        ch_cellmeasurement
    )
    ch_versions = ch_versions.mix(CELLMEASUREMENT.out.versions.first())

    ch_annotations = CELLMEASUREMENT.out.annotations

    //
    // Assemble the KRONOS input channel: stitched well image + whole-cell
    // mask. KRONOS itself is invoked once at the top level
    // (workflows/sp_segment.nf), same as every other segmenter.
    //
    AVITISTITCHWELL.out.image
        .join(AVITISTITCHWELL.out.cell_mask, by: 0)
        .set { ch_kronos_input }

    //
    // Optional SEGMENTATIONREPORT module. AVITI's channel names are fixed
    // ("Nucleus", "Cell-Membrane"[, "Actin"]), unlike the free-form
    // per-sample channel mapping the COMET/MIBI samplesheet supplies.
    // run_cellpose is reported as true since both AVITI segmenters are
    // Cellpose-family models; run_mesmer/run_cellsam are false.
    //
    ch_report = channel.empty()
    if (params.generate_report) {
        AVITISTITCHWELL.out.image
            .join(ch_annotations, by: 0)
            .map { meta, image, annotations ->
                [
                    meta,
                    annotations,
                    false, // run_mesmer
                    true,  // run_cellpose
                    false, // run_cellsam
                    'Nucleus',
                    meta.channel_mode == '2ch' ? 'Cell-Membrane' : 'Cell-Membrane:Actin',
                    image,
                ]
            }
            .set { ch_segmentation_report }

        SEGMENTATIONREPORT(
            ch_segmentation_report
        )
        ch_versions = ch_versions.mix(SEGMENTATIONREPORT.out.versions.first())
        ch_report = SEGMENTATIONREPORT.out.report
    }

    emit:
    nuclear_segmentation_mask   = AVITISTITCHWELL.out.nuclear_mask // channel: [ val(meta), *.tif ]
    wholecell_segmentation_mask = AVITISTITCHWELL.out.cell_mask    // channel: [ val(meta), *.tif ]
    annotations                 = ch_annotations                    // channel: [ val(meta), *.geojson ]
    kronos_input                = ch_kronos_input                   // channel: [ val(meta), tiff, whole_cell_mask ]
    report                      = ch_report                         // channel: [ val(meta), *.html ]

    versions = ch_versions                                          // channel: [ versions.yml ]
}
