#!/usr/bin/env nextflow

/*
========================================================================================
   Variant calling Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-dock
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/

def pipeline_name = """\
         S C B I R   D O C K I N G   P I P E L I N E
         ===========================================
         """.stripIndent()

if ( params.help ) {
   println """${pipeline_name}
         Nextflow pipeline to ....

         Usage:
            nextflow run scbirlab/nf-template --sample_sheet <csv> --inputs <dir>
            nextflow run scbirlab/nf-template -c <config-file>

         Required parameters:
            sample_sheet      Path to a CSV with information about the samples 
                                 to be processed

         Optional parameters (with defaults):  
            inputs             Directory containing inputs. Default: "./inputs".

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to sample_sheet")
}

log.info """${pipeline_name}
         inputs
            input dir.     : ${params.inputs}
            sample sheet   : ${params.sample_sheet}
         output            : ${params.outputs}
         """
         .stripIndent()

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/
include {
   DETECT_POCKET;
   SplitPockets
} from './modules/fpocket.nf'
include {
    SplitLigands;
    PREPARE_LIGANDS
} from './modules/ligands.nf'
include {
   GNINA_DOCK;
   Extract_Gnina_scores;
   AGGREGATE_SCORES;
} from './modules/gnina.nf'
include {
   FETCH_STRUCTURES
} from './modules/pdb.nf'
include {
    COMPUTE_FEATURES
} from './modules/features.nf'
include {
    TRAIN_SURROGATE;
    PREDICT_FULL_MATRIX
} from './modules/surrogate.nf'
include {
    VALIDATE;
    ASSIGN_TARGETS
} from './modules/annotate.nf'

workflow {

    // ── Parse protein metadata ──
    Channel.fromPath( 
        "${params.sample_sheet}", 
        checkIfExists: true,
    )
        .splitCsv(header: true)
        .map { row -> tuple( row.protein_id, row.pdb_id, row.chain ? row.chain : "A", row.uniprot_id ) }
        .set { proteins_ch }

    // ── Phase 1: Prepare ──
    ( params.test ? proteins_ch.take( 5 ) : proteins_ch ) 
        | FETCH_STRUCTURES
        | DETECT_POCKET

    SplitLigands(
        Channel.fromPath( 
            "${params.inputs}/${params.ligands_sdf}", 
            checkIfExists: true,
        ),
        Channel.value( params.test ? 2 : 10 ),
    )
    SplitLigands.out
        .flatten()
        .set { all_ligands }

    ( params.test ? all_ligands.take( 2 ) : all_ligands )
        | PREPARE_LIGANDS

    // ── Phase 2: Dock (all ligand chunks × all proteins) ──
    DETECT_POCKET.out.main
        | SplitPockets
    SplitPockets.out
        .transpose()
        .combine( 
            PREPARE_LIGANDS.out
        )
        .set { docking_inputs }

    ( params.test ? docking_inputs.take( 3 ) : docking_inputs )
        | GNINA_DOCK

    GNINA_DOCK.out.main
        | Extract_Gnina_scores

    // Collect all score files and aggregate
    Extract_Gnina_scores.out
        .map { v -> v[-1] }
        .collect() 
        | AGGREGATE_SCORES

    if ( params.surrogate_type ) {
        // ── Phase 3: Train surrogate ──
        // Collect protein sequences for ESM embeddings
        FETCH_STRUCTURES.out
            .map { v -> v[2] }
            .collect()
            .set { protein_structures }
            // In practice, you'd extract sequences from PDBs or use a pre-made FASTA

        COMPUTE_FEATURES(
            Channel.fromPath( 
                "${params.inputs}/${params.ligands_sdf}", 
                checkIfExists: true,
            ),
            protein_structures,
        )

        TRAIN_SURROGATE(
            AGGREGATE_SCORES.out.scores_long,
            COMPUTE_FEATURES.out.compound_fps,
            COMPUTE_FEATURES.out.protein_embeds,
        )

        // ── Phase 4: Predict full matrix ──
        PREDICT_FULL_MATRIX(
            TRAIN_SURROGATE.out.model,
            COMPUTE_FEATURES.out.compound_fps,
            COMPUTE_FEATURES.out.protein_embeds,
        )

        if ( params.ortholog_table ) {

            // ── Phase 5: Annotate + validate ──
            ASSIGN_TARGETS(
                PREDICT_FULL_MATRIX.out.full_matrix,
                Channel.fromPath( params.ortholog_table, checkIfExists: true ),
            )

        }

        // Optional: validate if known activity data provided
        if ( params.known_activities ) {
            VALIDATE(
                ASSIGN_TARGETS.out.selectivity,
                Channel.fromPath( params.known_activities, checkIfExists: true ),
            )
        }
    }
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
