process FETCH_STRUCTURES {
    /*
     * For each protein: fetch PDB structure (or AF2 model),
     * extract the relevant chain, remove waters/ligands.
     */
    tag "${id}:PDB=${pdb_id}:UniProt=${uniprot_id}"

    publishDir(
        "${params.outputs}/receptors", 
        mode: 'copy', 
        saveAs: { v -> "${id}-${uniprot_id}-${v}"},
    )

    input:
    tuple val( id ), val( pdb_id ), val( chain ), val( uniprot_id )

    output:
    tuple val( id ), val( uniprot_id ), path( "prepped.pdb" )

    script:
    if ( pdb_id ) {
        """
        set -euox pipefail
        # Fetch structure
        if [[ "${pdb_id}" == AF2_* ]]; then
            # AlphaFold DB
            af_id=\$(echo ${pdb_id} | sed 's/AF2_//')
            curl "https://alphafold.ebi.ac.uk/files/AF-\${af_id}-F1-model_v6.pdb" \
                -o prepped0.pdb
        else
            pdb_fetch "${pdb_id}" > fetched.pdb
            pdb_selchain -${chain} fetched.pdb > chain.pdb
            # Clean: remove waters, heteroatoms; add hydrogens; fix missing residues
            pdbfixer chain.pdb \
                --verbose \
                --replace-nonstandard \
                --add-residues \
                --output="prepped0.pdb"
        fi

        pdbfixer prepped0.pdb \
            --verbose \
            --add-atoms=all \
            --output="prepped00.pdb"

        grep -v '^MODEL\\|^ENDMDL' "prepped00.pdb" > prepped.pdb

        """
    } 
    else {
        """
        set -euox pipefail

        # AlphaFold DB
        curl "https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-model_v6.pdb" \
            -o prepped0.pdb

        pdbfixer prepped0.pdb \
            --verbose \
            --add-atoms=all \
            --output="prepped.pdb"

        grep -v '^MODEL\\|^ENDMDL' "prepped00.pdb" > prepped.pdb

        """
    }
    
}
