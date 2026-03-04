process FETCH_STRUCTURES {
    /*
     * For each protein: fetch PDB structure (or AF2 model),
     * extract the relevant chain, remove waters/ligands.
     */
    tag "${id}"

    input:
    tuple val( id ), val( pdb_id ), val( chain ), val( uniprot_id )

    output:
    tuple val( id ), val( uniprot_id ), path( "*_prepped.pdb" )

    script:
    """
    set -euox pipefail
    # Fetch structure
    if [[ "${pdb_id}" == "${pdb_id}" ]]; then #AF2_* ]]; then
        # AlphaFold DB
        af_id=\$(echo ${pdb_id} | sed 's/AF2_//')
        curl "https://alphafold.ebi.ac.uk/files/AF-${uniprot_id}-F1-model_v6.pdb" \
            -o prepped.pdb
    else
        pdb_fetch ${pdb_id} > complete.pdb
        pdb_selchain -${chain} complete.pdb > raw.pdb
        # Clean: remove waters, heteroatoms; add hydrogens; fix missing residues
        #pdbfixer raw.pdb \
        #    --verbose \
        #    --replace-nonstandard \
        #    --add-residues \
        #    --output="prepped.pdb"
        cp raw.pdb prepped.pdb
    fi

    #pdbfixer raw.pdb \
    #    --verbose \
    #    --add-atoms=all \
    #    --output="${id}_prepped.pdb"
    cp prepped.pdb "${id}_prepped.pdb"

    """
}
