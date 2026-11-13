#!/usr/bin/env bash
#SBATCH --job-name=LEAP
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=168:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
# Uncomment and edit if required by your cluster:
##SBATCH --partition=YOUR_PARTITION
##SBATCH --account=YOUR_ACCOUNT

set -euo pipefail

# -----------------------------------------------------------------------------
# EDIT: pipeline and input/output paths
# Relative paths below are resolved from LEAP_DIR, including the reference paths.
# Keep the database credentials inside the JSON; do not add them to this script.
# -----------------------------------------------------------------------------
LEAP_DIR="/path/to/ensembl-LEAP"
DB_SECRETS="/path/to/db_secrets.json"
READS_BAM="/path/to/reads.bam"
READS_INDEX="${READS_BAM}.csi"    # Change to the actual .csi or .bai if different.
CAGE_BED="/path/to/cage.bed"
POLYA_BED="/path/to/polya.bed"
WORK_DIR="/path/to/leap_work"
OUTPUT_DIR="/path/to/leap_results"

REFERENCE_FASTA="../reference_genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"

# -----------------------------------------------------------------------------
# EDIT: LEAP settings (retain your previous run's values where appropriate)
# These Slurm resources are for the Nextflow launcher. LEAP's hpc_slurm profile
# controls resources for the individual processing jobs.
# -----------------------------------------------------------------------------
PROFILE="hpc_slurm"
MINIMUM_READ_SUPPORT=2
OVERLAP_PENALTY=10
EXTEND_ON_CURRENT_ANNOTATION=true
RUN_DATACHECKS=true
PREPARE_MANUAL_REVIEW=true
# Every run also writes OUTPUT_DIR/impact/index.html and CSV/JSON statistics,
# including when manual review is disabled. No extra reporting parameter is needed.
RESUME=false

# -----------------------------------------------------------------------------
# EDIT: manual-review / HAVANA JBrowser settings
# FTP settings can stay blank for offline review, or be filled in later in the
# generated public launch script. The reference above is used for both modes.
# Manual review also prepares positive/negative/unsplit mapped-input BAMs with
# CSI indexes from READS_BAM. Allow space for these views and temporary copies.
# -----------------------------------------------------------------------------
REVIEW_GENOME_ID="hg38"
REVIEW_SPECIES="homo_sapiens"
JBROWSER_PIPELINE="/path/to/havana-jbrowser"
FTP_PROJECT="LEAP_review"
PUBLIC_BASE_URL=""               # e.g. https://ftp.example.org/root
FTP_TARGET_DIR=""                # FTP filesystem root, not the project folder.
HUB_EMAIL=""

# Load your site's Nextflow/Java/container modules here if needed.
# module load ...

# -----------------------------------------------------------------------------
# Launch: normally no edits needed below this point.
# Submit with: sbatch /path/to/run_leap.sbatch
# -----------------------------------------------------------------------------
[[ -f "$LEAP_DIR/main.nf" ]] || {
    echo 'Set LEAP_DIR to the LEAP checkout before submitting.' >&2
    exit 1
}
cd -- "$LEAP_DIR"
LEAP_DIR="$(pwd -P)"
command -v nextflow >/dev/null || {
    echo 'Load or install Nextflow before running this script.' >&2
    exit 1
}

require_file() {
    [[ -f "$2" && -r "$2" ]] || {
        printf 'Set %s to a readable file: %s\n' "$1" "$2" >&2
        exit 1
    }
}
require_file DB_SECRETS "$DB_SECRETS"
require_file READS_BAM "$READS_BAM"
require_file READS_INDEX "$READS_INDEX"
require_file CAGE_BED "$CAGE_BED"
require_file POLYA_BED "$POLYA_BED"
if [[ "$RUN_DATACHECKS" == true ]]; then
    require_file REFERENCE_FAI "$REFERENCE_FAI"
fi
if [[ "$PREPARE_MANUAL_REVIEW" == true ]]; then
    require_file REFERENCE_FASTA "$REFERENCE_FASTA"
fi

command=(nextflow run "$LEAP_DIR/main.nf"
    -profile "$PROFILE"
    -work-dir "$WORK_DIR"
    --db_secrets "$DB_SECRETS"
    --reads "$READS_BAM"
    --reads_index "$READS_INDEX"
    --cage "$CAGE_BED"
    --polya "$POLYA_BED"
    --minimum_read_support "$MINIMUM_READ_SUPPORT"
    --overlap_penalty "$OVERLAP_PENALTY"
    --extend_on_current_annotation "$EXTEND_ON_CURRENT_ANNOTATION"
    --output_directory "$OUTPUT_DIR"
    --run_datachecks "$RUN_DATACHECKS"
    --datacheck_reference_lengths "$REFERENCE_FAI"
    --prepare_manual_review "$PREPARE_MANUAL_REVIEW"
    --review_reference "$REFERENCE_FASTA"
    --review_genome_id "$REVIEW_GENOME_ID"
    --review_species "$REVIEW_SPECIES"
    --review_jbrowser_pipeline "$JBROWSER_PIPELINE"
    --review_ftp_project "$FTP_PROJECT"
    --review_public_base_url "$PUBLIC_BASE_URL"
    --review_ftp_target_dir "$FTP_TARGET_DIR"
    --review_hub_email "$HUB_EMAIL")

[[ "$RESUME" != true ]] || command+=(-resume)
"${command[@]}" "$@"

if [[ "$PREPARE_MANUAL_REVIEW" == true ]]; then
    printf 'Manual-review bundle: %s/manual_review/\n' "$OUTPUT_DIR"
    echo 'Read evidence: reads/{positive,negative,unsplit}.bam and their CSI indexes.'
    echo 'Choose run_jbrowser.sbatch (public FTP) or run_jbrowser_offline.sbatch.'
fi
