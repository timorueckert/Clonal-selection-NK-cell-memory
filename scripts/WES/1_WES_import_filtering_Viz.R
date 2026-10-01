# Import, filtering and visualization of mutations from WES data
# pre-processed with Imprint pipeline

### Unfiltered VCF files are only available in closed GHGA repository
### due to GDPR regulations; For re-analysis, utilize filtered 
### mutation list exported at the end of this script

#### Import and filtering of data processed with Ollie's IMPRINT pipeline ####

# Configure inputs
wes_dir <- path.expand("~/share/wes")

# Helper functions from imprint
source(file.path(wes_dir, "imprint/imprint/R/wes_helper_functions.R"))

bulk_dir <- file.path(wes_dir, "imprint/healthy_non_umi_outs")
reference_dir <- file.path(wes_dir, "reference_files")
blacklist_dir <- file.path(reference_dir, "blacklists")

cohorts <- c("cohort_1_cegat", "cohort_2_novogene", "cohort_4_cegat")
excluded_donors <- c("CH020") # Duplicated donor

# filtering thresholds
max_gnomad_af <- 1e-4 # above this a site is most likely germline, not somatic
min_read_position_pv <- 0.05 # vafator read-position rank-sum test
max_indel_length <- 10 # vafator counts alt reads by pileup; long events score badly
target_fdr <- 0.01 # varlociraptor filter-calls control-fdr --mode local-smart
min_reportable_vaf <- 0.025  # burden counting only, not for dN/dS

# load variants

manifests <- read_output_manifests(setNames(file.path(bulk_dir, cohorts), cohorts))
vcf_meta  <- build_manifest_vcf_meta(manifests, kind = "ensemble_vcf")

variants <- parallel::mclapply(
  seq_len(nrow(vcf_meta)),
  function(i) parse_somatic_vcf(vcf_meta$file[i], vcf_meta$tumor_id[i], vcf_meta[i, ], chip_genes = chip_genes),
  mc.cores = 16) %>%
  bind_rows() %>%
  left_join(
    vcf_meta %>% select(cohort, donor_id, tumor_id, pair_id),
    by = c("cohort", "donor_id", "tumor_id"), relationship = "many-to-one") %>%
  filter(!donor_id %in% excluded_donors) %>%
  mutate(
    variant_id = paste(cohort, donor_id, cell_type, chrom, pos, ref, alt, sep = ":"),
    indel_length = abs(nchar(ref) - nchar(alt)))

# cross-sample artefact evidence

# PoN is kit-scoped: novogene has its own, both cegat cohorts share
pon_cegat <- read_pon_vcf(file.path(bulk_dir, "pon_sweep/combined_pon_min2.vcf.gz"))

pon_combined <- bind_rows(
  read_pon_vcf(file.path(bulk_dir, "cohort_2_novogene/cohort/pon/pon.vcf.gz")) %>%
    mutate(cohort = "cohort_2_novogene"),
  pon_cegat %>% mutate(cohort = "cohort_1_cegat"),
  pon_cegat %>% mutate(cohort = "cohort_4_cegat")) %>%
  distinct(chrom, pos, ref, alt, cohort, .keep_all = TRUE)

variants <- mark_pon(variants, pon_combined)

# annotation only
difficult_region_beds <- c(in_blacklist = "ENCFF256LFZ.bed", in_segdup = "segdup_hg38.bed", in_repeat = "rmsk_hg38.bed")

variant_gr <- GRanges(
  seqnames = variants$chrom,
  ranges   = IRanges(start = variants$pos, width = pmax(nchar(variants$ref), 1L)))

for (region_col in names(difficult_region_beds)) {
  region_gr <- rtracklayer::import(file.path(blacklist_dir, difficult_region_beds[[region_col]]))
  variants[[region_col]] <- overlapsAny(variant_gr, region_gr, ignore.strand = TRUE)
}

# Recurrence counts DISTINCT DONORS carrying the site, not rows: one donor's two
# cell types are one observation of the site, not two
donors_per_site <- variants %>%
  distinct(cohort, donor_id, chrom, pos, ref, alt) %>%
  count(cohort, chrom, pos, ref, alt, name = "n_donors_with_variant_in_cohort")

variants <- variants %>%
  left_join(donors_per_site, by = c("cohort", "chrom", "pos", "ref", "alt")) %>%
  mutate(
    is_recurrent_in_cohort = n_donors_with_variant_in_cohort > 1,
    gnomad_max_af = pmax(
      coalesce(suppressWarnings(as.numeric(gnomADv4e_AF)), -Inf),
      coalesce(suppressWarnings(as.numeric(gnomADv4g_AF)), -Inf)),
    gnomad_max_af = if_else(is.infinite(gnomad_max_af), NA_real_, gnomad_max_af))

# candidate rule
variants <- variants %>%
  mutate(
    statistical_candidate = coalesce(
        # cross-sample artefactual evidence
        !in_pon &
        # recurrent sites must be robustly evidenced
        !(is_recurrent_in_cohort & n_callers == 1) &
        # significant germline evidence
        (is.na(gnomad_max_af) | gnomad_max_af <= max_gnomad_af) &
        # read-position bias testing from vafator
        (is.na(tumor_rspos_pv) | tumor_rspos_pv >= min_read_position_pv) &
        # every call needs real alt-read support; indels additionally need to be short
        (VARIANT_CLASS == "SNV" | indel_length <= max_indel_length),
      FALSE))

table(variants$statistical_candidate, variants$is_recurrent_in_cohort, useNA = 'ifany')

# Control FDR rate of considered candidates with varlociraptor
# Rescue variants called by more than one caller
primary_calls <- variants %>%
  filter(statistical_candidate) %>%
  local_smart_fdr(fdr = target_fdr) %>%
  filter(primary_somatic | n_callers > 1) %>%
  arrange(cohort, donor_id, cell_type, chrom, pos, ref, alt) %>%
  reorder_imprint_cols()


# Filter
muts_filtered <- primary_calls %>%
  filter(cell_type %in% c("NKG2C", "NKG2A")) %>% 
  filter(tumor_af > 0) %>% 
  filter(tumor_alt_reads > 2)

# Mark duplicated variants
muts_filtered <- muts_filtered %>% dplyr::group_by(cell_type, cohort, pos) %>%
  dplyr::mutate(duplicated = n() > 1) %>%
  ungroup()

table(muts_filtered$duplicated)

# Remove duplicated variants
muts_filtered <- muts_filtered %>%
  filter(!duplicated)

# Export
write.csv(muts_filtered, "data/WES/outs/2609_full_cohort_filtered.csv")