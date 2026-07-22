# Config for ML-probability extraction: samples 07–10 × genome alignment.
# Samples 07–10 are not transcriptome-aligned, so there's no transcriptome variant.

config <- list(
    basecaller   = "m6A_pseU_2OmeU",
    align        = "genome",

    bam_dir      = "/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/03_dorado_aligned_genome",
    bam_filename = function(sid) paste0(sid, "_genome_primary.bam"),

    sample_ids   = paste0("sample0", 7:10),
    conditions   = c("Treated", "Control", "Treated", "Control"),

    mod_names    = c("a" = "m6A", "17802" = "pseU", "19227" = "2OmeU"),
    mod_levels   = c("m6A", "pseU", "2OmeU")
)
