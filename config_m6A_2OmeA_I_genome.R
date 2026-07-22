# Config for ML-probability extraction: samples 01–04 × genome alignment.

config <- list(
    basecaller   = "m6A_2OmeA_I",
    align        = "genome",

    bam_dir      = "/dcs04/hicks/data/sparthib/drna_retinalrs/processed_data/03_dorado_aligned_genome",
    bam_filename = function(sid) paste0(sid, "_genome_primary.bam"),

    sample_ids   = paste0("sample0", 1:4),
    conditions   = c("Treated", "Control", "Treated", "Control"),

    mod_names    = c("a" = "m6A", "17596" = "Inosine", "69426" = "2OmeA"),
    mod_levels   = c("m6A", "Inosine", "2OmeA")
)
