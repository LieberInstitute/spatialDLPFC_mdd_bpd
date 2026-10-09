#   Gather sample info to form two tibbles:
#
#   - One with 34 rows (all H&E and IF samples) and columns for sample ID,
#     GUID, donor, sex, age, interveiw_date, image file, alignment file, fastq
#     files, spaceranger JSON, and loupe version
#   - One with 10 rows (one per donor) with sample ID, GUID, donor, sex, age,
#     and interview_date
#
#   This info will be borrowed by other R scripts in this directory

library(here)
library(tidyverse)
library(sessioninfo)

he_sample_info_path = here(
    "code", "13_submission", "spe_sample_info.csv"
)
he_sr_params_path = here(
    'code', '13_submission', 'spaceranger_parameters.txt'
)

out_dir = here('processed-data', '12_submission', 'nda-sub')

#guids = c(
#    "NDAR_INVRP206TL0", "NDAR_INVBD397ZME", "NDAR_INVZY012AKR", "NDAR_INVAA376RGT",
#    "NDAR_INVFD524YFZ", "NDAR_INVZP593AY8", "NDAR_INVAZ557YTF", "NDAR_INVVK615CK3",
#    "NDAR_INVEM171HHD", "NDAR_INVCG884PYU", "NDAR_INVHJ624JXT", "NDAR_INVDC739FJR",
#    "NDAR_INVFW122YUJ", "NDAR_INVAY229LWZ", "NDAR_INVTP473NFX", "NDAR_INVKC614AGH",
#    "NDAR_INVZY732JB8", "NDAR_INVPZ981ZZ2", "NDAR_INVFE904CRQ", "NDAR_INVRM926ZDX",
#    "NDAR_INVJB957KY1", "NDAR_INVZN731HBR", "NDAR_INVZN882LVK", "NDAR_INVTP875BX7",
#    "NDAR_INVCR898RRX", "NDAR_INVYA339AET", "NDAR_INVFV118XP1", "NDAR_INVYR277NZ1",
#    "NDAR_INVEM185UZJ", "NDAR_INVPF418GU7", "NDAR_INVJC255CPH", "NDAR_INVWW628YUJ",
#    "NDAR_INVFN378RDU", "NDAR_INVAR575HG5", "NDAR_INVWA002GHK", "NDAR_INVWV344CZ1",
#    "NDAR_INVWU602NXV", "NDAR_INVNU892JEJ", "NDAR_INVDW106VHW", "NDAR_INVHJ529LBR",
#    "NDAR_INVXV379VMU", "NDAR_INVNC002PAJ", "NDAR_INVDL298AEP", "NDAR_INVUZ538VBH",
#    "NDAR_INVVP296YRN", "NDAR_INVVL351RJA", "NDAR_INVXC538HR6", "NDAR_INVFW569JLK",
#    "NDAR_INVBF497DR8", "NDAR_INVLX496NWU", "NDAR_INVNE717BNC", "NDAR_INVFU506BA7",
#    "NDAR_INVKV771WVX", "NDAR_INVML104XA7", "NDAR_INVEM397THU", "NDAR_INVKF088EDN",
#    "NDAR_INVGK846MAM", "NDAR_INVUP036XFE", "NDAR_INVAN308YBL", "NDAR_INVXX654YB3",
#    "NDAR_INVEW232NT4", "NDAR_INVUD985ZPF", "NDAR_INVVG855PYQ"
#)

guids <- unlist(read.csv("guid-list_in-order.csv", header = FALSE), use.names = FALSE)

################################################################################
#   H&E samples
################################################################################

he_info1 = read_csv(he_sample_info_path, show_col_types = FALSE) |>
  rename(src_subject_id = subjects) |>
  mutate(
    interview_age = as.integer(round(age * 12, 0)),
    sample_id = str_replace(sample_id, '_2$', '')
  ) |>
  select(sample_id, src_subject_id, sex, interview_age)

he_info = read_tsv(
  he_sr_params_path, col_names = FALSE, show_col_types = FALSE
) |>
  mutate(
    sample_id_old = str_extract(X1, '^Br[0-9]{4}'),
    sample_id = sprintf('%s_%s', X2, X3),
    loupe_version = "8.0.0",
    fastq_files = sapply(
      X6,
      function(x) paste(list.files(x, full.names = TRUE), collapse = ';')
    ),
    stain = "HE",
    interview_date = '09/24/2024',
    spaceranger_json = here(
      "processed-data", "01_spaceranger", sample_id, "outs", "spatial",
      "scalefactors_json.json"
    )
  ) |>
    rename(image_file = X4, alignment_file = X5) |>
    select(
      sample_id_old, sample_id, loupe_version, fastq_files, image_file,
      alignment_file, interview_date, stain, spaceranger_json
    ) #|>
    # #   Add phenotype data
    # left_join(he_info, by = 'sample_id') |>
    # select(-sample_id_old) |>
    # #   Add GUIDs
    # left_join(
    #   tibble(
    #     src_subject_id = unique(he_info$src_subject_id),
    #     subjectkey = guids
    #   ),
    #   by = 'src_subject_id'
    # )
#he_info <- he_info1 %>% left_join(he_info, by = 'sample_id_old')

he_info <- left_join(he_info, he_info1, by ='sample_id') |>
  select(-sample_id_old) |>
  #   Add GUIDs
  left_join(
    tibble(
      src_subject_id = unique(he_info1$src_subject_id),
      subjectkey = guids
    ),
    by = 'src_subject_id'
  )



################
################

# he_info = read_csv(he_sample_info_path, show_col_types = FALSE) |>
#     rename(src_subject_id = subjects) |>
#     mutate(
#         interview_age = as.integer(round(age * 12, 0)),
#         sample_id_old = str_replace(sample_id, '_2$', '')
#     ) |>
#     select(sample_id_old, src_subject_id, sex, interview_age)
#
# all_json_paths = file.path(
#     list.files(
#         here("processed-data", "01_spaceranger", "spaceranger-all"), full.names = TRUE,
#         pattern = '^V'
#     ),
#     "outs", "spatial", "scalefactors_json.json"
# )
# stopifnot(all(file.exists(all_json_paths)))
#
# he_info = read_delim(
#         he_sr_params_path, col_names = FALSE, delim = " ",
#         show_col_types = FALSE
#     ) |>
#     mutate(
#         sample_id_old = str_extract(X1, 'Br[0-9]{4}'),
#         sample_id = sprintf('%s_%s', X2, X3),
#         loupe_version = "6.2.0",
#         fastq_files = str_replace_all(X6, ',', ';'),
#         #   Date of Visium sequencing
#         interview_date = '03/14/2022',
#         stain = "H&E",
#         spaceranger_json = sapply(
#             sample_id, function(x) all_json_paths[grep(x, all_json_paths)]
#         )
#     ) |>
#     rename(image_file = X4, alignment_file = X5) |>
#     select(
#         sample_id_old, sample_id, loupe_version, fastq_files, image_file,
#         alignment_file, interview_date, stain, spaceranger_json
#     ) |>
#     #   Add phenotype data
#     left_join(he_info, by = 'sample_id_old') |>
#     select(-sample_id_old) |>
#     #   Add GUIDs
#     left_join(
#         tibble(
#             src_subject_id = unique(he_info$src_subject_id),
#             subjectkey = guids
#         ),
#         by = 'src_subject_id'
#     )

################################################################################
#   Phenotype data that applies to many NDA data structures
################################################################################

#   Just one row per donor
pd = he_info |>
    group_by(src_subject_id) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(src_subject_id, subjectkey, interview_age, interview_date, sex)

write_csv(pd, file.path(out_dir, "pheno_data.csv"))

################################################################################
#   IF samples
################################################################################


################################################################################
#   Combine and export
################################################################################

#   At this point the same info should be collected for H&E and IF, and all
#   values should be non-NULL
stopifnot(!any(is.na(he_info)))

sample_info = he_info

#   All file paths should exist
stopifnot(all(file.exists(sample_info$spaceranger_json)))
stopifnot(all(file.exists(sample_info$image_file)))
stopifnot(all(file.exists(sample_info$alignment_file)))

write_csv(sample_info, file.path(out_dir, "imaging_sample_info.csv"))

session_info()
