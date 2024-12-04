library(tidyverse)

proc_jaspar_format <- function(fn, sep_ids = TRUE) {
    js_file <- tempfile()
    download.file(fn, js_file)
    jaspar <- PWMEnrich::readMotifs(js_file, remove.acc = TRUE)

    if (sep_ids) {
        # add the id prefix on redundant motif names
        names(jaspar) <- tibble(motif = names(jaspar)) %>%
            separate(motif, c("id", "motif"), sep = "\\t") %>%
            add_count(motif) %>%
            mutate(new_motif = ifelse(n == 1, motif, paste0(motif, ".", id))) %>%
            pull(new_motif)
    }

    motifs <- imap_dfr(jaspar, ~ {
        as.data.frame(t(.x) / rowSums(t(.x))) %>%
            mutate(pos = 1:n() - 1, motif = .y) %>%
            select(motif, pos, A, C, G, T)
    }) %>% as_tibble()

    return(motifs)
}

# JASPAR
JASPAR_motifs <- proc_jaspar_format("https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt", sep_ids = TRUE)
usethis::use_data(JASPAR_motifs, overwrite = TRUE, compress = "xz")

# HOCOMOCO
HOCOMOCO_motifs_human <- proc_jaspar_format("https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_jaspar_format.txt", sep_ids = FALSE)

HOCOMOCO_motifs_mouse <- proc_jaspar_format("https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/MOUSE/mono/HOCOMOCOv11_core_MOUSE_mono_jaspar_format.txt", sep_ids = FALSE)

HOCOMOCO_motifs <- bind_rows(HOCOMOCO_motifs_human, HOCOMOCO_motifs_mouse)
usethis::use_data(HOCOMOCO_motifs, overwrite = TRUE, compress = "xz")

# HOMER and Jolma
get_dataset_pssms <- function(datasets) {
    all_pssms <- mcATAC::get_available_pssms(datasets_of_interest = datasets)
    purrr::map_dfr(datasets, function(d) {
        all_pssms$datasets[[d]] %>%
            left_join(all_pssms$keys[[d]], by = "key") %>%
            mutate(
                dataset = d
            )
    })
}

misha.ext::gset_genome("mm10")

# JASPAR_motifs <- get_dataset_pssms("jaspar") %>%
#     select(motif = track, pos, A, C, G, T) %>%
#     as_tibble()
# usethis::use_data(JASPAR_motifs, overwrite = TRUE)

HOMER_motifs <- get_dataset_pssms("homer") %>%
    select(motif = track, pos, A, C, G, T) %>%
    as_tibble()
usethis::use_data(HOMER_motifs, overwrite = TRUE, compress = "xz")

JOLMA_motifs <- get_dataset_pssms("jolma") %>%
    select(motif = track, pos, A, C, G, T) %>%
    as_tibble()
usethis::use_data(JOLMA_motifs, overwrite = TRUE, compress = "xz")

MOTIF_DB <- create_motif_db(all_motif_datasets())
usethis::use_data(MOTIF_DB, overwrite = TRUE, compress = "xz")
