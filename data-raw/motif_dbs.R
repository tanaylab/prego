library(tidyverse)

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

JASPAR_motifs <- get_dataset_pssms("jaspar") %>%
    select(motif = track, pos, A, C, G, T) %>%
    as_tibble()
usethis::use_data(JASPAR_motifs, overwrite = TRUE)

HOMER_motifs <- get_dataset_pssms("homer") %>%
    select(motif = track, pos, A, C, G, T) %>%
    as_tibble()
usethis::use_data(HOMER_motifs, overwrite = TRUE)

JOLMA_motifs <- get_dataset_pssms("jolma") %>%
    select(motif = track, pos, A, C, G, T) %>%
    as_tibble()
usethis::use_data(JOLMA_motifs, overwrite = TRUE)
