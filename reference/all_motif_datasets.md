# Get a data frame of all the motif datasets bundled with prego

The data frame contain the PSSMs ('A', 'C', 'G' and 'T' columns), with
an additional column 'motif' containing the motif name. Individual
datasets are available within the package as `HOMER_motifs`,
`JASPAR_motifs`, `JOLMA_motifs`, and `HOCOMOCO_motifs`.

## Usage

``` r
all_motif_datasets()
```

## Value

a data frame which concatenates motifs from "HOMER", "JASPAR" and
"JOLMA". Motif names are prefixed with the dataset name, e.g.
"JASPAR.GATA4".

## References

- HOMER: :

  Heinz S, Benner C, Spann N, Bertolino E et al. Simple Combinations of
  Lineage-Determining Transcription Factors Prime cis-Regulatory
  Elements Required for Macrophage and B Cell Identities. Mol Cell 2010
  May 28;38(4):576-589. PMID: 20513432

- JASPAR: :

  Castro-Mondragon JA, Riudavets-Puig R, Rauluseviciute I, Berhanu Lemma
  R, Turchi L, Blanc-Mathieu R, Lucas J, Boddie P, Khan A, Manosalva
  Pérez N, Fornes O, Leung TY, Aguirre A, Hammal F, Schmelter D,
  Baranasic D, Ballester B, Sandelin A, Lenhard B, Vandepoele K,
  Wasserman WW, Parcy F, and Mathelier A JASPAR 2022: the 9th release of
  the open-access database of transcription factor binding profiles
  Nucleic Acids Res. 2022 Jan 7;50(D1):D165-D173.; doi:
  10.1093/nar/gkab1113

- JOLMA: :

  Jolma, A., Yin, Y., Nitta, K. et al. DNA-dependent formation of
  transcription factor pairs alters their binding specificity. Nature
  534, S15–S16 (2016).
  [doi:10.1038/nature18912](https://doi.org/10.1038/nature18912)

- HOCOMOCO: :

  Ivan V. Kulakovskiy; Ilya E. Vorontsov; Ivan S. Yevshin; Ruslan N.
  Sharipov; Alla D. Fedorova; Eugene I. Rumynskiy; Yulia A. Medvedeva;
  Arturo Magana-Mora; Vladimir B. Bajic; Dmitry A. Papatsenko; Fedor A.
  Kolpakov; Vsevolod J. Makeev: HOCOMOCO: towards a complete collection
  of transcription factor binding models for human and mouse via
  large-scale ChIP-Seq analysis. Nucl. Acids Res., Database issue,
  gkx1106 (11 November 2017).
  [doi:10.1093/nar/gkx1106](https://doi.org/10.1093/nar/gkx1106)

## Examples

``` r
all_motif_datasets()
#> # A tibble: 49,269 × 8
#>    motif        pos     A     C     G     T dataset motif_orig
#>    <chr>      <dbl> <dbl> <dbl> <dbl> <dbl> <chr>   <chr>     
#>  1 HOMER.AP_1     0 0.419 0.275 0.277 0.028 HOMER   AP_1      
#>  2 HOMER.AP_1     1 0.001 0.001 0.001 0.997 HOMER   AP_1      
#>  3 HOMER.AP_1     2 0.01  0.002 0.965 0.023 HOMER   AP_1      
#>  4 HOMER.AP_1     3 0.984 0.003 0.001 0.012 HOMER   AP_1      
#>  5 HOMER.AP_1     4 0.062 0.579 0.305 0.054 HOMER   AP_1      
#>  6 HOMER.AP_1     5 0.026 0.001 0.001 0.972 HOMER   AP_1      
#>  7 HOMER.AP_1     6 0.043 0.943 0.001 0.012 HOMER   AP_1      
#>  8 HOMER.AP_1     7 0.98  0.005 0.001 0.014 HOMER   AP_1      
#>  9 HOMER.AP_1     8 0.05  0.172 0.307 0.471 HOMER   AP_1      
#> 10 HOMER.AP_1     9 0.149 0.444 0.211 0.195 HOMER   AP_1      
#> # ℹ 49,259 more rows
```
