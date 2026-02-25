# Generate kmers

This function generates all possible kmers considering the gap length.
Gaps are represented by 'N'.

## Usage

``` r
generate_kmers(kmer_length, max_gap = 0, min_gap = 0)
```

## Arguments

- kmer_length:

  The number of non-gap characters in motifs that will be screened.

- max_gap:

  The maximum length of a gap to be considered in the pattern. Default:
  0

- min_gap:

  The minimum length of a gap to be considered in the pattern. Default:
  0

## Value

A vector of all possible kmers considering the gap length.

## Examples

``` r
# Generate kmers of length 2 without any gaps
generate_kmers(2)
#>  [1] "TT" "CT" "GT" "AT" "TC" "CC" "GC" "AC" "TG" "CG" "GG" "AG" "TA" "CA" "GA"
#> [16] "AA"

# Generate kmers of length 3 with a single gap (1 'N') at any position
generate_kmers(3, min_gap = 1, max_gap = 1)
#>   [1] "NTT" "NTT" "NTT" "NTT" "NCT" "NCT" "NCT" "NCT" "NGT" "NGT" "NGT" "NGT"
#>  [13] "NAT" "NAT" "NAT" "NAT" "NTC" "NTC" "NTC" "NTC" "NCC" "NCC" "NCC" "NCC"
#>  [25] "NGC" "NGC" "NGC" "NGC" "NAC" "NAC" "NAC" "NAC" "NTG" "NTG" "NTG" "NTG"
#>  [37] "NCG" "NCG" "NCG" "NCG" "NGG" "NGG" "NGG" "NGG" "NAG" "NAG" "NAG" "NAG"
#>  [49] "NTA" "NTA" "NTA" "NTA" "NCA" "NCA" "NCA" "NCA" "NGA" "NGA" "NGA" "NGA"
#>  [61] "NAA" "NAA" "NAA" "NAA" "TNT" "CNT" "GNT" "ANT" "TNT" "CNT" "GNT" "ANT"
#>  [73] "TNT" "CNT" "GNT" "ANT" "TNT" "CNT" "GNT" "ANT" "TNC" "CNC" "GNC" "ANC"
#>  [85] "TNC" "CNC" "GNC" "ANC" "TNC" "CNC" "GNC" "ANC" "TNC" "CNC" "GNC" "ANC"
#>  [97] "TNG" "CNG" "GNG" "ANG" "TNG" "CNG" "GNG" "ANG" "TNG" "CNG" "GNG" "ANG"
#> [109] "TNG" "CNG" "GNG" "ANG" "TNA" "CNA" "GNA" "ANA" "TNA" "CNA" "GNA" "ANA"
#> [121] "TNA" "CNA" "GNA" "ANA" "TNA" "CNA" "GNA" "ANA" "TTN" "CTN" "GTN" "ATN"
#> [133] "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN"
#> [145] "TTN" "CTN" "GTN" "ATN" "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN"
#> [157] "TAN" "CAN" "GAN" "AAN" "TTN" "CTN" "GTN" "ATN" "TCN" "CCN" "GCN" "ACN"
#> [169] "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN" "TTN" "CTN" "GTN" "ATN"
#> [181] "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN"

# Generate kmers of length 3 with a gap of 1 to 2 'N's at any position
generate_kmers(3, min_gap = 1, max_gap = 2)
#>   [1] "NTT" "NTT" "NTT" "NTT" "NCT" "NCT" "NCT" "NCT" "NGT" "NGT" "NGT" "NGT"
#>  [13] "NAT" "NAT" "NAT" "NAT" "NTC" "NTC" "NTC" "NTC" "NCC" "NCC" "NCC" "NCC"
#>  [25] "NGC" "NGC" "NGC" "NGC" "NAC" "NAC" "NAC" "NAC" "NTG" "NTG" "NTG" "NTG"
#>  [37] "NCG" "NCG" "NCG" "NCG" "NGG" "NGG" "NGG" "NGG" "NAG" "NAG" "NAG" "NAG"
#>  [49] "NTA" "NTA" "NTA" "NTA" "NCA" "NCA" "NCA" "NCA" "NGA" "NGA" "NGA" "NGA"
#>  [61] "NAA" "NAA" "NAA" "NAA" "TNT" "CNT" "GNT" "ANT" "TNT" "CNT" "GNT" "ANT"
#>  [73] "TNT" "CNT" "GNT" "ANT" "TNT" "CNT" "GNT" "ANT" "TNC" "CNC" "GNC" "ANC"
#>  [85] "TNC" "CNC" "GNC" "ANC" "TNC" "CNC" "GNC" "ANC" "TNC" "CNC" "GNC" "ANC"
#>  [97] "TNG" "CNG" "GNG" "ANG" "TNG" "CNG" "GNG" "ANG" "TNG" "CNG" "GNG" "ANG"
#> [109] "TNG" "CNG" "GNG" "ANG" "TNA" "CNA" "GNA" "ANA" "TNA" "CNA" "GNA" "ANA"
#> [121] "TNA" "CNA" "GNA" "ANA" "TNA" "CNA" "GNA" "ANA" "TTN" "CTN" "GTN" "ATN"
#> [133] "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN"
#> [145] "TTN" "CTN" "GTN" "ATN" "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN"
#> [157] "TAN" "CAN" "GAN" "AAN" "TTN" "CTN" "GTN" "ATN" "TCN" "CCN" "GCN" "ACN"
#> [169] "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN" "TTN" "CTN" "GTN" "ATN"
#> [181] "TCN" "CCN" "GCN" "ACN" "TGN" "CGN" "GGN" "AGN" "TAN" "CAN" "GAN" "AAN"
#> [193] "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT"
#> [205] "NNT" "NNT" "NNT" "NNT" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC"
#> [217] "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNG" "NNG" "NNG" "NNG"
#> [229] "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG"
#> [241] "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA"
#> [253] "NNA" "NNA" "NNA" "NNA" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [265] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [277] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [289] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [301] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [313] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"

# Generate kmers of length 3 with a gap of 2 'N's at any position
generate_kmers(3, min_gap = 2, max_gap = 2)
#>   [1] "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT" "NNT"
#>  [13] "NNT" "NNT" "NNT" "NNT" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC"
#>  [25] "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNC" "NNG" "NNG" "NNG" "NNG"
#>  [37] "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG" "NNG"
#>  [49] "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA" "NNA"
#>  [61] "NNA" "NNA" "NNA" "NNA" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#>  [73] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#>  [85] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#>  [97] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [109] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
#> [121] "TNN" "CNN" "GNN" "ANN" "TNN" "CNN" "GNN" "ANN"
```
