# Decoding and Encoding Phylogenetic Classification Annotations

Functions to represent, decode and encode phylogenetic classification
annotations used in FASTA files by RDP and the Greengenes project.

## Usage

``` r
decode_Greengenes(annotation)

GenClass16S(
  Kingdom = NA,
  Phylum = NA,
  Class = NA,
  Order = NA,
  Family = NA,
  Genus = NA,
  Species = NA,
  Otu = NA,
  Org_name = NA,
  Id = NA
)

encode_Greengenes(classification)

decode_RDP(annotation)

encode_RDP(classification)
```

## Arguments

- annotation:

  Annotation from a FASTA file containing the classification
  information.

- Kingdom:

  Name of the kingdom to which the organism belongs.

- Phylum:

  Name of the phylum to which the organism belongs.

- Class:

  Name of the class to which the organism belongs.

- Order:

  Name of the order to which the organism belongs.

- Family:

  Name of the family to which the organism belongs.

- Genus:

  Name of the genus to which the organism belongs.

- Species:

  Name of the species to which the organism belongs.

- Otu:

  Name of the otu to which the organism belongs.

- Org_name:

  Name of the organism.

- Id:

  ID of the sequence.

- classification:

  A `data.frame` created with `GenClass16S()` with the classification
  information.

## Value

`GenClass16S()` and `decodeX()` return a `data.frame`. `encodeX()`
returns a string with the corresponding annotation.

## Examples

``` r

seq <- readRNAStringSet(system.file("examples/RNA_example.fasta",
    package = "rRDP"
))

### the FASTA annotation is read as names. This data has a Greengenes format
### annotation
names(seq)
#> [1] "1675 AB015560.1 deep-sea sediment clone BD4-10 k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfobacterales; f__Nitrospinaceae; g__Nitrospina; otu_3187"                                 
#> [2] "4399 D14432.1 Rhodovibrio salinarum str. NCIMB2243 k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Rhodospirillaceae; g__Rhodovibrio; s__Rhodovibrio salinarum; otu_2816"
#> [3] "4403 X72908.1 Roseococcus thiosulfatophilus str. RB-3 Yurkov strain Drews k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Acetobacteraceae; g__Roseococcus; otu_2785"    
#> [4] "4404 AF173825.1 Antarctic clone LB3-94 k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Acetobacteraceae; g__Roseococcus; otu_2785"                                       
#> [5] "4411 Y07647.2 Drentse grassland soil clone vii k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Acetobacteraceae; Unclassified; otu_2752"                                 

classification <- decode_Greengenes(names(seq))
classification
#>    Kingdom         Phylum               Class             Order
#> 1 Bacteria Proteobacteria Deltaproteobacteria Desulfobacterales
#> 2 Bacteria Proteobacteria Alphaproteobacteria  Rhodospirillales
#> 3 Bacteria Proteobacteria Alphaproteobacteria  Rhodospirillales
#> 4 Bacteria Proteobacteria Alphaproteobacteria  Rhodospirillales
#> 5 Bacteria Proteobacteria Alphaproteobacteria  Rhodospirillales
#>                           Family       Genus               Species  Otu
#> 1                 Nitrospinaceae  Nitrospina               unknown 3187
#> 2              Rhodospirillaceae Rhodovibrio Rhodovibrio salinarum 2816
#> 3               Acetobacteraceae Roseococcus               unknown 2785
#> 4               Acetobacteraceae Roseococcus               unknown 2785
#> 5 Acetobacteraceae; Unclassified     unknown               unknown 2752
#>                                                               Org_name   Id
#> 1                            AB015560.1_deep-sea_sediment_clone_BD4-10 1675
#> 2                        D14432.1_Rhodovibrio_salinarum_str._NCIMB2243 4399
#> 3 X72908.1_Roseococcus_thiosulfatophilus_str._RB-3_Yurkov_strain_Drews 4403
#> 4                                    AF173825.1_Antarctic_clone_LB3-94 4404
#> 5                            Y07647.2_Drentse_grassland_soil_clone_vii 4411

### look at the Genus of all sequences
classification[, "Genus"]
#> [1] "Nitrospina"  "Rhodovibrio" "Roseococcus" "Roseococcus" "unknown"    

### to train the RDP classifier, the annotations need to be in RDP format
annotation <- encode_RDP(classification)
names(seq) <- annotation
seq
#> RNAStringSet object of length 5:
#>     width seq                                               names               
#> [1]  1481 AGAGUUUGAUCCUGGCUCAGAAC...GGUGAAGUCGUAACAAGGUAACC 1675 Root;Bacteri...
#> [2]  1404 GCUGGCGGCAGGCCUAACACAUG...CACGGUAAGGUCAGCGACUGGGG 4399 Root;Bacteri...
#> [3]  1426 GGAAUGCUNAACACAUGCAAGUC...AACAAGGUAGCCGUAGGGGAACC 4403 Root;Bacteri...
#> [4]  1362 GCUGGCGGAAUGCUUAACACAUG...UACCUUAGGUGUCUAGGCUAACC 4404 Root;Bacteri...
#> [5]  1458 AGAGUUUGAUUAUGGCUCAGAGC...UGAAGUCGUAACAAGGUAACCGU 4411 Root;Bacteri...

### now we can train the classifier
customRDP <- trainRDP(seq, dir = "sample_classifier")
customRDP
#> RDPClassifier
#> Location: /home/runner/work/rRDP/rRDP/docs/reference/sample_classifier 

## clean up
removeRDP(customRDP)
```
