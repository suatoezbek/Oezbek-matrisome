# Molecular dynamics of the matrisome across sea anemone life history 
 
B. Gideon Bergheim, Alison G. Cole, Mandy Rettel, Frank Stein, Stefan Redl, Michael W. Hess, Aissam Ikmi and Suat Özbek* 

*Correspondence: suat.oezbek@cos.uni-heidelberg.de

---

## Supplementary Files

### input_proteomes
FASTA files of complete proteomes used as basis for *in silico* matrisome generation

### in_silico_matrisomes
#### from_literature
In silico matrisomes taken from the [Matrisome Project](https://sites.google.com/uic.edu/matrisome/home) webpage. IDs were looked up in the respectie sequence data bases. For Bos Taurus (where no IDs were provided) genes were searched by gene name.<br>
- *Homo sapiens*
- *Mus musculus*
- *Bos taurus*
- *Danio rerio*
- *Coturnix japonica*
- *Drosophila melanogaster*
- *C. elegans*
- *Schmidtea mediterranea*

#### in_this_publication
The following *in silico* matrisomes where generated with the methods described in the Materials and Methods and used as representative members of different taxa.

Chordata
- *Branchiostoma belcheri*
- *Branchiostoma floridae*
- *Branchiostoma lanceolatum*

Choanozoa
 - *Monosiga brevicollis*
 - *Salpingoeca rosetta*

 Cnidaria
 - *Acropora digitifera*
 - *Exaiptasia diaphana*
 - *Nematostella vectensis*
 - *Porites asteroides*
 - *Stylophora pistillata*
 - *Xenia spec.*
 - *Morbakka virulenta*
 - *Clytia hemispherica*
 - *Hydra vulgaris*
 - *Thelohanellus_kitauei*
 - *Aurelia aurita*
 - *Calvadosia cruxmelitensis*

Ctenophora
- *Mnemiopsis leidyi*
- *Pleurobrachia bachei*
- *Beroe ovata*

Placozoa
- *Trichoplax spec. (H2)*
- *Tricoplax adhaerens*

Porifera
- *Amphimedon queenslandica*
- *Ephydatia muelleri*

#### other_in_silico_matrisomes
In addition to the *in silico* matrisomes mentioned in this publication we also generated a number of additional *in silico* matrisomes using the same methods. These are not included in any of the analyses in the publication.

Annelida
- *Capitella teleta* (was excluded due to severe fragmentation)

Cnidaria (where excluded to avoid a Cnidarian bias as much as possible)
- *Acropora cervicornis*
- *Acropora hyacinthus*
- *Acropora millepora*
- *Acropora muricata*
- *Actinia tenebrosa*
- *Montipora aequituberculata*
- *Montipora foliosa*
- *Pocillopora damicornis*
- *Pocillopora verrucosa*
- *Renilla reniformis*
- *Hydra viridissima*
- *Hydractinia echinata*
- *Hydractinia symbiolongicarpus*
- *Myxobolus honghuensis*
- *Cassiopea xamachana*

### orthofinder results
This directory contains most of the output of the orthofinder run used to identify orthogroups. The extensive overview ove r the files can be found on the [OrthoFinder Tutorial](https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html).

### Source Code
This directory contains the python code whichas used to generate the in silico matrisomes and the following analyses. Some of these rely on additional files found in suppl_code_files.

### mass_spec_run
This directory contains the MS files and R analysis markdown files written by Frank Stein.