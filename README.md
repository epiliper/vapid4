# Changes

## 2024-Feb-13: removed functionality to pull metadata from .csv

Removed do_meta_data() function, as it was only ever used to pull strain names already mentioned in fasta. Instead, the header from the fasta is used as the full name for generating .fsa files. 
