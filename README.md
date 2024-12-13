# Biocompatibility-yyc-codec
This is part of my undergraduate research project on a DNA storage model in Bacillus subtilis. This codec is based on the YING YANG codec developed by the BGI Research team ( https://github.com/ntpz870817/DNA-storage-YYC ). Here I will make some adaptations with the aim of improving the codec for in vivo archival.

### Last Commit Resume

13.12.24
-> Implemented a change in the Validation logic. Now the is thinking in veryfing Palindromes in dsDNA instead of thinking about the secondery structures formed by RNA/DNA fold.
-> More easty for me to understand comments in the code
-> NEED TO FIX DECODIFICATION AND IMPLEMENT A WORKING COMPARISION OF BITS FROM OG FILE AND DECODED FILE, ALSO THE CSV WITH THE NUCLEOTIDES AND CORRESPONDING DATA

## What to do Next
-> Enconding Log and Segments of DNA csv table

-> Index generations using more caracters than numbers

-> Validations of Incorporation // and how to deal with BAD/GOOD 0:1 data ratio

---> Nullomers check

---> Palindrome chek

---> Dual Strand DNA check

-> Error Correction Validation for Decoding process

-> Bits/base math

#### Extras

-> Simulations wiht https://dmci.xmu.edu.cn/dna/#/home ?? 

-> Organism B.subtillis mutations ratio ??

