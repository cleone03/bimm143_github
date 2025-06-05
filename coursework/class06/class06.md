# Class 6: R Functions
Christopher Leone (A16731724)

- [Introduction to Functions](#introduction-to-functions)
- [Writing a function to return a DNA sequence of a user specified
  length.](#writing-a-function-to-return-a-dna-sequence-of-a-user-specified-length)
- [Writing a function to return protein sequences of a user specified
  length.](#writing-a-function-to-return-protein-sequences-of-a-user-specified-length)

## Introduction to Functions

Let’s start writing our first function to add some numbers.

> Every R function has 3 important components:

- 1)  the name, we get to choose.
- 2)  input arguments (there could be any amount)
- 3)  function body (where the code is).

``` r
add <- function(x, y){
  x + y
}
```

I can now just use this function, as shown here:

``` r
x <- c(1:10)
y <- c(11:20)
add(x, y)
```

     [1] 12 14 16 18 20 22 24 26 28 30

Functions can have “required” input arguments and “optional” input
arguments. Optional arguments are defined with an “equals default value”
(`y=10`) in the function declination.

## Writing a function to return a DNA sequence of a user specified length.

1)  We can use the `sample()` function for ideas.

``` r
# genDNA <- function(size=5){}
students <- c("Jeff", "Jeremy", "Peter")
sample(students, size=5, replace = TRUE)
```

    [1] "Jeremy" "Jeff"   "Jeff"   "Peter"  "Peter" 

2)  Let’s try using this function to create nucleotide sequences.

``` r
# Must have "replace = TRUE" for repeats!
bases <- c("A", "T", "C", "G")
sample(bases, size=10, replace = TRUE)
```

     [1] "C" "C" "A" "A" "T" "A" "T" "T" "C" "C"

3)  Now we have a working snippet, so let’s use it to generate our
    `function()`.

``` r
genDNA <- function(size=10){
  bases <- c("A", "T", "C", "G")
sample(bases, size, replace = TRUE)
}
genDNA(50)
```

     [1] "C" "G" "G" "T" "A" "A" "C" "C" "G" "T" "C" "C" "G" "A" "T" "A" "T" "C" "T"
    [20] "A" "G" "T" "T" "A" "C" "G" "T" "T" "C" "T" "G" "A" "C" "C" "A" "G" "T" "G"
    [39] "T" "T" "A" "C" "G" "G" "T" "C" "G" "A" "A" "G"

4)  Finally, to polish the output so it returns a sequence such as
    “ATGCATA…” if we want, otherwise return the raw output. Here, I want
    the collapsed sequence:

``` r
genDNA <- function(size=10, collapsed = FALSE){
  bases <- c("A", "T", "C", "G")
  seq <- sample(bases, size, replace = TRUE)
  if(!collapsed) {
    paste(seq)
  } else {
    paste(seq, collapse="")
  }
}
genDNA(50, TRUE)
```

    [1] "ATCTTGGGATAAAATACCCATCGGAAATGCCGAGAATAGAAGATCTGCCC"

## Writing a function to return protein sequences of a user specified length.

Similar to `genDNA()`, let’s start using the sample function and build
up from there. Main difference includes bringing the amino acids into
play rather than 4 nucleotides.

We can get the set of 20 natural amino acids from the **bio3d** package.

``` r
head(bio3d::aa.table, 5)
```

        aa3 aa1    mass      formula          name
    ALA ALA   A  71.078   C3 H5 N O1       Alanine
    ARG ARG   R 157.194 C6 H13 N4 O1      Arginine
    ASN ASN   N 114.103  C4 H6 N2 O2    Asparagine
    ASP ASP   D 114.079   C4 H4 N O3 Aspartic Acid
    CYS CYS   C 103.143 C3 H5 N O1 S       Cystein

Once called, let’s use `aa.table` as our database for amino acids, and
complete the function `genProt()`.

``` r
genProt <- function(size=6, collapsed=TRUE){
  amino <- bio3d::aa.table$aa1[1:20]
  seq2 <- sample(amino, size, replace = TRUE)
  if(collapsed) {
    paste(seq2, collapse="")
  } else {
    paste(seq2)
  }
}
genProt(size=25, TRUE)
```

    [1] "HNQHHQPMGQLFREDYNEIVADIVH"

Now, what if I wanted to generate a sequence of random length 6-12?

``` r
genProt2 <- function(collapsed=TRUE){
  # to set a random size of protein each time of length 6-12.
  randomSize <- sample(c(6:12), 1)
  
  # our protein parameters
  amino <- bio3d::aa.table$aa1[1:20]
  seq3 <- sample(amino, size=randomSize, replace = TRUE)
  
  # do we want a collapsed AA sequence?
  if(collapsed) {
    paste0(seq3, collapse="")
  } else {
    paste(seq3)
  }
}
genProt2(T)
```

    [1] "IYQEPRGVN"

And if I wanted to print protein sequences of all lengths between 6-12?

``` r
prot <- sapply(6:12,FUN=genProt)
prot
```

    [1] "VAGFWS"       "PIFEITG"      "FTKCMEHA"     "HFTIMFEEW"    "SPEVIIHQRI"  
    [6] "NRGTMTNMTFE"  "HRDTGNWDSNRT"

It would also be cool and useful if I could get these in FASTA format
for easy searching. Let’s combine the functions of both the `cat()` and
`paste()` functions.

``` r
fasta <- paste(">ID.", 6:12, "\n", prot, sep="")
cat(fasta, sep="\n")
```

    >ID.6
    VAGFWS
    >ID.7
    PIFEITG
    >ID.8
    FTKCMEHA
    >ID.9
    HFTIMFEEW
    >ID.10
    SPEVIIHQRI
    >ID.11
    NRGTMTNMTFE
    >ID.12
    HRDTGNWDSNRT

Through a BLASTp query, I found that proteins 6 & 7 are NOT unique and
can be found matched in the database (100% coverage and identity), but
proteins 8-12 are in fact unique!
