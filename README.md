# Super Maximal Repeats Kernel

The algorithm we present is a strategy for finding a good reference for Relative Lempel Ziv by leveraging Super Maximal Repeats. Let $$S[1,n] in \Sigma^n$$ be a string. 

A supermaximal repeat (SMR) is a repeat that is not contained inside any other repeat. The notions of repeats and supermaximal repeats can naturally be extended to string sets. One way of doing this is to convert the set $${S1, ... , Sk}$$ to the string $$S' = S1 #1 S2 #2 ... #(k-1) Sk$$, where $$#1, ... , #k$$ are all distinct and not part of the original alphabet.
However, to construct a reference, we can define the following function:

We define k(S) for a string S as the following function: 

1. We define the set SMR' as all the Supermaximal Repeats of S
2. We join (concatenate) each element of SMR', while making sure that no suffix of SMR'[k-1] is a prefix of SMR'[k] for any k in the interval 0 <= k <= len(SMR'). For periodic strings, we write only the primitive root of a substring i.e., pr(xyxyxy) = xy

As such, we can make $$k(k(S), k(k(...(k(S))..)) => k^i(S) = eps$$, which corresponds to the i-eth recursion of the kernelization process.

## Bit-sizes

We use the LZ77/RLZ variant where each phrase is the shortest substring of length l that does not appear before. Each phrase can therefore be encoded as a triple (i,l-1,c) where i is the pointer pointing to the source of the phrase, l-1 is the number of copied characters and c is the extra character. In practice, i and c are encoded with no compression, while l is encoded with Elias-Fano.

To simplify notation, when writing log(x) we actually mean max(1, ceil[log2(x)])

The bit-size of RLZ is:

$$
Rlog(Sigma) + zlog(R) + zlog(n/z) + zlog(Sigma) + 2z 
$$

The algorithm runs k(S) until BITSIZE[RLZ(k^i+1(S))] > BITSIZE[RLZ(k^i(S))]. For now, we do not build the Relative Lempel Ziv data struture to memory, but we plan to.  

## Running the code

Clone the repository with git clone --recursive. To build the code, run make. Then, save the file you want to encode in data/ and run ./build/smr < data/file.txt.
