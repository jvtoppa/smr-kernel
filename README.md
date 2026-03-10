# Super Maximal Repeats Kernel

The algorithm we present is a strategy for finding a good reference for Relative Lempel Ziv by leveraging Super Maximal Repeats. Let $S[1,n] \in \Sigma^n$ be a string. 

A supermaximal repeat (SMR) is a repeat that is not contained inside any other repeat. The notions of repeats and supermaximal repeats can naturally be extended to string sets. One way of doing this is to convert the set ${S_1, ... , S_k}$ to the string $` S' = S_1 \$_1 S_2 \$_2 \dots \$_{k-1} S_k`$, where $`\$_1, ... , \$_k`$ are all distinct and not part of the original alphabet.
However, to construct a reference, we can define the following function:

We define k(S) for a string S as the following function: 

1. We define the set SMR' as all the Supermaximal Repeats of S
2. We join (concatenate) each element of SMR', while making sure that no suffix of $SMR'[k-1]$ is a prefix of $SMR'[k]$ for any $k$ in the interval $1 <= k <= len(SMR')$. For periodic strings, we write only the primitive root of a substring i.e., $pr(xyxyxy) = xy$

In other words, 

$$
\prod_{i=1}^{n}pr(SMR'_i) = k(S)
$$

As such, it is easy to notice that $k(k(S), k(k(...(k(S))..)) => k^m(S) = \epsilon$, which corresponds to the m-th (last) recursion step of the kernelization process.

## Bit-sizes

We use the LZ77/RLZ ($\bar{Z}$) variant where each phrase is the shortest substring of length $\ell$ that does not appear before. Each phrase can therefore be encoded as a triple $(i,\ell-1,c)$ where $i$ is the pointer pointing to the source of the phrase, $\ell-1$ is the number of copied characters and c is the extra character. In practice, i and c are encoded with no compression, while $\ell$ is encoded with Elias-Fano.

To simplify notation, when writing $\log(x)$ we actually mean $max(1, ceil[\log_2(x)])$

The bit-size of RLZ is:

$$
R\log(\Sigma) + z\log(R) + z\log(n/z) + z\log(\Sigma) + 2z = B_\bar{Z}
$$

We have noticed experimentally that there exists a global minimum for the size of RLZ with this strategy. As such, we can simply do a recursion call until $$B_\bar{Z}[\bar{Z}(k^{i+1}(S))] > B_\bar{Z}[RLZ[\bar{Z}(k^i(S))]$$. For now, we do not build the Relative Lempel Ziv data struture to memory, but we plan to.  

## Running the code

Clone the repository with git clone --recursive. To build the code, run make. Then, save the file you want to encode in data/ and run ./build/smr < data/file.txt.
