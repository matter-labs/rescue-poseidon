# Rescue and Poseidon
## Overview
This repo contains implementations of arithmetization oriented hash functions that constructed by using a sponge construction over prime field. These hash functions are optimal in the constraint systems while also supporting different finite fields which supported by bellman. 

## Sponge 
The sponge construction proceeds in two phases: the absorbing phase followed by the squeezing phase.
- In the absorbing phase, the `r` elements input message blocks are summed into the first `r` elements of the state, interleaved with applications of the permutation `f`. When all message blocks are processed, the sponge construction switches to the squeezing phase.

- In the squeezing phase, the first `r` elements of the state are returned as output blocks, interleaved with applications of the permutation `f` . The number of output blocks is chosen at will by the user.

The last `c` elements of the state are never directly affected by the input blocks and are never output during the squeezing phase.

## Padding Strategies
Padding prevents trivial collisions. Each hash function nearly uses same padding strategies. The only difference is that Rescue Prime requires no padding for fixed length input. Rescue and Poseidon require same padding rule for variable length input.

### Fixed Length
The capacity value is `length x (^264 ) + (o − 1)` where `o` the output length. The padding consists of the field elements being 0.

### Variable Length
Padding is necessary for variable-length inputs, even if the input is already a multiple of the rate in length. The capacity value is `2^64 + (o − 1)` where `o` the output length.The padding consists of one field element being 1, and the remaining elements being 0

## Hash Functions
### Poseidon
`let hasher = PoseidonHasher::default();`

### Rescue
`let hasher = RescueHasher::default();`

### Rescue Prime
`let hasher = RescuePrimeHasher::default();`

### Example
#### Stateful Sponge
```
    let preimage = vec![..];
    let mut hasher = <Hasher>::default();
    hasher.absorb_multi(preimage);
    let hash = hasher.squeeze();
```
#### Fixed Length Hashing
```
    let hash = poseidon_hash(&preimage);
```


## Benchmarks
`CPU: 3,1 GHz Intel Core i5`

| hashes    | 1x permutation runtime (μs) | 1x permutation gates | number of rounds |
| --- | -------- | -------- | -------- |
| Poseidon   | 13     | 166     | 8f + 33p     |
| Rescue   | 680     | 266     | 44f     |
| Rescue Prime   | 300     | 104     | 9f     |



## TODO 

## References
- [1] [Cryprographic sponge functions](https://keccak.team/files/CSF-0.1.pdf)
- [2] [The sponge and duplex constructions](https://keccak.team/sponge_duplex.html)
- [3] [STARK Friendly Hash – Survey and Recommendation](https://eprint.iacr.org/2020/948.pdf)
- [4] [MARVELlous: a STARK-Friendly Family of Cryptographic Primitives](https://eprint.iacr.org/2018/1098.pdf)
- [5] [POSEIDON: A New Hash Function for Zero-Knowledge Proof Systems](https://eprint.iacr.org/2019/458.pdf)
- [6] [Rescue-Prime: a Standard Specification (SoK)](https://eprint.iacr.org/2020/1143.pdf)