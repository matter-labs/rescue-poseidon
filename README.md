# Rescue and Poseidon
## Overview
This repo contains implementations of arithmetization oriented hash functions(Rescue, Poseidon, Rescue Prime) that constructed by a sponge construction over prime field for both out-of circuits and in-circuit usages. Each algebraic hash function uses same sponge construction with different round function or permutation function. Gadgets are optimal in the constraint systems while also supporting different scalar fields which supported by bellman. 

## Example
These are examples for Rescue but all of them are same for Poseidon and Rescue Prime as well.

### Fixed Length
```rust
    use rescue_poseidon::rescue_hash;

    let input = [Fr; 2]; // init fixed-length array

    let out = rescue_hash(&input);
```

### Variable Length
```rust
    use rescue_poseidon::rescue_hash_var_length;

    let input = [..]; // init some input values, input should be multiple of rate=2

    let out = rescue_hash_var_length(&input);    
```

### Gadget (Fixed Length)
```rust
    use rescue_poseidon::rescue_gadget;

    let input = [Num; 2]; // init fixed-length array

    let out = rescue_gadget(&input);    
```

### Gadget (Variable Length)
```rust
    use rescue_poseidon::rescue_gadget_var_length;

    let input = [..]; // init some input values, input should be multiple of rate=2

    let out = rescue_gadget_var_length(&input);    
```


## Benchmarks & Constraint System Costs
`CPU: 3,1 GHz Intel Core i5`

| hashes    | 1x permutation runtime (μs) | 1x permutation gates | number of rounds |
| --- | -------- | -------- | -------- |
| Poseidon   | 13     | 166     | 8f + 33p     |
| Rescue   | 680     | 266     | 44f     |
| Rescue Prime   | 300     | 104     | 9f     |



## References
- [1] [Cryprographic sponge functions](https://keccak.team/files/CSF-0.1.pdf)
- [2] [The sponge and duplex constructions](https://keccak.team/sponge_duplex.html)
- [3] [STARK Friendly Hash – Survey and Recommendation](https://eprint.iacr.org/2020/948.pdf)
- [4] [MARVELlous: a STARK-Friendly Family of Cryptographic Primitives](https://eprint.iacr.org/2018/1098.pdf)
- [5] [POSEIDON: A New Hash Function for Zero-Knowledge Proof Systems](https://eprint.iacr.org/2019/458.pdf)
- [6] [Rescue-Prime: a Standard Specification (SoK)](https://eprint.iacr.org/2020/1143.pdf)