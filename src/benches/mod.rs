#[cfg(test)]
mod test {
    use crate::{poseidon::PoseidonHasher, rescue::RescueHasher, rescue_prime::RescuePrimeHasher};
    use crate::{sponge::StatefulSponge, tests::init_rng};
    use franklin_crypto::bellman::{Field, pairing::bn256::{Bn256, Fr}};
    use rand::{Rand, Rng};
    use test::Bencher;

    fn test_inputs() -> Vec<Fr> {
        let rng = &mut init_rng();
        let mut inputs = vec![];
        for _ in 0..2 {
            inputs.push(Fr::rand(rng));
        }
        inputs
    }

    #[bench]
    fn bench_poseidon_round_function(b: &mut Bencher) {
        let mut hasher = PoseidonHasher::<Bn256, 3, 2>::default();

        b.iter(|| hasher.absorb(&test_inputs()))
    }

    #[bench]
    fn bench_rescue_round_function(b: &mut Bencher) {
        let mut hasher = RescueHasher::<Bn256, 3, 2>::default();

        b.iter(|| hasher.absorb(&test_inputs()))
    }
    #[bench]
    #[should_panic] // TODO
    fn bench_rescue_prime_round_function(b: &mut Bencher) {
        let mut hasher = RescuePrimeHasher::<Bn256, 3, 2>::default();

        b.iter(|| hasher.absorb(&test_inputs()))
    }

    #[bench]
    fn bench_old_rescue_round_function(b: &mut Bencher) {
        use franklin_crypto::rescue::bn256::Bn256RescueParams;
        use franklin_crypto::rescue::StatefulRescue;

        let params = Bn256RescueParams::new_checked_2_into_1();

        let mut rescue = StatefulRescue::<Bn256>::new(&params);

        b.iter(|| rescue.absorb(&test_inputs()))
    }

    #[bench]
    fn bench_old_poseidon_round_function(b: &mut Bencher) {
        use poseidon_hash::bn256::Bn256PoseidonParams;
        use poseidon_hash::StatefulSponge;

        let params = Bn256PoseidonParams::new_checked_2_into_1();

        let mut poseidon = StatefulSponge::<Bn256>::new(&params);

        b.iter(|| poseidon.absorb(&test_inputs()))
    }
}
