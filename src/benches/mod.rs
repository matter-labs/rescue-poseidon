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

        // b.iter(|| hasher.absorb_multi(&test_inputs()))
        b.iter(|| hasher.absorb(&test_inputs()))
    }

    #[bench]
    fn bench_rescue_round_function(b: &mut Bencher) {
        let mut hasher = RescueHasher::<Bn256, 3, 2>::default();

        // b.iter(|| hasher.absorb_multi(&test_inputs()))
        b.iter(|| hasher.absorb(&test_inputs()))
    }
    #[bench]
    fn bench_rescue_prime_round_function(b: &mut Bencher) {
        let mut hasher = RescuePrimeHasher::<Bn256, 3, 2>::default();

        // b.iter(|| hasher.absorb_multi(&test_inputs()))
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

    #[bench]
    fn bench_array_static_by_ref(b: &mut Bencher){
        let rng = &mut init_rng();
        const S : usize = 20000;
        let static_vec : [Fr; S] = [Fr::rand(rng); S];

        b.iter(|| call_static_by_ref(&static_vec, rng));

    }
    #[bench]
    fn bench_array_static_by_value(b: &mut Bencher){
        let rng = &mut init_rng();
        const S : usize = 20000;
        let static_vec : [Fr; S] = [Fr::rand(rng); S];
        
        b.iter(|| call_static_by_value(static_vec, rng));
    }

    #[bench]
    fn bench_array_dynamic_by_value(b: &mut Bencher){
        let rng = &mut init_rng();
        const S : usize = 20000;
        let vec: Vec<Fr>  = (0..S).map(|_| Fr::rand(rng) ).collect();
        
        b.iter(|| call_dynamic_by_value(&vec, rng));
    }
    #[bench]
    fn bench_array_dynamic_by_ref(b: &mut Bencher){
        let rng = &mut init_rng();
        const S : usize = 20000;
        let vec: Vec<Fr>  = (0..S).map(|_| Fr::rand(rng)).collect();
        
        b.iter(|| call_dynamic_by_ref(&vec, rng));
    }


    fn call_static_by_value<R: Rng, const L: usize>(v: [Fr; L], rng: &mut R){
        let mut tmp = Fr::rand(rng);
        v.iter().for_each(|el| tmp.mul_assign(el));
    }

    fn call_static_by_ref<R: Rng, const L: usize>(v: &[Fr; L], rng: &mut R){
        let mut tmp = Fr::rand(rng);
        v.iter().for_each(|el| tmp.mul_assign(el));
    }

    fn call_dynamic_by_value<R: Rng>(v: &[Fr], rng: &mut R){
        let v = v.to_vec();
        let mut tmp = Fr::rand(rng);
        v.iter().for_each(|el| tmp.mul_assign(el));
    }

    fn call_dynamic_by_ref<R: Rng>(v: &[Fr], rng: &mut R){
        let mut tmp = Fr::rand(rng);
        v.iter().for_each(|el| tmp.mul_assign(el));
    }
}
