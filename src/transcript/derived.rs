#[macro_export]
macro_rules! stateful_transcript {
    ($transcrit_name:ty, $hasher_path:expr) => {
        impl<E: Engine> Prng<E::Fr> for $transcrit_name {
            type Input = E::Fr;

            type InitializationParameters = ();

            fn new() -> Self {
                Self {
                    sponge: $hasher_path(),
                }
            }

            fn commit_input(&mut self, input: &Self::Input) {
                self.sponge.absorb(*input);
            }

            fn get_challenge(&mut self) -> E::Fr {
                self.sponge.squeeze_single()
            }
        }

        impl<E: Engine> Transcript<E::Fr> for $transcrit_name {
            fn commit_bytes(&mut self, _: &[u8]) {
                unimplemented!()
            }

            fn commit_field_element(&mut self, element: &E::Fr) {
                self.commit_input(element)
            }

            fn get_challenge_bytes(&mut self) -> Vec<u8> {
                let mut buf = vec![];
                let fe = self.get_challenge();
                let fe_as_repr = fe.into_repr();
                fe_as_repr.write_le(&mut buf).expect("filled with bytes");

                buf
            }

            fn commit_fe<FF: PrimeField>(&mut self, _: &FF) {
                unimplemented!()
            }
        }
    };
}
