#[macro_export]
macro_rules! sponge_impl {
    ($hasher_name:ty) => {
        impl<E: Engine, const S: usize, const R: usize> SpongeState<E, S> for $hasher_name {

            fn state_as_ref(&self) -> &[E::Fr; S] {
                &self.state
            }

            fn state_as_mut(&mut self) -> &mut [E::Fr;S] {
                &mut self.state
            }           
        }

        impl<E: Engine, const S: usize, const R: usize> StatefulSponge<E, S, R> for $hasher_name {}

        impl<E: Engine, const S: usize, const R: usize> SpongeMode<E> for $hasher_name {
            fn get_mode(&self) -> SpongeModes{
                self.sponge_mode.to_owned()
            }
            fn update_mode(&mut self, mode: SpongeModes){
                self.sponge_mode = mode;
            }
        }
    };
}
