#[macro_export]
macro_rules! sponge_impl {
    ($hasher_name:ty) => {
    // ($hasher_name:ty, $hasher_params:expr) => {
    // ($hasher_name:ty, $sponge_type:expr) => {
        impl<E: Engine> SpongeState<E> for $hasher_name {

            fn state_as_ref(&self) -> &[E::Fr] {
                self.state.as_ref()
            }
            fn storage_as_ref(&self) -> &[E::Fr] {
                self.tmp_storage.as_ref()
            }            
            fn state_as_mut(&mut self) -> &mut [E::Fr] {
                self.state.as_mut()
            }
            fn storage_as_mut(&mut self) -> &mut Vec<E::Fr> {
                self.tmp_storage.as_mut()
            }
        }

        impl<E: Engine> SpongeParams for $hasher_name {
            fn rate(&self) -> usize {
                self.params.rate
            }
        }

        impl<E: Engine> StatefulSponge<E> for $hasher_name {}
        
        impl<E: Engine> SpongeMode<E> for $hasher_name {   
            fn get_mode(&self) -> SpongeModes<E>{
                self.sponge_mode.to_owned()
            }         
            fn update_mode(&mut self, mode: SpongeModes<E>){
                self.sponge_mode = mode;
            }         
        }
    };
}