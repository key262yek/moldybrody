extern crate proc_macro;
use proc_macro::TokenStream;
use proc_quote::{quote, quote_spanned};
use syn::spanned::Spanned;
use syn::{parse_macro_input, DeriveInput, Data, Fields, Index};


/// Derive macros autocomplete the trait Topology.
///
/// Topology trait을 autocomplete 해주는 매크로이다.
/// Topology trait은 대응하는 Point generic type이 특정되어야하는데,
/// 이를 위해서 PointName이라는 attributes를 받는다.
#[proc_macro_derive(Topology, attributes(PointName))]
pub fn derive_topology(item: TokenStream) -> TokenStream {
    let x = parse_macro_input!(item as DeriveInput);

    let name = x.ident;
    let point_name : syn::Ident = x.attrs[0].parse_args().unwrap();

    let tokens = match x.data{
        Data::Struct(ref data) => {
            match data.fields{
                Fields::Unnamed(ref fields) =>{
                    let recurse = fields.unnamed.iter().enumerate().map(|(i, f)| {
                        let index = Index::from(i);
                        quote_spanned! {f.span()=>
                            self.#index.check_move(&pos.#index, &movement.#index)
                        }
                    });
                    quote! {
                        impl Topology<#point_name> for #name{
                            fn check_move(&self, pos : &#point_name, movement : &#point_name) -> bool{
                                false #(|| #recurse)*
                            }
                        }
                    }
                },
                _ => unimplemented!(),
            }
        },
        Data::Enum(_) |
        Data::Union(_) => unimplemented!(),
    };
    // panic!("{}", tokens.to_string());
    tokens.into()
}

/// Declarative macro to generate derive macro
#[allow(unused_macros)]
macro_rules! simple_trait{
    ($trait_name:ident, $func_name:ident) => {
        #[proc_macro_derive($trait_name)]
        pub fn $func_name(item : TokenStream) -> TokenStream{
            let x = parse_macro_input!(item as DeriveInput);

            let struct_name = x.ident;
            let (impl_generics, ty_generics, where_clause) = x.generics.split_for_impl();
            let tokens = quote!{
                impl #impl_generics $trait_name for #struct_name #ty_generics #where_clause {}
            };
            tokens.into()
        }
    }
}

simple_trait!(Point, derive_point);
simple_trait!(Force, derive_force);
simple_trait!(Displacement, derive_disp);
simple_trait!(State, derive_state);

#[proc_macro_derive(Mass)]
pub fn derive_mass(item : TokenStream) -> TokenStream{
    let x = parse_macro_input!(item as DeriveInput);

    let struct_name = x.ident;
    let (impl_generics, ty_generics, where_clause) = x.generics.split_for_impl();
    let tokens = quote!{
        impl #impl_generics Mass for #struct_name #ty_generics #where_clause {
            fn mass(&self) -> f64{
                self.mass
            }
        }
    };
    tokens.into()
}


