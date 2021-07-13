extern crate proc_macro;
use proc_macro::TokenStream;
use proc_quote::{quote, quote_spanned};
use syn::spanned::Spanned;
use syn::{parse_macro_input, DeriveInput, Data, Fields, Index};

/// Derive macros autocomplete the trait Point.
///
/// MoldyBrody에서는 Point trait이 존재하는데, 이를 자동으로 구성해주는 매크로이다.
#[proc_macro_derive(Point)]
pub fn derive_point(item: TokenStream) -> TokenStream {
    let x = syn::parse_macro_input!(item as DeriveInput);
    let name = x.ident;

    let tokens = proc_quote::quote!{
        impl Point for #name {}
    };
    tokens.into()
}


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


