extern crate proc_macro;
use proc_macro::{Span, TokenStream, TokenTree};
use proc_quote::{quote, quote_spanned};
use syn::{parse_macro_input, DeriveInput, Data, Fields, Ident, Type};
use std::iter;


#[proc_macro_derive(State)]
pub fn derive_state(item : TokenStream) -> TokenStream{
    fn from_string(name : &str) -> Ident{
        Ident::new(name, proc_macro2::Span::call_site())
    }
    let x = parse_macro_input!(item as DeriveInput);

    let struct_name = x.ident;
    let (impl_generics, ty_generics, where_clause) = x.generics.split_for_impl();

    let names = match x.data{
        Data::Struct(ds) => {
            match ds.fields{
                Fields::Named(n) => {
                    n.named
                },
                _ => {
                    return err(Span::call_site(), "expected named fields, found unnamed or unit");
                }
            }
        },
        _ => {
            return err(Span::call_site(), "expected struct, found enum or union");
        },
    };

    let mut ty_pos : Option<Type> = None;
    let mut ty_vel : Option<Type> = None;

    let mut token = quote!{};

    for f in names.iter(){
        let ident = f.ident.clone().unwrap();
        if ident == from_string("mass"){
            let ty = f.ty.clone();

            token = quote!{
                #token

                impl #impl_generics Mass<#ty> for #struct_name #ty_generics #where_clause {
                    fn mass(&self) -> #ty{
                        self.mass
                    }
                }
            };
        } else if ident == from_string("charge"){
            let ty = f.ty.clone();

            token = quote!{
                #token

                impl #impl_generics Charge<#ty> for #struct_name #ty_generics #where_clause {
                    fn charge(&self) -> #ty{
                        self.charge
                    }
                }
            };
        } else if ident == from_string("diff_const"){
            let ty = f.ty.clone();

            token = quote!{
                #token

                impl #impl_generics Diffusion<#ty> for #struct_name #ty_generics #where_clause {
                    fn diff_const(&self) -> #ty{
                        self.diff_const
                    }
                }
            };
        } else if ident == from_string("orientation"){
            let ty = f.ty.clone();

            token = quote!{
                #token

                impl #impl_generics Orientation<#ty> for #struct_name #ty_generics #where_clause {
                    fn orientation(&self) -> #ty{
                        self.orientation
                    }
                }
            };
        } else if ident == from_string("pos"){
            ty_pos = Some(f.ty.clone());
        } else if ident == from_string("vel"){
            ty_vel = Some(f.ty.clone());
        }
    }

    match ty_pos{
        Some(p) => {
            match ty_vel{
                Some(v) => {
                    token = quote!{
                        #token

                        impl #impl_generics State for #struct_name #ty_generics #where_clause {
                            type Movement = (#p, #v);
                            type Position = #p;

                            fn pos(&self) -> &Self::Position{
                                &self.pos
                            }
                            fn pos_mut(&mut self) -> &mut Self::Position{
                                &mut self.pos
                            }

                            fn disp<'a>(&self, movement : &'a Self::Movement) -> &'a Self::Position{
                                &movement.0
                            }

                            fn renew_state(&mut self, movement : &Self::Movement){
                                self.pos += &movement.0;
                                self.vel += &movement.1;
                            }
                        }

                        impl #impl_generics HasVelocity for #struct_name #ty_generics #where_clause {
                            fn vel(&self) -> &<Self as State>::Position {
                                &self.vel
                            }

                            fn vel_mut(&mut self) -> &mut <Self as State>::Position {
                                &mut self.vel
                            }
                        }
                    };
                },
                None => {
                    token = quote!{
                        #token

                        impl #impl_generics State for #struct_name #ty_generics #where_clause {
                            type Movement = #p;
                            type Position = #p;

                            fn pos(&self) -> &Self::Position{
                                &self.pos
                            }
                            fn pos_mut(&mut self) -> &mut Self::Position{
                                &mut self.pos
                            }

                            fn disp<'a>(&self, movement : &'a Self::Movement) -> &'a Self::Position{
                                movement
                            }

                            fn renew_state(&mut self, movement : &Self::Movement){
                                self.pos += movement;
                            }
                        }
                    };
                }
            }
        },
        None => {
            return err(Span::call_site(), "There is no fields of name 'pos' in struct. State derive macro only work for struct has fields pos, (optional) vel, mass, charge");
        }
    }
    return token.into();
}





/// Expands into multiple `iterator.next(),` statements. E.g.
/// `many_greetings!(3);` will expand into three `println`s.
#[proc_macro]
pub fn iter_to_array(input: TokenStream) -> TokenStream {
    let tokens = input.into_iter().collect::<Vec<_>>();

    // Make sure at least one token is provided.
    if tokens.is_empty() {
        return err(Span::call_site(), "expected integer, found no input");
    }

    // Make sure we don't have too many tokens.
    // first : identifier
    // second : comma
    // third : count
    if tokens.len() > 3 {
        return err(tokens[3].span(), "unexpected second token");
    }

    // Get iterator name


    // Get the number from our token.
    let count = match &tokens[1] {
        TokenTree::Literal(lit) => {
            // Unfortunately, `Literal` doesn't have nice methods right now, so
            // the easiest way for us to get an integer out of it is to convert
            // it into string and parse it again.
            if let Ok(count) = lit.to_string().parse::<usize>() {
                count
            } else {
                let msg = format!("expected unsigned integer, found `{}`", lit);
                return err(lit.span(), msg);
            }
        }
        other => {
            let msg = format!("expected integer literal, found `{}`", other);
            return err(other.span(), msg);
        }
    };

    // Return multiple `println` statements.
    iter::repeat(quote! { ,  })
        .map(TokenStream::from)
        .take(count)
        .collect()
}

/// Report an error with the given `span` and message.
fn err(span: Span, msg: impl Into<String>) -> TokenStream {
    let msg = msg.into();
    quote_spanned!(span.into()=> {
        compile_error!(#msg);
    }).into()
}


