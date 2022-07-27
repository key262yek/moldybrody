// use clap::builder::Command;
use clap::Arg;

pub trait CommandBuilder<'h, const N: usize> {
    const SUBCOMMAND: &'h str;
    const NARGS: usize = N;

    // fn command() -> Command<'h>;
    fn args() -> [Arg<'h>; N];
}
