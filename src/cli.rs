use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_pcoa::run;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-pcoa", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Square symmetric distance-matrix TSV (skbio DistanceMatrix form); reads stdin when "-" or omitted.
    #[arg(default_value = "-")]
    input: PathBuf,

    /// Parse the input as comma-separated instead of tab-separated.
    #[arg(long, default_value_t = false)]
    csv: bool,

    /// Output path; writes stdout when "-".
    #[arg(short = 'o', long, default_value = "-")]
    output: String,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;

        let delim = if self.csv { ',' } else { '\t' };

        let reader: Box<dyn std::io::BufRead> = if self.input.as_os_str() == "-" {
            Box::new(BufReader::new(std::io::stdin().lock()))
        } else {
            Box::new(BufReader::new(File::open(&self.input).map_err(|e| {
                RsomicsError::InvalidInput(format!("{}: {e}", self.input.display()))
            })?))
        };
        let mut out: Box<dyn Write> = if self.output == "-" {
            Box::new(BufWriter::new(std::io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.output).map_err(RsomicsError::Io)?,
            ))
        };
        run(reader, &mut out, delim)?;
        out.flush().map_err(RsomicsError::Io)
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Principal Coordinates Analysis of a symmetric distance matrix.",
    origin: Some(Origin {
        upstream: "scikit-bio skbio.stats.ordination.pcoa",
        upstream_license: "BSD-3-Clause",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/biomet/53.3-4.325"),
    }),
    usage_lines: &["[distance_matrix.tsv] [-o ordination.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "csv",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: Some("false"),
                description: "Parse the distance matrix as comma-separated.",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Output path (- for stdout).",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "PCoA of a Bray-Curtis distance matrix",
            command: "rsomics-pcoa dm.tsv",
        },
        Example {
            description: "Pipe a distance matrix straight from rsomics-beta-diversity",
            command: "rsomics-beta-diversity counts.tsv | rsomics-pcoa",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
