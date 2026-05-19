pub mod fusions;
pub mod mhguide;

use regex::Regex;
use serde::{Deserialize, Deserializer};
use std::fmt::Display;
use std::str::FromStr;

pub use fusions::*;
pub use mhguide::*;

#[derive(Debug, PartialEq)]
pub enum RefGenomeVersion {
    Hg19,
    Hg38,
}

impl<'de> Deserialize<'de> for RefGenomeVersion {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value: u8 = Deserialize::deserialize(deserializer)?;
        match value {
            37 => Ok(RefGenomeVersion::Hg19),
            38 => Ok(RefGenomeVersion::Hg38),
            _ => Err(serde::de::Error::custom(format!(
                "Invalid RefGenomeVersion: {value}"
            ))),
        }
    }
}

impl Display for RefGenomeVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RefGenomeVersion::Hg19 => write!(f, "HG19"),
            RefGenomeVersion::Hg38 => write!(f, "HG38"),
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct PatientIdentifier {
    pub h_number: String,
    pub pid: String,
}

impl<'de> Deserialize<'de> for PatientIdentifier {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value: String = Deserialize::deserialize(deserializer)?;
        let parts: Vec<&str> = value.split('_').collect();
        if parts.len() != 2 {
            return Err(serde::de::Error::custom("Invalid PatientIdentifier format"));
        }
        let h_number = parts[0].to_string();
        let pid = parts[1].replace("PID", "");
        Ok(PatientIdentifier { h_number, pid })
    }
}

#[allow(clippy::expect_used)]
pub fn three_letter_protein_modification(input: &str) -> String {
    let mapping = [
        ("A", "Ala"),
        ("C", "Cys"),
        ("G", "Gly"),
        ("I", "Ile"),
        ("L", "Leu"),
        ("M", "Met"),
        ("P", "Pro"),
        ("S", "Ser"),
        ("T", "Thr"),
        ("V", "Val"),
        ("F", "Phe"),
        ("Y", "Tyr"),
        ("W", "Trp"),
        ("H", "His"),
        ("Q", "Gln"),
        ("R", "Arg"),
        ("N", "Asn"),
        ("K", "Lys"),
        ("D", "Asp"),
        ("E", "Glu"),
    ];

    let three_letter_codes = Regex::from_str(
        r"(Phe|Leu|Ser|Tyr|Cys|Trp|Pro|His|Gln|Arg|Ile|Met|Thr|Asn|Lys|Val|Ala|Asp|Glu|Gly)",
    )
    .expect("Invalid regex");

    if three_letter_codes.is_match(input) {
        return input.to_string();
    }

    let mut result = input.to_string();
    for (old, new) in mapping {
        result = result.replace(old, new);
    }

    result
}

#[derive(Debug, Default, PartialEq)]
pub struct DnaChange {
    pub start: String,
    pub end: String,

    pub ref_allele: String,
    pub alt_allele: String,
}

impl FromStr for DnaChange {
    type Err = String;

    #[allow(clippy::expect_used)]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let regexes = [
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?P<ref>[ACGT])>(?P<alt>[ACGT])$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?del$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)(?:_(?P<end>-?\d+))?dup$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)_-?(?P<end>-?\d+)ins(?P<alt>[ACGT]+)$")
                .expect("Invalid regex"),
            Regex::new(r"(?P<type>[cg])\.(?P<start>-?\d+)_-?(?P<end>-?\d+)delins(?P<alt>[ACGT]+)$")
                .expect("Invalid regex"),
        ];

        for regex in &regexes {
            if let Some(captures) = regex.captures(s) {
                let start = captures["start"].parse::<i128>().unwrap_or_default();
                let end = captures
                    .name("end")
                    .map_or(0, |m| m.as_str().parse::<i128>().unwrap_or_default());
                let ref_allele = captures
                    .name("ref")
                    .map(|m| m.as_str().into())
                    .unwrap_or_default();
                let alt_allele = captures
                    .name("alt")
                    .map(|m| m.as_str().into())
                    .unwrap_or_default();

                let start = if start == 0 {
                    String::new()
                } else {
                    start.to_string()
                };
                let end = if end == 0 {
                    String::new()
                } else {
                    end.to_string()
                };
                return Ok(DnaChange {
                    start,
                    end,
                    ref_allele,
                    alt_allele,
                });
            }
        }
        Err("Invalid DNA change format".to_string())
    }
}
