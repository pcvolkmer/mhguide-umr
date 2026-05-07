pub mod fusions;
pub mod mhguide;

use regex::Regex;
use serde::{Deserialize, Deserializer};
use std::fmt::Display;
use std::str::FromStr;

pub use mhguide::*;
pub use fusions::*;

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
        let pid = parts[1].to_string();
        Ok(PatientIdentifier { h_number, pid })
    }
}

#[allow(clippy::expect_used)]
pub fn three_letter_protein_modification(short: &str) -> String {
    fn map_value(value: &str) -> String {
        match value {
            "*" => "*",
            "=" => "=",
            "fs" => "fs",
            "F" => "Phe",
            "L" => "Leu",
            "S" => "Ser",
            "Y" => "Tyr",
            "C" => "Cys",
            "W" => "Trp",
            "P" => "Pro",
            "H" => "His",
            "Q" => "Gln",
            "R" => "Arg",
            "I" => "Ile",
            "M" => "Met",
            "T" => "Thr",
            "N" => "Asn",
            "K" => "Lys",
            "V" => "Val",
            "A" => "Ala",
            "D" => "Asp",
            "E" => "Glu",
            "G" => "Gly",
            _ => value,
        }
            .to_string()
    }

    let regex = Regex::new(r"^p\.(?<refA>[*FLSYCWPHQRIMTNKVADEG])?(?<posA>\d+)?(?<sep>_)?(?<refB>[*FLSYCWPHQRIMTNKVADEG])(?<posB>\d+)(?<type>del|ins|delins|dup)?(?<alt>[*=FLSYCWPHQRIMTNKVADEG]+|fs)?$")
        .expect("Invalid regex");

    if let Some(captures) = regex.captures(short) {
        let refa_capture = match captures.name("refA") {
            Some(m) => m.as_str(),
            None => "",
        };
        let posa_capture = match captures.name("posA") {
            Some(m) => m.as_str(),
            None => "",
        };
        let sep_capture = match captures.name("sep") {
            Some(m) => m.as_str(),
            None => "",
        };
        let refb_capture = match captures.name("refB") {
            Some(m) => m.as_str(),
            None => "",
        };
        let posb_capture = match captures.name("posB") {
            Some(m) => m.as_str(),
            None => "",
        };
        let type_capture = match captures.name("type") {
            Some(m) => m.as_str(),
            None => "",
        };
        let alt_capture = match captures.name("alt") {
            Some(m) => m
                .as_str()
                .chars()
                .collect::<Vec<_>>()
                .iter()
                .map(|c| map_value(&c.to_string()))
                .collect::<String>(),
            None => String::new(),
        };
        return format!(
            "p.{}{}{}{}{}{}{}",
            map_value(refa_capture),
            posa_capture,
            map_value(sep_capture),
            map_value(refb_capture),
            posb_capture,
            type_capture,
            alt_capture
        );
    }

    short.to_string()
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



