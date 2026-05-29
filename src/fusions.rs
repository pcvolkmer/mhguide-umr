use regex::Regex;
use std::fmt::Display;
use std::str::FromStr;

#[derive(Debug, PartialEq)]
pub enum Fusion {
    RnaFusion {
        partner_3: String,
        partner_5: String,
        transcript_id_3: String,
        transcript_id_5: String,
        transcript_position_3: u32,
        transcript_position_5: u32,
        exon_id_3: String,
        exon_id_5: String,
        strand_3: String,
        strand_5: String,
        number_reported_reads: u32,
    },
    #[allow(unused)]
    DnaFusion {
        partner_3: String,
        partner_5: String,
        chromosome_3: String,
        chromosome_5: String,
        transcript_position_3: u32,
        transcript_position_5: u32,
    },
}

impl Display for Fusion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Fusion::RnaFusion { .. } => {
                write!(f, "RNA Fusion")
            }
            Fusion::DnaFusion { .. } => {
                write!(f, "DNA Fusion")
            }
        }
    }
}

impl FromStr for Fusion {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let rna_regex = Regex::new(r"(?<partner_5>[A-Z0-9_\\-]+)\(ex (?<exon_5>\d+)\)::(?<partner_3>[A-Z0-9_\\-]+)\(ex (?<exon_3>\d+)\)[,;]\s[Tt]ranscript\sID:\s(?<transcript_id_5>NM_\d+\.\d+)/(?<transcript_id_3>NM_\d+\.\d+)[,;]\s([Ss]trand:\s(?<strand_5>[+-]?)/(?<strand_3>[+-]?)[,;]\s)?[Bb]reakpoint:\schr\d+:(?<transcript_position_5>\d+)/chr\d+:(?<transcript_position_3>\d+)[,;]\s[Ss]upporting\sread\spairs:\s(?<number_reported_reads>\d+)").map_err(|_| ())?;

        match rna_regex.captures(s) {
            Some(captures) => {
                let partner_3 = match captures.name("partner_3") {
                    Some(value) => value.as_str().to_owned(),
                    _ => return Err(()),
                };
                let partner_5 = match captures.name("partner_5") {
                    Some(value) => value.as_str().to_owned(),
                    _ => return Err(()),
                };
                let transcript_id_3 = match captures.name("transcript_id_3") {
                    Some(value) => value.as_str().to_owned(),
                    _ => return Err(()),
                };
                let transcript_id_5 = match captures.name("transcript_id_5") {
                    Some(value) => value.as_str().to_owned(),
                    _ => return Err(()),
                };
                let transcript_position_3 = match captures.name("transcript_position_3") {
                    Some(value) => match value.as_str().parse::<u32>() {
                        Ok(value) => value,
                        Err(_) => return Err(()),
                    },
                    _ => return Err(()),
                };
                let transcript_position_5 = match captures.name("transcript_position_5") {
                    Some(value) => match value.as_str().parse::<u32>() {
                        Ok(value) => value,
                        Err(_) => return Err(()),
                    },
                    _ => return Err(()),
                };
                let exon_id_3 = match captures.name("exon_3") {
                    Some(value) => value.as_str().to_string(),
                    _ => return Err(()),
                };
                let exon_id_5 = match captures.name("exon_5") {
                    Some(value) => value.as_str().to_string(),
                    _ => return Err(()),
                };
                let strand_3 = match captures.name("strand_3") {
                    Some(value) => value.as_str().to_owned(),
                    _ => String::new(),
                };
                let strand_5 = match captures.name("strand_5") {
                    Some(value) => value.as_str().to_owned(),
                    _ => String::new(),
                };
                let number_reported_reads = match captures.name("number_reported_reads") {
                    Some(value) => match value.as_str().parse::<u32>() {
                        Ok(value) => value,
                        _ => return Err(()),
                    },
                    _ => return Err(()),
                };

                Ok(Fusion::RnaFusion {
                    partner_3,
                    partner_5,
                    transcript_id_3,
                    transcript_id_5,
                    transcript_position_3,
                    transcript_position_5,
                    exon_id_3,
                    exon_id_5,
                    strand_3,
                    strand_5,
                    number_reported_reads,
                })
            }
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::fusions::Fusion;
    use crate::fusions::Fusion::RnaFusion;
    use crate::mhguide::MhGuide;
    use rstest::rstest;
    use std::str::FromStr;

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_fusion_examples_from_report_specification() {
        static INPUT: &str = "AKAP8L(ex 12)::BRD4(ex 2), \
            transcript ID: NM_014371.4/NM_014299.2, \
            breakpoint: chr19:15507961/chr19:15383944, \
            supporting read pairs: 1158";

        let value = Fusion::from_str(INPUT).unwrap();
        assert_eq!(
            value,
            RnaFusion {
                partner_3: "BRD4".to_string(),
                partner_5: "AKAP8L".to_string(),
                transcript_id_3: "NM_014299.2".to_string(),
                transcript_id_5: "NM_014371.4".to_string(),
                transcript_position_3: 15383944,
                transcript_position_5: 15507961,
                exon_id_3: "2".to_string(),
                exon_id_5: "12".to_string(),
                strand_3: String::new(),
                strand_5: String::new(),
                number_reported_reads: 1158,
            }
        );
    }

    #[rstest]
    #[case(
        "ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Strand: -/-; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234",
        RnaFusion {
                partner_3: "ABCD2".to_string(),
                partner_5: "ABCD1".to_string(),
                transcript_id_3: "NM_012456.2".to_string(),
                transcript_id_5: "NM_012345.4".to_string(),
                transcript_position_3: 13456789,
                transcript_position_5: 12345678,
                exon_id_3: "2".to_string(),
                exon_id_5: "1".to_string(),
                strand_3: "-".to_string(),
                strand_5: "-".to_string(),
                number_reported_reads: 1234,
            }
    )]
    #[case(
        "ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Strand: -/+; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234",
        RnaFusion {
                partner_3: "ABCD2".to_string(),
                partner_5: "ABCD1".to_string(),
                transcript_id_3: "NM_012456.2".to_string(),
                transcript_id_5: "NM_012345.4".to_string(),
                transcript_position_3: 13456789,
                transcript_position_5: 12345678,
                exon_id_3: "2".to_string(),
                exon_id_5: "1".to_string(),
                strand_3: "+".to_string(),
                strand_5: "-".to_string(),
                number_reported_reads: 1234,
            }
    )]
    #[case(
        "ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Strand: +/-; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234",
        RnaFusion {
                partner_3: "ABCD2".to_string(),
                partner_5: "ABCD1".to_string(),
                transcript_id_3: "NM_012456.2".to_string(),
                transcript_id_5: "NM_012345.4".to_string(),
                transcript_position_3: 13456789,
                transcript_position_5: 12345678,
                exon_id_3: "2".to_string(),
                exon_id_5: "1".to_string(),
                strand_3: "-".to_string(),
                strand_5: "+".to_string(),
                number_reported_reads: 1234,
            }
    )]
    #[case(
        "ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234",
        RnaFusion {
                partner_3: "ABCD2".to_string(),
                partner_5: "ABCD1".to_string(),
                transcript_id_3: "NM_012456.2".to_string(),
                transcript_id_5: "NM_012345.4".to_string(),
                transcript_position_3: 13456789,
                transcript_position_5: 12345678,
                exon_id_3: "2".to_string(),
                exon_id_5: "1".to_string(),
                strand_3: String::new(),
                strand_5: String::new(),
                number_reported_reads: 1234,
            }
    )]
    #[case(
        "ABCD1(ex 1)::ABCD2(ex 2); Transcript ID: NM_012345.4/NM_012456.2; Strand: /; Breakpoint: chr19:12345678/chr19:13456789; Supporting read pairs: 1234",
        RnaFusion {
                partner_3: "ABCD2".to_string(),
                partner_5: "ABCD1".to_string(),
                transcript_id_3: "NM_012456.2".to_string(),
                transcript_id_5: "NM_012345.4".to_string(),
                transcript_position_3: 13456789,
                transcript_position_5: 12345678,
                exon_id_3: "2".to_string(),
                exon_id_5: "1".to_string(),
                strand_3: String::new(),
                strand_5: String::new(),
                number_reported_reads: 1234,
            }
    )]
    #[allow(clippy::expect_used)]
    fn test_extract_rna_fusion_from_string(#[case] input: &str, #[case] expected: Fusion) {
        let value = Fusion::from_str(input).expect("Failed to parse fusion string");
        assert_eq!(value, expected);
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    #[allow(clippy::panic)]
    fn test_extract_rna_fusion_from_report() {
        static MHGUIDE: &str = include_str!("../testfiles/rnafusion-mhguide.json");

        let value = serde_json::from_str::<MhGuide>(MHGUIDE).unwrap().fusions();
        assert_eq!(value.len(), 1);

        match value.first() {
            Some(fusion) => {
                assert_eq!(
                    fusion,
                    &RnaFusion {
                        partner_3: "ABCD2".to_string(),
                        partner_5: "ABCD1".to_string(),
                        transcript_id_3: "NM_012456.2".to_string(),
                        transcript_id_5: "NM_012345.4".to_string(),
                        transcript_position_3: 13456789,
                        transcript_position_5: 12345678,
                        exon_id_3: "2".to_string(),
                        exon_id_5: "1".to_string(),
                        strand_3: "-".to_string(),
                        strand_5: "-".to_string(),
                        number_reported_reads: 1234,
                    }
                );
            }
            _ => panic!("No RNA fusion found"),
        };
    }
}
