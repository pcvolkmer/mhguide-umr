use crate::fusions::Fusion;
use crate::{DnaChange, PatientIdentifier, RefGenomeVersion};
use rayon::prelude::*;
use regex::Regex;
use serde::{Deserialize, Deserializer};
use std::fmt::Display;
use std::str::FromStr;

#[derive(Debug, Deserialize, PartialEq)]
pub struct MhGuide {
    #[serde(rename = "GENERAL")]
    pub general: General,
    #[serde(rename = "VARIANT_LONG_LIST")]
    pub variants: Vec<Variant>,
    #[serde(rename = "BIOMARKERS")]
    pub biomarkers: Biomarkers,
    #[serde(rename = "REPORT_NARRATIVE")]
    pub report_narrative: String,
}

impl MhGuide {
    /// Retrieves all the variants.
    ///
    /// This method collects all variants using parallel iteration, enabling efficient processing for large datasets.
    ///
    /// # Returns
    ///
    /// A `Vec` containing references to all the `Variant` instances in `self.variants`.
    ///
    /// # Example
    /// ```rust
    /// let result = mh_guide.all_variants();
    /// for variant in result {
    ///     println!("{:?}", variant);
    /// }
    /// ```
    pub fn all_variants(&self) -> Vec<&Variant> {
        self.variants.par_iter().collect()
    }

    /// ```rust
    /// Retrieves a list of relevant genetic variants based on oncogenic properties and matching criteria.
    ///
    /// This function collects the relevant variants for a report by:
    /// 1. Starting with a collection of oncogenic variants.
    /// 2. Adding CNA variants from Biomarkers
    /// 3. Adding variants mentioned in `REPORT_NARRATIVE` without being artifacts.
    /// 4. Removing variants that are mentioned as "Artifacts" if requested
    ///
    /// The function ensures that the resulting list is deduplicated before being returned.
    ///
    /// # Returns
    /// A vector of references to relevant `Variant` objects.
    ///
    /// # Examples
    /// ```rust
    /// let variants = mh_guide.relevant_variants();
    /// for variant in variants {
    ///     println!("{:?}", variant);
    /// }
    /// ```
    pub fn relevant_variants(&self, no_artifacts: bool) -> Vec<&Variant> {
        let mut result = self.oncogenic_variants();

        let cnv_biomarker_variant_ids = self
            .biomarkers
            .notable_biomarkers
            .iter()
            .flat_map(|nb| nb.biomarkers.iter())
            .filter(|b| {
                b.variant_effect == Some(VariantEffect::CopyGain)
                    || b.variant_effect == Some(VariantEffect::CopyLoss)
            })
            .map(|b| b.id.to_owned())
            .collect::<Vec<_>>();

        result.extend(
            self.variants
                .par_iter()
                .filter(|v| matches!(v.display_variant_type, Some(ResultType::CopyNumberVariant)))
                .filter(|v| cnv_biomarker_variant_ids.contains(&v.id))
                .collect::<Vec<_>>(),
        );

        let report_narrative_simple_variants = self.report_narrative_simple_variants();

        result.extend(
            self.variants
                .par_iter()
                .filter(|v| match v.protein_modification {
                    Some(ref protein_modification) => {
                        report_narrative_simple_variants
                            .iter()
                            .any(|(gene, modification)| {
                                modification.starts_with("p.")
                                    && gene.clone() == v.gene_symbol.clone().unwrap_or_default()
                                    && modification == protein_modification
                            })
                    }
                    _ => false,
                })
                .collect::<Vec<_>>(),
        );

        result.extend(
            self.variants
                .par_iter()
                .filter(|v| match v.transcript_hgvs_modified_object {
                    Some(ref transcript_hgvs_modified_object) => report_narrative_simple_variants
                        .iter()
                        .any(|(gene, modification)| {
                            modification.starts_with("c.")
                                && gene.clone() == v.gene_symbol.clone().unwrap_or_default()
                                && modification == transcript_hgvs_modified_object
                        }),
                    _ => false,
                })
                .collect::<Vec<_>>(),
        );

        result.extend(
            self.variants
                .par_iter()
                .filter(|v| match v.protein_variant_type {
                    Some(ref protein_variant_type) => self
                        .report_narrative_copy_variants()
                        .iter()
                        .any(|(gene, gnc)| {
                            gene.clone() == v.gene_symbol.clone().unwrap_or_default()
                                && protein_variant_type == &ResultType::CopyNumberVariant
                                && match v.copy_number {
                                    Some(copy_number) => copy_number.eq(gnc),
                                    None => false,
                                }
                        }),
                    _ => false,
                })
                .collect::<Vec<_>>(),
        );

        if no_artifacts {
            result.retain(|v| {
                !self
                    .removable_report_narrative_variants()
                    .iter()
                    .any(|(gene, modification)| {
                        gene == &v.gene_symbol.clone().unwrap_or_default()
                            && modification.starts_with("p.")
                            && modification == &v.protein_modification.clone().unwrap_or_default()
                    })
            });

            result.retain(|v| {
                !self
                    .removable_report_narrative_variants()
                    .iter()
                    .any(|(gene, modification)| {
                        gene == &v.gene_symbol.clone().unwrap_or_default()
                            && modification.starts_with("c.")
                            && modification
                                == &v
                                    .transcript_hgvs_modified_object
                                    .clone()
                                    .unwrap_or_default()
                    })
            });
        }

        result.sort_by_key(|v| v.protein_modification.clone());
        result.sort_by_key(|v| v.transcript_hgvs_modified_object.clone());
        result.sort_by_key(|v| v.gene_symbol.clone());
        result.dedup();

        result
    }

    /// ```rust
    /// Filters and retrieves all oncogenic variants
    ///
    /// This function searches through the `variants` field and identifies those
    /// variants that are classified as "oncogenic". The search is case-insensitive
    /// and checks if the `ONCOGENIC_CLASSIFICATION_NAME` contains the substring
    /// "oncogenic".
    ///
    /// # Returns
    ///
    /// A `Vec` containing references to the variants that are classified as oncogenic.
    ///
    /// # Example
    ///
    /// ```rust
    /// let variants = my_struct.oncogenic_variants();
    ///
    /// for variant in variants {
    ///     println!("Oncogenic Variant: {:?}", variant);
    /// }
    /// ```
    pub fn oncogenic_variants(&self) -> Vec<&Variant> {
        self.variants
            .par_iter()
            .filter(|v| match v.oncogenic_classification_name {
                Some(ref name) => name.to_ascii_lowercase().contains("oncogenic"),
                _ => false,
            })
            .collect()
    }

    /// Retrieves the Tumor Mutational Burden (TMB) value from the notable biomarkers.
    ///
    /// This function iterates through the notable biomarkers and their corresponding
    /// biomarkers to find a biomarker with a `display_variant_type` equal to `VariantType::TMB`.
    /// If a matching biomarker is found, its `TMB_VARIANT_COUNT_PER_MEGABASE` value is fetched,
    /// parsed from a `String` to an `f32`, and returned as an `Option<f32>`.
    ///
    /// ### Returns
    /// - `Some(f32)` if a valid TMB value is found and successfully parsed.
    /// - `None` if no matching TMB variant is found, the `TMB_VARIANT_COUNT_PER_MEGABASE` is missing,
    ///   or its value cannot be parsed into an `f32`.
    ///
    /// ### Example
    /// ```rust
    /// let tmb_value = some_instance.tmb_value();
    /// match tmb_value {
    ///     Some(value) => println!("TMB Value: {}", value),
    ///     None => println!("No TMB Value found."),
    /// }
    /// ```
    pub fn tmb_value(&self) -> Option<f32> {
        for notable_biomarker in &self.biomarkers.notable_biomarkers {
            for biomarker in &notable_biomarker.biomarkers {
                let Some(display_variant_type) = &biomarker.display_variant_type else {
                    return None;
                };
                if display_variant_type == &ResultType::TMB {
                    let score = match biomarker.tmb_variant_count_per_megabase {
                        Some(ref value) => value.clone(),
                        None => String::new(),
                    };
                    return f32::from_str(&score).ok();
                }
            }
        }
        None
    }

    /// Retrieves the Homologous recombination deficiency (HRD) value from the notable biomarkers.
    ///
    /// This function iterates through the notable biomarkers and their corresponding
    /// biomarkers to find a biomarker with a `DISPLAY_VARIANT_TYPE` equal to `VariantType::HRD`.
    /// If a matching biomarker is found, its `SCORE` value is fetched,
    /// parsed from a `String` to an `f32`, and returned as an `Option<f32>`.
    ///
    /// ### Returns
    /// - `Some(f32)` if a valid TMB value is found and successfully parsed.
    /// - `None` if no matching TMB variant is found, the `SCORE` is missing,
    ///   or its value cannot be parsed into an `f32`.
    ///
    /// ### Example
    /// ```rust
    /// let hrd_score = some_instance.hrd_score();
    /// match hrd_score {
    ///     Some(value) => println!("HRD Score: {}", value),
    ///     None => println!("No HRD Score found."),
    /// }
    /// ```
    pub fn hrd_score(&self) -> Option<f32> {
        self.biomarker_score_value(&ResultType::HRD)
    }

    /// Retrieves the Microsatellite Instability (MSI) value from the notable biomarkers.
    ///
    /// This function iterates through the notable biomarkers and their corresponding
    /// biomarkers to find a biomarker with a `DISPLAY_VARIANT_TYPE` equal to `VariantType::MSI`.
    /// If a matching biomarker is found, its `SCORE` value is fetched,
    /// parsed from a `String` to an `f32`, and returned as an `Option<f32>`.
    ///
    /// ### Returns
    /// - `Some(f32)` if a valid TMB value is found and successfully parsed.
    /// - `None` if no matching TMB variant is found, the `SCORE` is missing,
    ///   or its value cannot be parsed into an `f32`.
    ///
    /// ### Example
    /// ```rust
    /// let msi_score = some_instance.msi_score();
    /// match msi_score {
    ///     Some(value) => println!("MSI Score: {}", value),
    ///     None => println!("No MSI Score found."),
    /// }
    /// ```
    pub fn msi_score(&self) -> Option<f32> {
        self.biomarker_score_value(&ResultType::MSI)
    }

    /// Extracts a list of `Fusion` objects from the `report_narrative` text.
    ///
    /// This method processes the `report_narrative` by splitting its contents
    /// into individual lines, attempting to parse each line into a `Fusion` instance,
    /// and collecting successfully parsed `Fusion` objects into a vector.
    ///
    /// # Returns
    ///
    /// A `Vec<Fusion>` containing all `Fusion` objects that were successfully
    /// parsed from the `report_narrative`.
    ///
    /// # Example
    ///
    /// ```rust
    /// let result = instance.fusions();
    /// for fusion in result {
    ///     println!("{:?}", fusion);
    /// }
    /// ```
    pub fn fusions(&self) -> Vec<Fusion> {
        self.report_narrative
            .split('\n')
            .filter_map(|line| Fusion::from_str(line).ok())
            .collect::<Vec<_>>()
    }

    fn biomarker_score_value(&self, variant_type: &ResultType) -> Option<f32> {
        for notable_biomarker in &self.biomarkers.notable_biomarkers {
            for biomarker in &notable_biomarker.biomarkers {
                let Some(display_variant_type) = &biomarker.display_variant_type else {
                    return None;
                };
                if display_variant_type == variant_type {
                    let score = match biomarker.score {
                        Some(ref score) => score.clone(),
                        None => String::new(),
                    };
                    return f32::from_str(&score).ok();
                }
            }
        }
        None
    }

    fn report_narrative_simple_variants(&self) -> Vec<(String, String)> {
        let mut result = Self::find_report_narrative_simple_variants(&self.report_narrative);
        let removable = self.removable_report_narrative_variants();

        result.retain(|(gene, modification)| {
            !removable
                .iter()
                .any(|(g, m)| g == gene && m == modification)
        });

        result
    }

    #[allow(clippy::expect_used)]
    fn removable_report_narrative_variants(&self) -> Vec<(String, String)> {
        self.report_narrative
            .split('\n')
            // Only one variant per line
            .filter(|&s| Self::find_report_narrative_simple_variants(s).len() == 1)
            // Exclusion string(s)
            .filter(|&s| {
                (s.contains("möglich") || s.contains("wahrscheinlich")) && s.contains("Artefakt")
            })
            .flat_map(Self::find_report_narrative_simple_variants)
            .collect()
    }

    #[allow(clippy::expect_used)]
    fn find_report_narrative_simple_variants(s: &str) -> Vec<(String, String)> {
        fn collect(s: &str, re: &Regex) -> Vec<(String, String)> {
            re.find_iter(s)
                .filter_map(|m| {
                    let parts = m.as_str().split(' ').collect::<Vec<_>>();
                    if parts.len() == 2 {
                        Some((parts[0].trim().to_owned(), parts[1].trim().to_owned()))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        }

        let protein_regex = Regex::new(
            r"[A-Z0-9_\\-]+\s+p\.[*FLSYCWPHQRIMTNKVADEG]?(\d+)?_?[*FLSYCWPHQRIMTNKVADEG]\d+(del|ins|delins|dup)?(\*|=|fs|[*=FLSYCWPHQRIMTNKVADEG]+)?"
        )
            .expect("Invalid regex");
        let mut result = collect(s, &protein_regex);

        let cdna_regex = Regex::new(
            r"[A-Z0-9_\\-]+\s+c\.(-?\d+)(?:_(-?\d+))?([ACGT]>|dup|del|ins|delins)([ACGT]+)?",
        )
        .expect("Invalid regex");
        let cdna_result = collect(s, &cdna_regex);
        result.extend(cdna_result);

        result
    }

    #[allow(clippy::expect_used)]
    fn report_narrative_copy_variants(&self) -> Vec<(String, f32)> {
        let regex = Regex::new(r"(?<gene>[A-Z0-9_\\-]+)\s*.*GCN\s*[:=]\s*(?<gcn>\d+\.\d+)")
            .expect("Invalid regex");

        self.report_narrative
            .split('\n')
            .filter_map(|line| {
                let captures = regex.captures(line)?;
                let gene = captures.name("gene");
                let gcn = captures.name("gcn");
                if gene.is_none() || gcn.is_none() {
                    return None;
                }
                let gene = gene.expect("Missing gene").as_str().to_owned();
                let gcn = gcn
                    .expect("Missing GNC")
                    .as_str()
                    .parse::<f32>()
                    .unwrap_or_default();
                Some((gene, gcn))
            })
            .collect::<Vec<_>>()
    }
}

#[derive(Debug, Deserialize, PartialEq)]
pub struct General {
    #[serde(rename = "ORDER_DATE")]
    pub order_date: String,
    #[serde(rename = "REF_GENOME_VERSION")]
    pub ref_genome_version: RefGenomeVersion,
    #[serde(rename = "PATIENT_IDENTIFIER")]
    pub patient_identifier: PatientIdentifier,
}

#[derive(Debug, PartialEq)]
pub enum VariantEffect {
    CopyGain,
    CopyLoss,
    Other(String),
}

impl<'de> Deserialize<'de> for VariantEffect {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        match Deserialize::deserialize(deserializer)? {
            "Copy gain" => Ok(Self::CopyGain),
            "Copy loss" => Ok(Self::CopyLoss),
            other => Ok(Self::Other(other.to_string())),
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum ResultType {
    SimpleVariant(String),
    CopyNumberVariant,
    DnaFusion,
    RnaFusion,
    TMB,
    HRD,
    MSI,
    Other(String),
}

impl<'de> Deserialize<'de> for ResultType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        match Deserialize::deserialize(deserializer)? {
            "SNV" => Ok(Self::SimpleVariant("SNV".to_string())),
            "ins" => Ok(Self::SimpleVariant("ins".to_string())),
            "del" => Ok(Self::SimpleVariant("del".to_string())),
            "CNA" => Ok(Self::CopyNumberVariant),
            "DNA fusion" => Ok(Self::DnaFusion),
            "RNA fusion" => Ok(Self::RnaFusion),
            "TMB" => Ok(Self::TMB),
            "HRD" => Ok(Self::HRD),
            "MSI" => Ok(Self::MSI),
            other => Ok(Self::Other(other.to_string())),
        }
    }
}

impl Default for ResultType {
    fn default() -> Self {
        Self::Other("nicht angegeben".to_string())
    }
}

impl Display for ResultType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SimpleVariant(variant_type) => write!(f, "Einfache Variante ({variant_type})"),
            Self::CopyNumberVariant => write!(f, "Copy Number Variation"),
            Self::DnaFusion => write!(f, "DNA Fusion"),
            Self::RnaFusion => write!(f, "RNA Fusion"),
            Self::TMB => write!(f, "TMB"),
            Self::HRD => write!(f, "HRD"),
            Self::MSI => write!(f, "MSI/MMR"),
            Self::Other(other) => write!(f, "Anderes Ergebnis ({other})"),
        }
    }
}

#[derive(Debug, Deserialize, PartialEq)]
pub struct Variant {
    #[serde(rename = "DETECTED_VAR_ID")]
    id: u32,
    #[serde(rename = "GENE_SYMBOL")]
    pub gene_symbol: Option<String>,
    #[serde(rename = "PROTEIN_MODIFICATION")]
    pub protein_modification: Option<String>,
    #[serde(rename = "PROTEIN_VARIANT_TYPE")]
    pub protein_variant_type: Option<ResultType>,
    #[serde(rename = "DISPLAY_VARIANT_TYPE")]
    pub display_variant_type: Option<ResultType>,
    #[serde(rename = "CHROMOSOMAL_MODIFIED_OBJECT")]
    pub chromosome: Option<String>,
    #[serde(rename = "CHROMOSOMAL_MODIFICATION")]
    pub chromosome_modification: Option<String>,
    #[serde(rename = "TRANSCRIPT_HGVS_MODIFIED_OBJECT")]
    pub transcript_hgvs_modified_object: Option<String>,
    #[serde(rename = "TOTAL_READS_IN_TUMOR")]
    pub total_reads_in_tumor: Option<u32>,
    #[serde(rename = "VARIANT_ALLELE_FREQUENCY_IN_TUMOR")]
    pub variant_allele_frequency_in_tumor: Option<f32>,
    #[serde(rename = "DBSNP")]
    pub db_snp: Option<String>,
    #[serde(rename = "COPY_NUMBER")]
    pub copy_number: Option<f32>,
    #[serde(rename = "CLASSIFICATION_NAME")]
    pub classification_name: Option<String>,
    #[serde(rename = "ONCOGENIC_CLASSIFICATION_NAME")]
    oncogenic_classification_name: Option<String>,
}

impl Variant {
    pub fn dna_change(&self) -> DnaChange {
        DnaChange::from_str(
            self.chromosome_modification
                .as_ref()
                .unwrap_or(&String::new()),
        )
        .unwrap_or_default()
    }
}

#[derive(Debug, Deserialize, PartialEq)]
pub struct Biomarkers {
    #[serde(rename = "NOTABLE_BIOMARKERS")]
    notable_biomarkers: Vec<NotableBiomarker>,
}

#[derive(Debug, Deserialize, PartialEq)]
pub struct NotableBiomarker {
    #[serde(rename = "BIOMARKERS")]
    biomarkers: Vec<Biomarker>,
}

#[derive(Debug, Deserialize, PartialEq)]
pub struct Biomarker {
    #[serde(rename = "DETECTED_VAR_ID")]
    id: u32,
    #[serde(rename = "DISPLAY_MODIFIED_OBJECT")]
    pub display_modified_object: Option<String>,
    #[serde(rename = "DISPLAY_VARIANT_TYPE")]
    pub display_variant_type: Option<ResultType>,
    #[serde(rename = "VARIANT_EFFECT")]
    pub variant_effect: Option<VariantEffect>,
    #[serde(rename = "TMB_VARIANT_COUNT_PER_MEGABASE")]
    pub tmb_variant_count_per_megabase: Option<String>,
    #[serde(rename = "COPY_NUMBER")]
    pub copy_number: Option<String>,
    #[serde(rename = "SCORE")]
    pub score: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mhguide::ResultType::{CopyNumberVariant, SimpleVariant, HRD, MSI, TMB};
    use crate::mhguide::{Biomarker, Biomarkers, General, MhGuide, NotableBiomarker, Variant};
    use crate::three_letter_protein_modification;
    use rstest::rstest;

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_sv_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    id: 12345678,
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: Some("p.A123V".to_string()),
                    protein_variant_type: Some(SimpleVariant("SNV".to_string())),
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345678G>A".to_string()),
                    transcript_hgvs_modified_object: Some("c.123C>T".to_string()),
                    total_reads_in_tumor: Some(567),
                    variant_allele_frequency_in_tumor: Some(42.42),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None,
                    classification_name: Some("Likely benign".to_string()),
                    oncogenic_classification_name: None
                }],
                biomarkers: Biomarkers {
                    notable_biomarkers: vec![NotableBiomarker {
                        biomarkers: vec![Biomarker {
                            id: 12345678,
                            display_modified_object: Some("TMB-L".to_string()),
                            display_variant_type: Some(TMB),
                            variant_effect: None,
                            tmb_variant_count_per_megabase: Some("0.1900".to_string()),
                            copy_number: None,
                            score: None,
                        }]
                    }]
                },
                report_narrative: String::new()
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_sv_del_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv_del-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    id: 12345678,
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: None,
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345670_12345678del".to_string()),
                    transcript_hgvs_modified_object: Some("c.120-1_128_1del".to_string()),
                    total_reads_in_tumor: Some(567),
                    variant_allele_frequency_in_tumor: Some(42.42),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None,
                    classification_name: Some("Likely benign".to_string()),
                    oncogenic_classification_name: None
                }],
                biomarkers: Biomarkers {
                    notable_biomarkers: vec![NotableBiomarker {
                        biomarkers: vec![Biomarker {
                            id: 12345678,
                            display_modified_object: Some("TMB-L".to_string()),
                            display_variant_type: Some(TMB),
                            variant_effect: None,
                            tmb_variant_count_per_megabase: Some("0.1900".to_string()),
                            copy_number: None,
                            score: None,
                        }]
                    }]
                },
                report_narrative: String::new()
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_cnv_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/cnv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![Variant {
                    id: 12345678,
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: None,
                    protein_variant_type: None,
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: None,
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: None
                }],
                biomarkers: Biomarkers {
                    notable_biomarkers: vec![NotableBiomarker {
                        biomarkers: vec![Biomarker {
                            id: 12345678,
                            display_modified_object: Some("EGFR".to_string()),
                            display_variant_type: Some(CopyNumberVariant),
                            variant_effect: Some(VariantEffect::CopyGain),
                            tmb_variant_count_per_megabase: None,
                            copy_number: Some("12.34".to_string()),
                            score: None,
                        }]
                    }]
                },
                report_narrative: String::new()
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_ref_allele() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(mhguide.variants.len(), 1);
        assert_eq!(
            mhguide.variants.first().unwrap().dna_change().ref_allele,
            "G"
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_alt_allele() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/sv-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(mhguide.variants.len(), 1);
        assert_eq!(
            mhguide.variants.first().unwrap().dna_change().alt_allele,
            "A"
        );
    }

    #[rstest]
    #[case("c.123C>T",
        DnaChange{ start: "123".to_string(), end: String::new(), ref_allele: "C".to_string(), alt_allele: "T".to_string() }
    )]
    #[case("c.-123C>T",
        DnaChange{ start: "-123".to_string(), end: String::new(), ref_allele: "C".to_string(), alt_allele: "T".to_string() }
    )]
    #[case("c.123_124insA",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: "A".to_string() }
    )]
    #[case("c.123_124del",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    #[case("c.-123_123del",
        DnaChange{ start: "-123".to_string(), end: "123".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    #[case("c.123_124delinsCTGA",
        DnaChange{ start: "123".to_string(), end: "124".to_string(), ref_allele: String::new(), alt_allele: "CTGA".to_string() }
    )]
    #[case("g.41149933A>G",
        DnaChange{ start: "41149933".to_string(), end: String::new(), ref_allele: "A".to_string(), alt_allele: "G".to_string() }
    )]
    #[case("g.41149933_41150000dup",
        DnaChange{ start: "41149933".to_string(), end: "41150000".to_string(), ref_allele: String::new(), alt_allele: String::new() }
    )]
    fn test_dna_change_parsing(#[case] case: &str, #[case] expected: DnaChange) {
        let actual = DnaChange::from_str(case);
        assert_eq!(actual, Ok(expected));
    }

    #[rstest]
    #[case("p.F123G", "p.Phe123Gly")]
    #[case("p.L123F", "p.Leu123Phe")]
    #[case("p.S123L", "p.Ser123Leu")]
    #[case("p.Y123S", "p.Tyr123Ser")]
    #[case("p.C123Y", "p.Cys123Tyr")]
    #[case("p.W123C", "p.Trp123Cys")]
    #[case("p.P123W", "p.Pro123Trp")]
    #[case("p.H123P", "p.His123Pro")]
    #[case("p.Q123H", "p.Gln123His")]
    #[case("p.R123Q", "p.Arg123Gln")]
    #[case("p.I123R", "p.Ile123Arg")]
    #[case("p.M123I", "p.Met123Ile")]
    #[case("p.T123M", "p.Thr123Met")]
    #[case("p.N123T", "p.Asn123Thr")]
    #[case("p.K123N", "p.Lys123Asn")]
    #[case("p.V123K", "p.Val123Lys")]
    #[case("p.A123V", "p.Ala123Val")]
    #[case("p.D123A", "p.Asp123Ala")]
    #[case("p.E123D", "p.Glu123Asp")]
    #[case("p.G123E", "p.Gly123Glu")]
    #[case("p.Y123=", "p.Tyr123=")]
    #[case("p.Y123fs", "p.Tyr123fs")]
    #[case("p.S123_I125delinsF", "p.Ser123_Ile125delinsPhe")]
    #[case("p.S123_I125delinsFE", "p.Ser123_Ile125delinsPheGlu")]
    #[case("p.S123_I125del", "p.Ser123_Ile125del")]
    #[case("p.Y123dup", "p.Tyr123dup")]
    // Examples from Onkostar Notices
    #[case("p.L858R", "p.Leu858Arg")]
    #[case("p.*del*", "p.*del*")]
    #[case("p.V600*", "p.Val600*")]
    // Not mappable - keep as is
    #[case("p.X123X", "p.X123X")]
    // Keep existing three letter codes
    #[case("p.Gly123Glu", "p.Gly123Glu")]
    #[case("p.Ser123_Ile125delinsPhe", "p.Ser123_Ile125delinsPhe")]
    // Examples from UKR
    #[case("p.E123Rfs*14", "p.Glu123Argfs*14")]
    #[case("p.E123Rfs*?", "p.Glu123Argfs*?")]
    fn test_three_letter_protein_modification(#[case] short: &str, #[case] long: &str) {
        assert_eq!(three_letter_protein_modification(short), long);
    }

    #[rstest]
    #[case("", 1)]
    #[case("KMT2C p.K1234fs laut XYZ oncogenic", 2)]
    #[case("KMT2C p.K1234fs laut XYZ oncogenic; FANCA p.S1234F noch dazu", 3)]
    #[case(
        "KMT2C p.K1234fs laut XYZ oncogenic; FANCA p.S1234F noch dazu; BRAF p.K1234F soll nicht doppelt sein",
        3
    )]
    #[case(
        "KMT2C p.K1234fs, ENST00000123456.12 c.1234dup, ins Nonsense, 25.23, Likely oncogenic",
        2
    )]
    fn test_should_add_protein_modification_report_narrative_matches(
        #[case] report_narrative: &str,
        #[case] expected_variants: usize,
    ) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345678,
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("KMT2C".to_string()),
                    protein_modification: Some("p.K1234fs".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("FANCA".to_string()),
                    protein_modification: Some("p.S1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: report_narrative.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), expected_variants);
    }

    #[rstest]
    #[case("", 1)]
    #[case("KMT2C c.123_124del laut XYZ oncogenic", 2)]
    #[case("KMT2C c.123_124del laut XYZ oncogenic; FANCA c.123A>G noch dazu", 3)]
    #[case(
        "KMT2C c.123_124del laut XYZ oncogenic; FANCA c.123A>G noch dazu; BRAF c.123T>C soll nicht doppelt sein",
        3
    )]
    fn test_should_add_cdna_modification_report_narrative_matches(
        #[case] report_narrative: &str,
        #[case] expected_variants: usize,
    ) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345678,
                    gene_symbol: Some("BRAF".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("KMT2C".to_string()),
                    protein_modification: Some("p.K1234fs".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123_124del".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("FANCA".to_string()),
                    protein_modification: Some("p.S1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123A>G".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: report_narrative.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), expected_variants);
    }

    #[rstest]
    #[case("A1BG-AS1 p.K1234F should be allowed, too", 1)]
    #[case("A1BG-AS1 c.123T>C should be allowed, too", 1)]
    #[case("APOBEC3A_B p.K1234F should be allowed, too", 1)]
    #[case("APOBEC3A_B c.123T>C should be allowed, too", 1)]
    fn test_allow_hyphen_underscore_in_symbol(
        #[case] report_narrative: &str,
        #[case] expected_variants: usize,
    ) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345678,
                    gene_symbol: Some("A1BG-AS1".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr19".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("APOBEC3A_B".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr22".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("benign".to_string()),
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: report_narrative.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), expected_variants);
    }

    #[rstest]
    #[case(
        "A1BG-AS1 p.K1234F liegt auf einem Homopolymer; mögliches Artefakt",
        false,
        2
    )]
    #[case(
        "A1BG-AS1 c.123T>C liegt auf einem Homopolymer; mögliches Artefakt",
        false,
        2
    )]
    #[case(
        "A1BG-AS1 p.K1234F liegt auf einem Homopolymer; wahrscheinlich ein Artefakt",
        false,
        2
    )]
    #[case(
        "A1BG-AS1 c.123T>C liegt auf einem Homopolymer; wahrscheinlich ein Artefakt",
        false,
        2
    )]
    #[case(
        "A1BG-AS1 p.K1234F liegt auf einem Homopolymer; mögliches Artefakt",
        true,
        1
    )]
    #[case(
        "A1BG-AS1 c.123T>C liegt auf einem Homopolymer; mögliches Artefakt",
        true,
        1
    )]
    #[case(
        "A1BG-AS1 p.K1234F liegt auf einem Homopolymer; wahrscheinlich ein Artefakt",
        true,
        1
    )]
    #[case(
        "A1BG-AS1 c.123T>C liegt auf einem Homopolymer; wahrscheinlich ein Artefakt",
        true,
        1
    )]
    fn test_remove_artifacts_if_required(
        #[case] report_narrative: &str,
        #[case] required: bool,
        #[case] expected_variants: usize,
    ) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345678,
                    gene_symbol: Some("A1BG-AS1".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr19".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("APOBEC3A_B".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr22".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: report_narrative.to_string(),
        };

        let actual = mh_guide.relevant_variants(required);

        assert_eq!(actual.len(), expected_variants);
    }

    #[rstest]
    #[case("KMT2C Copy number LOSS GCN = 0.00", 2)]
    #[case("KMT2C, GCN = 0.00", 2)]
    #[case("KMT2C, chr9, GCN: 0.00", 2)]
    fn test_add_cnv_report_narrative(
        #[case] report_narrative: &str,
        #[case] expected_variants: usize,
    ) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345678,
                    gene_symbol: Some("A1BG-AS1".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr19".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("KMT2C".to_string()),
                    protein_modification: Some("p.K1234fs".to_string()),
                    protein_variant_type: Some(CopyNumberVariant),
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(0.00),
                    classification_name: None,
                    oncogenic_classification_name: Some("Unclassified".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("KMT2C".to_string()),
                    protein_modification: Some("p.K1234fs".to_string()),
                    protein_variant_type: Some(CopyNumberVariant),
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: None,
                    classification_name: None,
                    oncogenic_classification_name: Some("Unclassified".to_string()),
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: report_narrative.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), expected_variants);
    }

    #[rstest]
    #[case(Biomarkers {
                notable_biomarkers: vec![],
        }, 1)]
    #[case(Biomarkers {
                notable_biomarkers: vec![
                    NotableBiomarker {
                        biomarkers: vec![Biomarker {
                            id: 12345678,
                            display_modified_object: None,
                            display_variant_type: Some(CopyNumberVariant),
                            variant_effect: Some(VariantEffect::CopyGain),
                            tmb_variant_count_per_megabase: None,
                            copy_number: Some("12.34".to_string()),
                            score: None,
                        }]
                    },
                ],
        }, 2)]
    fn test_add_cnv_biomarkers(#[case] biomarkers: Biomarkers, #[case] expected_variants: usize) {
        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![
                Variant {
                    id: 12345600,
                    gene_symbol: Some("A1BG-AS1".to_string()),
                    protein_modification: Some("p.K1234F".to_string()),
                    protein_variant_type: None,
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr19".to_string()),
                    chromosome_modification: None,
                    transcript_hgvs_modified_object: Some("c.123T>C".to_string()),
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(12.34),
                    classification_name: None,
                    oncogenic_classification_name: Some("oncogenic".to_string()),
                },
                Variant {
                    id: 12345678,
                    gene_symbol: Some("EGFR".to_string()),
                    protein_modification: Some("Copy number Gain".to_string()),
                    protein_variant_type: Some(CopyNumberVariant),
                    display_variant_type: Some(CopyNumberVariant),
                    chromosome: Some("chr7".to_string()),
                    chromosome_modification: Some("Chr7:12345_54321gain".to_string()),
                    transcript_hgvs_modified_object: None,
                    total_reads_in_tumor: Some(123),
                    variant_allele_frequency_in_tumor: None,
                    db_snp: None,
                    copy_number: Some(87.),
                    classification_name: Some("Unclassified".to_string()),
                    oncogenic_classification_name: Some("Unclassified".to_string()),
                },
            ],
            biomarkers,
            report_narrative: String::new(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), expected_variants);
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_biomarker_deserialization() {
        static SV_MHGUIDE: &str = include_str!("../testfiles/biomarkers-mhguide.json");

        let mhguide = serde_json::from_str::<MhGuide>(SV_MHGUIDE).unwrap();
        assert_eq!(
            mhguide,
            MhGuide {
                general: General {
                    order_date: "2026-02-11".to_string(),
                    ref_genome_version: RefGenomeVersion::Hg19,
                    patient_identifier: PatientIdentifier {
                        h_number: "H10000-26".to_string(),
                        pid: "PID0123456".to_string()
                    }
                },
                variants: vec![],
                biomarkers: Biomarkers {
                    notable_biomarkers: vec![
                        NotableBiomarker {
                            biomarkers: vec![Biomarker {
                                id: 12345670,
                                display_modified_object: Some("TMB-L".to_string()),
                                display_variant_type: Some(TMB),
                                variant_effect: None,
                                tmb_variant_count_per_megabase: Some("0.1900".to_string()),
                                copy_number: None,
                                score: None,
                            }]
                        },
                        NotableBiomarker {
                            biomarkers: vec![Biomarker {
                                id: 12345671,
                                display_modified_object: Some("HRD-positive".to_string()),
                                display_variant_type: Some(HRD),
                                variant_effect: None,
                                tmb_variant_count_per_megabase: None,
                                copy_number: None,
                                score: Some("12.00".to_string()),
                            }]
                        },
                        NotableBiomarker {
                            biomarkers: vec![Biomarker {
                                id: 12345672,
                                display_modified_object: Some("MS-stable".to_string()),
                                display_variant_type: Some(MSI),
                                variant_effect: None,
                                tmb_variant_count_per_megabase: None,
                                copy_number: None,
                                score: Some("0.21".to_string()),
                            }]
                        }
                    ]
                },
                report_narrative: String::new()
            }
        );
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_extract_hrd_from_biomarkers() {
        static MHGUIDE: &str = include_str!("../testfiles/biomarkers-mhguide.json");

        let value = serde_json::from_str::<MhGuide>(MHGUIDE)
            .unwrap()
            .hrd_score();
        assert_eq!(value, Some(12.0));
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_extract_msi_from_biomarkers() {
        static MHGUIDE: &str = include_str!("../testfiles/biomarkers-mhguide.json");

        let value = serde_json::from_str::<MhGuide>(MHGUIDE)
            .unwrap()
            .msi_score();
        assert_eq!(value, Some(0.21));
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    fn test_extract_tmb_from_biomarkers() {
        static MHGUIDE: &str = include_str!("../testfiles/biomarkers-mhguide.json");

        let value = serde_json::from_str::<MhGuide>(MHGUIDE)
            .unwrap()
            .tmb_value();
        assert_eq!(value, Some(0.19));
    }

    #[test]
    #[allow(clippy::unwrap_used)]
    #[allow(clippy::panic)]
    fn test_sv_examples_from_report_specification() {
        static INPUT: &str = "Variant, Coding DNA, Type and effect, VAF (%), Classification\n\
            KIAA0586 p.R116fs, ENST00000619416.4 c.347del, del Frameshift, 47.25, Likely oncogenic\n\
            KMT2C p.Y816*, ENST00000262189.11 c.2447dup, ins Nonsense, 25.23, Likely oncogenic\n\
            LMF1 p.V63fs, ENST00000262301.16 c.187del, del Frameshift, 46.84, Likely oncogenic";

        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },

            variants: vec![
                Variant {
                    id: 10000000,
                    gene_symbol: Some("KIAA0586".to_string()),
                    protein_modification: Some("p.R116fs".to_string()),
                    protein_variant_type: Some(SimpleVariant("SNV".to_string())),
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345678del".to_string()),
                    transcript_hgvs_modified_object: Some("c.347del".to_string()),
                    total_reads_in_tumor: Some(567),
                    variant_allele_frequency_in_tumor: Some(47.25),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None,
                    classification_name: Some("Likely benign".to_string()),
                    oncogenic_classification_name: None,
                },
                Variant {
                    id: 20000000,
                    gene_symbol: Some("KMT2C".to_string()),
                    protein_modification: Some("p.Y816*".to_string()),
                    protein_variant_type: Some(SimpleVariant("SNV".to_string())),
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345678dup".to_string()),
                    transcript_hgvs_modified_object: Some("c.2447dup".to_string()),
                    total_reads_in_tumor: Some(567),
                    variant_allele_frequency_in_tumor: Some(47.25),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None,
                    classification_name: Some("Likely benign".to_string()),
                    oncogenic_classification_name: None,
                },
                Variant {
                    id: 30000000,
                    gene_symbol: Some("LMF1".to_string()),
                    protein_modification: Some("p.V63fs".to_string()),
                    protein_variant_type: Some(SimpleVariant("SNV".to_string())),
                    display_variant_type: Some(SimpleVariant("SNV".to_string())),
                    chromosome: Some("chr1".to_string()),
                    chromosome_modification: Some("g.12345678del".to_string()),
                    transcript_hgvs_modified_object: Some("c.187del".to_string()),
                    total_reads_in_tumor: Some(567),
                    variant_allele_frequency_in_tumor: Some(47.25),
                    db_snp: Some("rs202602111".to_string()),
                    copy_number: None,
                    classification_name: Some("Likely benign".to_string()),
                    oncogenic_classification_name: None,
                },
            ],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: INPUT.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), 3);
    }

    #[test]
    fn test_cnv_examples_from_report_specification() {
        static INPUT: &str = "Gene, Chromosome number, Gene Copy Number (GCN)\n\
            ELAVL2, chr9, GCN: 12.2";

        let mh_guide = MhGuide {
            general: General {
                order_date: "2026-02-11".to_string(),
                ref_genome_version: RefGenomeVersion::Hg19,
                patient_identifier: PatientIdentifier {
                    h_number: "H10000-26".to_string(),
                    pid: "PID0123456".to_string(),
                },
            },
            variants: vec![Variant {
                id: 12345678,
                gene_symbol: Some("ELAVL2".to_string()),
                protein_modification: Some("p.K1234fs".to_string()),
                protein_variant_type: Some(CopyNumberVariant),
                display_variant_type: Some(CopyNumberVariant),
                chromosome: Some("chr9".to_string()),
                chromosome_modification: None,
                transcript_hgvs_modified_object: None,
                total_reads_in_tumor: Some(123),
                variant_allele_frequency_in_tumor: None,
                db_snp: None,
                copy_number: Some(12.20),
                classification_name: None,
                oncogenic_classification_name: Some("Unclassified".to_string()),
            }],
            biomarkers: Biomarkers {
                notable_biomarkers: vec![],
            },
            report_narrative: INPUT.to_string(),
        };

        let actual = mh_guide.relevant_variants(false);

        assert_eq!(actual.len(), 1);
    }
}
