use statrs::distribution::DiscreteCDF;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::collections::HashSet;
use tempfile::tempdir;

#[test]
fn runs_binary_and_reports_metrics() {
    let temp_dir = tempdir().expect("tempdir");
    let fasta_path = temp_dir.path().join("sample.fa");

    write_fasta(&fasta_path, ">seq1\nACGTACG\n");

    let assert = assert_cmd::cargo::cargo_bin_cmd!("motif_rater")
        .args([
            "--motif",
            "ACG",
            "--genome-fasta-file",
            fasta_path.to_str().unwrap(),
            "--header=false",
        ])
        .assert();

    let output = String::from_utf8(assert.get_output().stdout.clone()).expect("utf8 output");
    let fields: Vec<&str> = output.trim().split('\t').collect();

    assert_eq!(fields.len(), 8);
    assert_eq!(fields[0], "sample.fa");
    assert_eq!(fields[1], "7");
    assert_eq!(fields[2], "0.571429");
    assert_eq!(fields[3], "ACG");
    assert_eq!(fields[4], "3");
    assert_eq!(fields[5], "0.175");
    assert_eq!(fields[6], "0.600000");

    let p_value: f64 = fields[7].parse().expect("p-value is float");
    let poisson = statrs::distribution::Poisson::new(0.175).expect("lambda");
    let expected_p = 1.0 - poisson.cdf(2);
    assert!((p_value - expected_p).abs() < 1e-6);

    temp_dir.close().unwrap();
}

#[test]
fn motif_length_deduplicates_palindromes() {
    let temp_dir = tempdir().expect("tempdir");
    let fasta_path = temp_dir.path().join("sample.fa");

    write_fasta(&fasta_path, ">seq1\nACGTACGT\n");

    let assert = assert_cmd::cargo::cargo_bin_cmd!("motif_rater")
        .args([
            "--motif-length",
            "4",
            "--genome-fasta-file",
            fasta_path.to_str().unwrap(),
            "--header=false",
        ])
        .assert()
        .success();

    let output = String::from_utf8(assert.get_output().stdout.clone()).expect("utf8 output");
    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 136);

    let mut canonical = HashSet::new();
    for line in lines {
        let motif = line
            .split('\t')
            .nth(3)
            .expect("motif column present")
            .to_string();
        let rc = simple_reverse_complement(&motif);
        let canonical_form = motif.min(rc);
        assert!(canonical.insert(canonical_form));
    }

    temp_dir.close().unwrap();
}

fn simple_reverse_complement(motif: &str) -> String {
    motif
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}

fn write_fasta(path: &Path, contents: &str) {
    let mut file = File::create(path).expect("create fasta");
    file.write_all(contents.as_bytes()).expect("write fasta");
}
