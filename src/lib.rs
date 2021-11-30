use core::str;
use std::fmt;
use std::{collections::HashMap, fs};
use std::io;
use structopt::StructOpt;
//use tfbs::{self, Sequence};

/// A DNA sequence have bases and indices. These indices are used to save base position and
/// they are necessary when the sequences have gaps (ie after alignmnent).
pub struct Sequence {
    idxs: Vec<usize>,
    bases: Vec<char>,
}

impl From<&str> for Sequence {
    /// Create a Sequence from a `seq`. Gaps '-' are removed from bases and indices.
    /// There is no test for valid bases (ACGT or acgt).
    fn from(seq: &str) -> Self {
        seq.char_indices()
            .filter(|&(_,x)| !x.eq(&'-'))
            .fold(Sequence { bases: Vec::<char>::new(), idxs: Vec::<usize>::new()}, | mut s, b| {
                s.idxs.push(b.0);
                s.bases.push(b.1);
                s
            })
    }
}

/// DNA strand forward or reverse.
#[derive(Debug)]
pub enum Strand {
    Forward, 
    Reverse,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
        }
    }
}


/// Score are assigned to a range, so there is a `start` and an `end` besides the `score` value.
/// Reverse strand scores have `start` > `end`. 
/// The score length is the same as matrix length
#[derive(Debug)]
pub struct Score {
    seq_start: usize,
    seq_end: usize,
    algn_start: usize,
    algn_end: usize,
    score: f64,
}

impl fmt::Display for Score {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // write!(f, "{}\t{}\t{:.3}", self.start, self.end, self.score )
        write!(f, "{}\t{}\t{}\t{}\t{:.3}", self.seq_start, self.seq_end, self.algn_start, self.algn_end, self.score )
    }
}

/// DNA matrix 
/// # Atributes
///
#[derive(Debug)]
pub struct DNAMatrix {
    pub name: String,
    length: usize,
    probs: Vec<Vec<f64>>,
    conservation: Vec<f64>,
    max_score: f64,
    threshold: f64,
    pub strand: Strand
}

impl DNAMatrix {
    pub fn new(name: &str, threshold: f64, counts: &Vec<Vec<f64>>, strand: Strand) -> Self {
        let probs = match strand {
            Strand::Forward => Self::calculate_probs(counts),
            Strand::Reverse => Self::calculate_probs(&Self::comp_rev_counts(counts)),
        };

        let mut matrix = DNAMatrix {
            name: name.to_string(),
            length: counts.len(),
            probs,
            conservation: vec![],
            max_score: 0.0,
            threshold,
            strand,
        };
        matrix.calculate_conservation();
        matrix.calculate_max_score();
        matrix
    }

    fn comp_rev_counts(v: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let mut r = v.clone();
        r.reverse();
        for b in r.iter_mut() {
            b.reverse();
        }
        r
    }

    fn calculate_probs(counts: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let mut probs = vec![];
        let max_count: f64 = counts
            .iter()
            .map(|x| x.iter().sum())
            .fold(f64::NEG_INFINITY, f64::max);
        for position in counts {
            if position.len() != 4 {
                panic!("Matrix has {} values when 4 is expected.", position.len())
            }
            let mut v: Vec<f64> = position.iter().map(|x| x / max_count).collect();
            let sum: f64 = position.iter().sum();
            v.push(1.0 - sum / max_count);
            probs.push(v);
        }
        probs
    }

    fn calculate_conservation(&mut self) {
        for position in self.probs.iter() {
            let sum: f64 = position
                .iter()
                .fold(0.0, |acc, x| if x > &0.0 { acc + x * x.ln() } else { acc });
            let c = sum / f64::ln(position.len() as f64) + 1.0;
            self.conservation.push(c);
        }
    }

    fn calculate_max_score(&mut self) {
        self.max_score = self
            .probs
            .iter()
            .map(|p| {
                p.iter()
                    .take(4)
                    .fold(f64::NEG_INFINITY, |state, x| state.max(*x))
            })
            .enumerate()
            .fold(0.0, |acc, (i, v)| acc + self.conservation[i] * v);
    }

    pub fn scan(&self, seq: &Sequence) -> Vec<Score> {
        let (s,e) = match self.strand {
            Strand::Forward => (1, self.length),
            Strand::Reverse => (self.length, 1),
        };
        let scores: Vec<Score> = seq.bases
            .windows(self.length)
            .map(
                |w| {
                    w.iter().enumerate().fold(0.0, |acc, (i, b)| {
                        acc + self.probs[i][Self::lookup(b)] * self.conservation[i]
                    })
                }, 
            )
            .enumerate()
            .map(|v: (usize, f64)| Score {seq_start: v.0 + s, seq_end: v.0 + e, algn_start: seq.idxs[v.0 + s - 1] + 1, algn_end: seq.idxs[v.0 + e - 1] + 1, score: v.1 / self.max_score})
            .filter(|v| v.score >= self.threshold)
            .collect();
        scores
    }

    fn lookup(b: &char) -> usize {
        match b {
            'A' | 'a' => 0,
            'C' | 'c' => 1,
            'G' | 'g' => 2,
            'T' | 't' => 3,
            _ => panic!("DNA base unknown {}", b),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn test_calculate_probs() {
    //     let v = super::DNAMatrix::calculate_probs(&[
    //         &[2.0, 0.0, 0.0, 0.0],
    //         &[1.0, 0.0, 0.0, 1.0],
    //         &[1.0, 0.0, 0.0, 0.0],
    //         &[0.50, 0.50, 0.50, 0.50],
    //     ]);
    //     println!("{:?}", v);
    // }

    // #[test]
    // fn test_create_dna_matrix() {
    //     let c: &[&[f64]] = &[
    //         &[2.0, 0.0, 0.0, 0.0],
    //         &[1.0, 0.0, 0.0, 1.0],
    //         &[1.0, 0.0, 0.0, 0.0],
    //         &[0.50, 0.50, 0.50, 0.50],
    //     ];
    //     let m = super::DNAMatrix::new("teste", 0.1, c);
    //     println!("{:?}", m.conservation);
    //     println!("{:?}", m.max_score);
    // }

    // #[test]
    // fn test_scan() {
    //     let c: &[&[f64]] = &[
    //         &[2.0, 0.0, 0.0, 0.0],
    //         &[1.0, 0.0, 0.0, 1.0],
    //         &[1.0, 0.0, 0.0, 0.0],
    //         &[0.50, 0.50, 0.50, 0.50],
    //     ];
    //     let m = super::DNAMatrix::new("teste", 0.5, c, Strand::Forward);
    //     let seq = b"ACGTACGTACGTAGATGTCTAGTACGTACGCTAGCTAGCTGAGACTGACTAGTACGTAAGCTAGCACG";
    //     let seqb = b"AAAAACCCCCGGGGTTTTTCGTACGTACGTAGATGTCTAGTACGTACGCTAGCTAGCTGAGACTGACTAGTACGTAAGCTAGCACG";
    //     let scores = m.scan(seq);
    //     let scoresb = m.scan(seqb);
    //     println!("{:?}", scores);
    //     println!("{:?}", scoresb);
    // }

    #[test]
    fn test_scan_forward() {
        let c = &vec![
            vec![2.0, 0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.50, 0.50, 0.50, 0.50],
        ];
        let s = Strand::Forward;
        let m = super::DNAMatrix::new("teste", 0.5, c, s);
        let seq = Sequence::from("-ACG-TACGTACGTAGATGTCTAGTACGTACGCTAGCTAGCTGAGACTGACTAGTACGTAAGCTAGCACG");
        let scores = m.scan(&seq);
        println!("{:?}", scores);
    }

    #[test]
    fn test_scan_reverse() {
        let c = &vec![
            vec![2.0, 0.0, 0.0, 0.0],
            vec![1.0, 0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.50, 0.50, 0.50, 0.50],
        ];
        let s = Strand::Reverse;
        let m = super::DNAMatrix::new("teste", 0.5, c, s);
        let seq = Sequence::from("-ACG-TACGTACGTAGATGTCTAGTACGTACGCTAGCTAGCTGAGACTGACTAGTACGTAAGCTAGCACG");
        let scores = m.scan(&seq);
        println!("{:?}", scores);
    }

    // #[test]
    // fn create_dna_matrix() {
    //     let v = vec![vec![837.0, 1889.0, 1280.0, 718.0], vec![193.0, 0.0, 0.0, 4725.0], vec![4725.0, 65.0, 275.0, 1232.0]];
    //     let m = super::DNAMatrix {
    //         name: String::from("teste"),
    //         length: v.len(),
    //         probs: v,
    //         alphabet: bio::alphabets::dna::alphabet(),
    //     };
    //     println!("Matriz {} has length {} and probabilities {:?}", m.name, m.length, m.probs);
    // }

    // #[test]
    // fn it_works() {
    //     let x = 0.0 * 0.00_f64.ln();
    //     println!("{}", x);
    //     assert_eq!(2 + 2, 4);
    // }

    // #[test]
    // fn dnamotif() {
    //     use bio::pattern_matching::pssm::{DNAMotif, Motif};
    //     let pssm = DNAMotif::from_seqs(
    //         vec![
    //             b"AAAAA".to_vec(),
    //             b"CAATA".to_vec(),
    //             b"TAAGA".to_vec(),
    //             b"CAAAA".to_vec(),
    //         ]
    //         .as_ref(),
    //         None,
    //     )
    //     .unwrap();
    //     println!("matrix scores {:?}", pssm.scores);
    //     let start_pos = pssm.score(b"CCCCCAATA").unwrap().scores;
    //     println!("motif found at position {:?}", start_pos);
    // }

    // #[test]
    // fn dnamotif_from_array() {
    //     use bio::pattern_matching::pssm::DNAMotif;
    //     use ndarray::Array2;
    //     // BARHL
    //     let p: Array2<f64> = ndarray::array![
    //         [837.0, 1889.0, 1280.0, 718.0],
    //         [193.0, 0.0, 0.0, 4725.0],
    //         [4725.0, 65.0, 275.0, 1232.0],
    //         [4725.0, 55.0, 0.0, 0.0],
    //         [4725.0, 199.0, 662.0, 2114.0],
    //         [0.0, 4725.0, 0.0, 2008.0],
    //         [622.0, 0.0, 4725.0, 259.0],
    //         [1065.0, 1001.0, 1508.0, 1151.0],
    //     ];
    //     let pssm = DNAMotif::from(p);
    //     println!("matrix scores {:?}", pssm.scores);
    // }

    // #[test]
    // fn dnamotif_from_array_with_gap() {
    //     // normalization problem
    //     use bio::pattern_matching::pssm::DNAMotif;
    //     use ndarray::Array2;
    //     // BARHL
    //     let p: Array2<f64> = ndarray::array![
    //         [837.0, 1889.0, 1280.0, 718.0, 1.0],
    //         [193.0, 0.0, 0.0, 4725.0, 1.0],
    //         [4725.0, 65.0, 275.0, 1232.0, 1.0],
    //         [4725.0, 55.0, 0.0, 0.0, 1.0],
    //         [4725.0, 199.0, 662.0, 2114.0, 1.0],
    //         [0.0, 4725.0, 0.0, 2008.0, 1.0],
    //         [622.0, 0.0, 4725.0, 259.0, 1.0],
    //         [1065.0, 1001.0, 1508.0, 1151.0, 1.0],
    //     ];
    //     let pssm = DNAMotif::from(p);
    //     println!("matrix scores {:?}", pssm.scores);
    // }

    // #[test]
    // fn pssm_info() {
    //     use bio::pattern_matching::pssm::{DNAMotif, Motif};
    //     use ndarray::Array2;
    //     let p: Array2<f64> = ndarray::array![
    //         [837.0, 1889.0, 1280.0, 718.0, 1.0],
    //         [193.0, 0.0, 0.0, 4725.0, 1.0],
    //         [4725.0, 65.0, 275.0, 1232.0, 1.0],
    //         [4725.0, 55.0, 0.0, 0.0, 1.0],
    //         [4725.0, 199.0, 662.0, 2114.0, 1.0],
    //         [0.0, 4725.0, 0.0, 2008.0, 1.0],
    //         [622.0, 0.0, 4725.0, 259.0, 1.0],
    //         [1065.0, 1001.0, 1508.0, 1151.0, 1.0],
    //     ];
    //     let pssm = DNAMotif::from(p);
    //     let info = pssm.info_content();
    //     println!("info content {:?}", info);
    // }

    // #[test]
    // fn pssm_search() {
    //     use bio::pattern_matching::pssm::{DNAMotif, Motif};
    //     use ndarray::Array2;
    //     let p: Array2<f64> = ndarray::array![
    //         [837.0, 1889.0, 1280.0, 718.0, 1.0],
    //         [193.0, 0.0, 0.0, 4725.0, 1.0],
    //         [4725.0, 65.0, 275.0, 1232.0, 1.0],
    //         [4725.0, 55.0, 0.0, 0.0, 1.0],
    //         [4725.0, 199.0, 662.0, 2114.0, 1.0],
    //         [0.0, 4725.0, 0.0, 2008.0, 1.0],
    //         [622.0, 0.0, 4725.0, 259.0, 1.0],
    //         [1065.0, 1001.0, 1508.0, 1151.0, 1.0],
    //     ];
    //     let pssm = DNAMotif::from(p);
    //     let start_pos = pssm.score(b"CCCCCAATAAGTCAGTAGCTCGTAGCTACGTAGTAGCTACGAT").unwrap().sum;
    //     println!("motif found at position {:?}", start_pos);
    // }
}



pub fn tfbs_print_git() {
        println!("Hello, git!");
}

#[derive(Debug)]
struct JasparMatrix {
    id: String,
    counts: Vec<Vec<f64>>,
}

#[derive(Debug)]
pub enum Threshold {
    Fixed(f64),
    Specific(HashMap<String, f64>)
    // Specific(String)
}

pub fn read_specific_thresholds(filename: &str) -> HashMap<String, f64> {
    let data = fs::read_to_string(filename).expect("Unable to read file");
    let thresholds: HashMap<String, f64> = data.lines()
    .fold(HashMap::<String,f64>::new(), |mut h, line|{
        if let Some((name, cutoff)) = line.split_once("\t") {
            h.insert(name.to_string(), cutoff.parse::<f64>().expect("threshold should be a number"));
        };
        h
    });
    thresholds 
}

#[derive(Debug)]
pub fn read_matrix_from_jaspar_file(filename: &str, threshold: &Threshold) -> Vec<DNAMatrix> {
    let data = fs::read_to_string(filename).expect("Unable to read file");


    let matrices: Vec<DNAMatrix> = data
        .lines()
        .fold(Vec::<JasparMatrix>::new(), |mut m, line| -> Vec<JasparMatrix> {
            if line.starts_with("ID") {
                // Some(PFM::Name(line.trim_start_matches("ID").trim().to_string()))
                let name = line.trim_start_matches("ID").trim().to_string();
                m.push(JasparMatrix {id: name, counts: vec![]});
                m
            } else if line.starts_with(char::is_numeric) {
                m.last_mut().unwrap().counts.push(line.trim_start_matches(char::is_numeric)
                .trim()
                .split("\t")
                .map(|x| x.parse::<f64>().expect("count should be a number"))
                .collect());
                m
                // Some(PFM::Counts(v))
            } else {
                m
            }
        })
        .into_iter()
        .fold(Vec::<DNAMatrix>::new(), |mut acc, m| {
            let t = match threshold {
                Threshold::Fixed(t) => t,
                Threshold::Specific(s) => s.get(&m.id).unwrap(),
            };
            acc.push(DNAMatrix::new(&m.id, *t, &m.counts, Strand::Forward));
            acc.push(DNAMatrix::new(&m.id, *t, &m.counts, Strand::Reverse));
            acc
        });
    matrices
}

#[derive(StructOpt)]
struct Opt {
    #[structopt(name = "matrices", short, long)]
    jaspar_matrix_file: String, 
    #[structopt(name = "threshold", short, long)]
    threshold: String,
    #[structopt(name = "sequence", short, long)]
    sequence: String,
}


//let th:&str="0.8";
//let ma:&str="/home/helder/Projetos/INSECTM/ze_rust_scan/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt";
//let se:&str="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

//pub fn run<'a>(th: &'a str, ma: &'a str, se: &'a str)  {
//    //let opt = Opt::from_args(); 
//    let threshold = match th.parse::<f64>() {
//        Ok(number) => Threshold::Fixed(number),
//        Err(_) => Threshold::Specific(read_specific_thresholds(&th)),
//    };
//    let matrices = read_matrix_from_jaspar_file(&ma, &threshold);
//    let seq = Sequence::from(se);
//
//    for m in matrices {
//        let scores =  m.scan(&seq);
//        print!("{}", scores.iter().fold(String::new(), |acc, arg| acc + &format!("{}\t{}\t{}\n", m.name, arg, m.strand)));
//    }
//}

//run!("0.8", "/home/helder/Projetos/INSECTM/ze_rust_scan/JASPAR2020_CORE_vertebrates_non-redundant_pfms_transfac.txt", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
