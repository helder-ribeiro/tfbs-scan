use std::{collections::HashMap, fs};
use std::io;
use structopt::StructOpt;
use tfbs::{self, Sequence};


#[derive(Debug)]
struct JasparMatrix {
    id: String,
    counts: Vec<Vec<f64>>,
}

enum Threshold {
    Fixed(f64),
    Specific(HashMap<String, f64>)
    // Specific(String)
}

fn read_specific_thresholds(filename: &str) -> HashMap<String, f64> {
    let data = fs::read_to_string(filename).expect("Unable to read file");
    let thresholds: HashMap<String, f64> = data.lines()
    .fold(HashMap::<String,f64>::new(), |mut h, line|{
        if let Some((name, cutoff)) = line.split_once("\t") {
            h.insert(name.to_string(), cutoff.parse::<f64>().expect("threshould should be a number"));
        };
        h
    });
    thresholds 
}

fn read_matrix_from_jaspar_file(filename: &str, threshold: &Threshold) -> Vec<tfbs::DNAMatrix> {
    let data = fs::read_to_string(filename).expect("Unable to read file");


    let matrices: Vec<tfbs::DNAMatrix> = data
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
        .fold(Vec::<tfbs::DNAMatrix>::new(), |mut acc, m| {
            let t = match threshold {
                Threshold::Fixed(t) => t,
                Threshold::Specific(s) => s.get(&m.id).unwrap(),
            };
            acc.push(tfbs::DNAMatrix::new(&m.id, *t, &m.counts, tfbs::Strand::Forward));
            acc.push(tfbs::DNAMatrix::new(&m.id, *t, &m.counts, tfbs::Strand::Reverse));
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

fn main() {
    let opt = Opt::from_args(); 
    let threshold = match opt.threshold.parse::<f64>() {
        Ok(number) => Threshold::Fixed(number),
        Err(_) => Threshold::Specific(read_specific_thresholds(&opt.threshold)),
    };
    let matrices = read_matrix_from_jaspar_file(&opt.jaspar_matrix_file, &threshold);
    let seq = Sequence::from(opt.sequence.as_str());

    for m in matrices {
        let scores =  m.scan(&seq);
        print!("{}", scores.iter().fold(String::new(), |acc, arg| acc + &format!("{}\t{}\t{}\n", m.name, arg, m.strand)));
    }
}


