use crate::utilities::parse::parse_mzml::{MzML, parse_mzml};
use std::process; // needed if you want to exit

pub fn mzml_to_json(bytes: &[u8]) -> Result<String, serde_json::Error> {
    let mzml: MzML = parse_mzml(bytes, true).unwrap_or_else(|e| {
        eprintln!("parse_mzml error: {e}");
        process::exit(1);
    });
    serde_json::to_string(&mzml) // pass a reference!
}
