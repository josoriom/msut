use crate::utilities::parse::decode::decode;
use crate::utilities::parse::parse_mzml::MzML;

pub fn bin_to_json(bin: &[u8], pretty: bool) -> Result<String, String> {
    let mzml: MzML = decode(bin)?;
    if pretty {
        serde_json::to_string_pretty(&mzml).map_err(|e| e.to_string())
    } else {
        serde_json::to_string(&mzml).map_err(|e| e.to_string())
    }
}
