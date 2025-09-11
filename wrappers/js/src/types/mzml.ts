export interface CvEntry {
  id: string | null;
  fullName: string | null;
  version: string | null;
  uri: string | null;
}
export interface CvPair {
  name: string | null;
  value: string | null;
  accession: string | null;
  cvRef: string | null;
  unitName: string | null;
}
export interface SourceFile {
  id: string | null;
  name: string | null;
  location: string | null;
  cvParams: CvPair[];
}
export interface FileDescription {
  fileContent: CvPair[];
  sourceFiles: SourceFile[];
}
export interface RefParamGroup {
  id: string | null;
  cvParams: CvPair[];
}
export interface Sample {
  id: string | null;
  name: string | null;
  cvParams: CvPair[];
}
export interface UserParam {
  name: string | null;
  value: string | null;
}
export interface Software {
  id: string | null;
  version: string | null;
  cvParams: CvPair[];
  userParams: UserParam[];
}
export interface Component {
  kind: string | null;
  order: number | null;
  cvParams: CvPair[];
}
export interface InstrumentConfiguration {
  id: string | null;
  refParamGroup: string | null;
  cvParams: CvPair[];
  components: Component[];
  softwareRef: string | null;
}
export interface ProcessingMethod {
  order: number | null;
  softwareRef: string | null;
  cvParams: CvPair[];
}
export interface DataProcessing {
  id: string | null;
  methods: ProcessingMethod[];
}
export interface AcquisitionSettings {
  id: string | null;
  instrumentRef: string | null;
  cvParams: CvPair[];
}
export interface SpectrumSummary {
  index: number;
  id: string | null;
  arrayLength: number;
  msLevel: number | null;
  scanType: string | null;
  polarity: string | null;
  spectrumType: string | null;
  retentionTime: number | null;
  scanWindowLowerLimit: number | null;
  scanWindowUpperLimit: number | null;
  totalIonCurrent: number | null;
  basePeakIntensity: number | null;
  basePeakMZ: number | null;
  precursor?: {
    isolationWindowTargetMz?: number | null;
    isolationWindowLowerOffset?: number | null;
    isolationWindowUpperOffset?: number | null;
    selectedIonMz?: number | null;
  } | null;
  mzArray: Float64Array | Float32Array | null;
  intensityArray: Float64Array | Float32Array | null;
}
export interface ChromatogramSummary {
  index: number;
  id: string | null;
  arrayLength: number;
  timeArray: Float64Array | Float32Array | null;
  intensityArray: Float64Array | Float32Array | null;
}
export interface Run {
  id: string | null;
  startTimeStamp: string | null;
  defaultInstrumentConfigurationRef: string | null;
  spectrumListCount: number | null;
  chromatogramListCount: number | null;
  spectra: SpectrumSummary[];
  chromatograms: ChromatogramSummary[];
}
export interface IndexOffset {
  idRef: string | null;
  offset: bigint;
}
export interface IndexList {
  spectrum: IndexOffset[];
  chromatogram: IndexOffset[];
  indexListOffset: bigint | null;
  fileChecksum: string | null;
}
export interface MzML {
  cvList: CvEntry[];
  fileDescription: FileDescription | null;
  referenceableParamGroups: RefParamGroup[];
  sampleList: Sample[];
  instrumentConfigurations: InstrumentConfiguration[];
  softwareList: Software[];
  dataProcessingList: DataProcessing[];
  acquisitionSettingsList: AcquisitionSettings[];
  run: Run | null;
  indexList: IndexList | null;
}
