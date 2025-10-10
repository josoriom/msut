import * as fs from "fs";
import * as path from "path";

function firstExisting(...candidates: string[]) {
  for (const p of candidates) if (fs.existsSync(p)) return p;
  return candidates[0];
}

function platformLibPath(process: NodeJS.Process): string {
  const base = path.join(__dirname, "..", "native");
  const { platform, arch } = process;
  const file =
    platform === "win32"
      ? "msut.dll"
      : platform === "darwin"
      ? "libmsut.dylib"
      : "libmsut.so";
  let dir: string;
  if (platform === "darwin") {
    dir =
      arch === "arm64"
        ? firstExisting(
            path.join(base, "darwin-arm64"),
            path.join(base, "macos-arm64")
          )
        : firstExisting(
            path.join(base, "darwin-x64"),
            path.join(base, "macos-x86_64")
          );
  } else if (platform === "linux") {
    dir =
      arch === "arm64"
        ? firstExisting(
            path.join(base, "linux-arm64-gnu"),
            path.join(base, "linux-arm64")
          )
        : firstExisting(
            path.join(base, "linux-x64-gnu"),
            path.join(base, "linux-x86_64")
          );
  } else if (platform === "win32") {
    dir = firstExisting(
      path.join(base, "win32-x64"),
      path.join(base, "windows-x86_64")
    );
  } else {
    throw new Error(`Unsupported ${platform}/${arch}`);
  }
  return path.join(dir, file);
}

const addonPath = path.join(__dirname, "..", "build", "Release", "msut.node");
const native = require(addonPath);
if (typeof native.bind === "function") native.bind(platformLibPath(process));

export type PeakOptions = Partial<{
  integralThreshold: number;
  intensityThreshold: number;
  widthThreshold: number;
  noise: number;
  autoNoise: boolean | number;
  autoBaseline: boolean | number;
  baselineWindow: number;
  baselineWindowFactor: number;
  allowOverlap: boolean | number;
  windowSize: number;
  snRatio: number;
}>;

export type BaselineOptions = Partial<{
  baselineWindow: number;
  baselineWindowFactor: number;
}>;

export type Peak = {
  from: number;
  to: number;
  rt: number;
  integral: number;
  intensity: number;
  ratio: number;
  np: number;
};

export type Target = {
  id?: string;
  rt: number;
  mz: number;
  ranges: number;
};

export type ChromItem = {
  id?: string;
  idx?: number;
  index?: number;
  rt: number;
  window?: number;
  range?: number;
};

export type ChromPeak = {
  index?: number;
  id?: string;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

export type FindFeaturesOptions = {
  eic?: { ppmTolerance?: number; mzTolerance?: number };
  grid?: { start?: number; end?: number; stepSize?: number };
  findPeak?: PeakOptions;
  cores?: number;
};

export function packPeakOptions(opts?: PeakOptions): Buffer | undefined {
  if (!opts) return undefined;
  const b = Buffer.alloc(64);
  const f64 = (v: unknown, off: number) =>
    b.writeDoubleLE(Number.isFinite(v as number) ? Number(v) : NaN, off);
  const i32 = (v: unknown, off: number) =>
    b.writeInt32LE(
      typeof v === "boolean"
        ? v
          ? 1
          : 0
        : Number.isFinite(v as number)
        ? (v as number) | 0
        : 0,
      off
    );
  f64(opts.integralThreshold, 0);
  f64(opts.intensityThreshold, 8);
  i32(opts.widthThreshold, 16);
  f64(opts.noise, 24);
  i32(opts.autoNoise, 32);
  i32(opts.autoBaseline, 36);
  i32(opts.baselineWindow, 40);
  i32(opts.baselineWindowFactor, 44);
  i32(opts.allowOverlap, 48);
  i32(opts.windowSize, 52);
  f64(opts.snRatio, 56);
  return b;
}

function toBuffer(v: Uint8Array | ArrayBuffer): Buffer {
  return v instanceof Uint8Array
    ? Buffer.from(v)
    : Buffer.from(new Uint8Array(v));
}

export function parseMzML(data: Uint8Array | ArrayBuffer): Buffer {
  const buf = toBuffer(data);
  const fn = native.parseMzml || native.parseMzML;
  return fn(buf) as Buffer;
}

export function binToJson(bin: Uint8Array | ArrayBuffer): string {
  const b = toBuffer(bin);
  return native.binToJson(b) as string;
}

export function calculateEic(
  bin: Uint8Array | ArrayBuffer,
  targets: number,
  from: number,
  to: number,
  ppmTol = 20,
  mzTol = 0.005
) {
  const b = toBuffer(bin);
  const fn = native.calculateEic;
  return fn(b, +targets, from, to, ppmTol, mzTol) as {
    x: Float64Array;
    y: Float64Array;
  };
}

export function findPeaks(
  x: Float64Array,
  y: Float64Array,
  opts?: PeakOptions
) {
  const s = native.findPeaks(x, y, packPeakOptions(opts));
  return JSON.parse(s) as Peak[];
}

export function getPeak(
  x: Float64Array,
  y: Float64Array,
  rt: number,
  range: number,
  opts?: PeakOptions
) {
  const s = native.getPeak(x, y, rt, range, packPeakOptions(opts));
  return JSON.parse(s) as Peak;
}

export const findNoiseLevel = native.findNoiseLevel as (
  y: Float64Array
) => number;

export function getPeaksFromEic(
  bin: Uint8Array | ArrayBuffer,
  targets: Target[],
  fromLeft = 0.5,
  toRight = 0.5,
  options?: PeakOptions,
  cores = 1
) {
  const n = targets.length;
  const rts = new Float64Array(n);
  const mzs = new Float64Array(n);
  const rng = new Float64Array(n);
  const ids: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const t = targets[i];
    rts[i] = +t.rt;
    mzs[i] = +t.mz;
    rng[i] = +t.ranges;
    ids[i] = t.id ?? "";
  }
  const b = toBuffer(bin);
  const optBuf = packPeakOptions(options);
  const json = native.getPeaksFromEic(
    b,
    rts,
    mzs,
    rng,
    ids,
    +fromLeft,
    +toRight,
    optBuf,
    cores | 0
  ) as string;
  return JSON.parse(json) as Array<{
    id?: string;
    mz: number;
    ort: number;
    rt: number;
    from: number;
    to: number;
    intensity: number;
    integral: number;
    noise: number;
  }>;
}

export function getPeaksFromChrom(
  bin: Uint8Array | ArrayBuffer,
  items: ChromItem[],
  options?: PeakOptions,
  cores = 1
) {
  const n = items.length;
  const idxs = new Uint32Array(n);
  const rts = new Float64Array(n);
  const rng = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const it = items[i];
    const idx = Number.isFinite(it.idx)
      ? (it.idx as number)
      : Number.isFinite(it.index)
      ? (it.index as number)
      : -1;
    idxs[i] = idx != null && idx >= 0 ? idx >>> 0 : 0xffffffff;
    rts[i] = +it.rt;
    const win = it.window ?? it.range ?? 0;
    rng[i] = +win;
  }
  const b = toBuffer(bin);
  const optBuf = packPeakOptions(options);
  const json = native.getPeaksFromChrom(
    b,
    idxs,
    rts,
    rng,
    optBuf,
    cores | 0
  ) as string;
  return JSON.parse(json) as ChromPeak[];
}

export function calculateBaseline(
  y: Float64Array | ArrayLike<number>,
  options?: BaselineOptions
): Float64Array {
  const y64 =
    y instanceof Float64Array ? y : new Float64Array(y as ArrayLike<number>);
  const win = (options && (options.baselineWindow as any)) | 0;
  const fac = (options && (options.baselineWindowFactor as any)) | 0;
  return native.calculateBaseline(y64, win, fac) as Float64Array;
}

export type Feature = {
  mz: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
  ratio: number;
  np: number;
};

export function findFeatures(
  data: Uint8Array | ArrayBuffer,
  fromTo: { from: number; to: number },
  options: FindFeaturesOptions = {}
): Feature[] {
  const {
    eic = { mzTolerance: 0.0025, ppmTolerance: 5.0 },
    grid = { start: 20, end: 700, stepSize: 0.005 },
    findPeak = {},
    cores = 1,
  } = options;
  const { from, to } = fromTo;

  const b = toBuffer(data);

  const eicPpm =
    typeof eic.ppmTolerance === "number" &&
    Number.isFinite(eic.ppmTolerance) &&
    eic.ppmTolerance >= 0
      ? eic.ppmTolerance
      : NaN;

  const eicMz =
    typeof eic.mzTolerance === "number" &&
    Number.isFinite(eic.mzTolerance) &&
    eic.mzTolerance >= 0
      ? eic.mzTolerance
      : NaN;

  const gridStart =
    typeof grid.start === "number" && Number.isFinite(grid.start)
      ? grid.start
      : NaN;

  const gridEnd =
    typeof grid.end === "number" && Number.isFinite(grid.end) ? grid.end : NaN;

  const gridStep =
    typeof grid.stepSize === "number" &&
    Number.isFinite(grid.stepSize) &&
    grid.stepSize > 0
      ? grid.stepSize
      : NaN;

  console.log("", { stepSize: grid.stepSize, gridStep });

  const peakBuf = packPeakOptions(findPeak);

  const s = native.findFeatures(
    b,
    from,
    to,
    eicPpm,
    eicMz,
    gridStart,
    gridEnd,
    gridStep,
    peakBuf ?? null,
    cores
  ) as string;

  return JSON.parse(s) as Feature[];
}

module.exports = {
  parseMzML,
  binToJson,
  calculateEic,
  getPeak,
  findPeaks,
  findNoiseLevel,
  getPeaksFromEic,
  getPeaksFromChrom,
  calculateBaseline,
  findFeatures,
};
