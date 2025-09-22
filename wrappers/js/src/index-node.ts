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
  allowOverlap: boolean | number;
  windowSize: number;
  snRatio: number;
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

export function packPeakOptions(opts?: PeakOptions): Buffer | undefined {
  if (!opts) return undefined;
  const b = Buffer.alloc(56);
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
  i32(opts.allowOverlap, 36);
  i32(opts.windowSize, 40);
  f64(opts.snRatio, 48);
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
  targets: string,
  from: number,
  to: number,
  ppmTol = 20,
  mzTol = 0.005
) {
  const b = toBuffer(bin);
  const fn = native.calculateEic;
  return fn(b, targets, from, to, ppmTol, mzTol) as {
    x: Float64Array;
    y: Float32Array;
  };
}

export function findPeaks(
  x: Float64Array,
  y: Float32Array,
  opts?: PeakOptions
) {
  const s = native.findPeaks(x, y, packPeakOptions(opts));
  return JSON.parse(s) as Peak[];
}

export function getPeak(
  x: Float64Array,
  y: Float32Array,
  rt: number,
  range: number,
  opts?: PeakOptions
) {
  const s = native.getPeak(x, y, rt, range, packPeakOptions(opts));
  return JSON.parse(s) as Peak;
}

export const findNoiseLevel = native.findNoiseLevel as (
  y: Float32Array
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

module.exports = {
  parseMzML,
  binToJson,
  calculateEic,
  getPeak,
  findPeaks,
  findNoiseLevel,
  getPeaksFromEic,
  getPeaksFromChrom,
};
