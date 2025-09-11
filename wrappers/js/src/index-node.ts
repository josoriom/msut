import * as fs from "fs";
import * as path from "path";

const addonPath = path.join(__dirname, "..", "build", "Release", "msut.node");
const native = require(addonPath);

function firstExisting(...candidates: string[]) {
  for (const p of candidates) if (fs.existsSync(p)) return p;
  return candidates[0];
}
function platformLibPath(): string {
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

if (typeof native.bind === "function") native.bind(platformLibPath());

export type Peak = {
  id?: string;
  mz: number;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
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
  ranges: number;
};

export type ChromPeak = Peak & { index?: number };

export function parseMzML(
  data: Uint8Array | ArrayBuffer,
  opts: { slim?: boolean; json?: boolean } = {}
) {
  const buf =
    data instanceof Uint8Array
      ? Buffer.from(data)
      : Buffer.from(new Uint8Array(data));
  if (opts.json) {
    const { json, blob } = native.parseMzMLToJson(buf, !!opts.slim);
    const obj = JSON.parse(json);
    Object.defineProperty(obj, "__bin1", { value: blob, enumerable: false });
    return obj;
  }
  return native.parseMzML(buf, !!opts.slim) as Buffer;
}

export function binToJson(
  bin: Uint8Array | ArrayBuffer,
  pretty = false
): string {
  const b =
    bin instanceof Uint8Array
      ? Buffer.from(bin)
      : Buffer.from(new Uint8Array(bin));
  return native.binToJson(b, !!pretty);
}

export function binToEic(
  bin: Uint8Array | ArrayBuffer,
  targets: string,
  from: number,
  to: number,
  ppmTol = 0,
  mzTol = 0
) {
  const b =
    bin instanceof Uint8Array
      ? Buffer.from(bin)
      : Buffer.from(new Uint8Array(bin));
  return native.binToEic(b, targets, from, to, ppmTol, mzTol) as {
    x: Float64Array;
    y: Float32Array;
  };
}

export function getPeak(
  x: Float64Array,
  y: Float32Array,
  rt: number,
  range: number,
  opts?: Record<string, unknown>
) {
  const s = native.getPeak(x, y, rt, range, opts);
  return JSON.parse(s);
}

export function findPeaks(
  x: Float64Array,
  y: Float32Array,
  opts?: Record<string, unknown>
) {
  const s = native.findPeaks(x, y, opts);
  return JSON.parse(s);
}

export const findNoiseLevel = native.findNoiseLevel as (
  y: Float32Array
) => number;

export function getPeaksFromEic(
  bin: Uint8Array | ArrayBuffer,
  targets: Target[],
  fromTo: { from: number; to: number } = { from: 0.5, to: 5 },
  options?: Record<string, unknown>,
  cores = 1
): Peak[] {
  const { from, to } = fromTo;
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
  const b =
    bin instanceof Uint8Array
      ? Buffer.from(bin)
      : Buffer.from(new Uint8Array(bin));
  const json = native.getPeaksFromEic(
    b,
    rts,
    mzs,
    rng,
    ids,
    from,
    to,
    options,
    cores
  ) as string;
  return JSON.parse(json) as Peak[];
}

export function getPeaksFromChrom(
  bin: Uint8Array | ArrayBuffer,
  items: ChromItem[],
  options?: Record<string, unknown>,
  cores = 1
): ChromPeak[] {
  const n = items.length;
  const idxs = new Uint32Array(n);
  const rts = new Float64Array(n);
  const rng = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const it = items[i];
    const idx = Number.isFinite(it.idx)
      ? it.idx!
      : Number.isFinite(it.index)
      ? (it.index as number)
      : -1;
    idxs[i] = idx != null && idx >= 0 ? idx >>> 0 : 0xffffffff;
    rts[i] = +it.rt;
    rng[i] = +it.ranges;
  }
  const b =
    bin instanceof Uint8Array
      ? Buffer.from(bin)
      : Buffer.from(new Uint8Array(bin));
  const json = native.getPeaksFromChrom(
    b,
    idxs,
    rts,
    rng,
    options,
    cores
  ) as string;
  return JSON.parse(json) as ChromPeak[];
}
module.exports = {
  parseMzML,
  binToJson,
  binToEic,
  getPeak,
  findPeaks,
  findNoiseLevel: native.findNoiseLevel as (y: Float32Array) => number,
  getPeaksFromEic,
  getPeaksFromChrom,
};
