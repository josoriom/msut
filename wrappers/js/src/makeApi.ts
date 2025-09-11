import { makeParseMzML, ParseMzML } from "./utilities/parseMzML.js";
import { calculateEic } from "./utilities/calculateEic.js";

export type FindPeaksOptions = {
  integralThreshold?: number;
  intensityThreshold?: number;
  widthThreshold?: number;
  noise?: number;
  autoNoise?: boolean;
  allowOverlap?: boolean;
  windowSize?: number;
  snRatio?: number;
};

export type Peak = {
  from: number;
  to: number;
  rt: number;
  integral: number;
  intensity: number;
  ratio: number;
  np: number;
};

export type ChromPeakRow = {
  index: number;
  id: string;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

export type EicPeakRow = {
  id: string;
  mz: number;
  ort: number;
  rt: number;
  from: number;
  to: number;
  intensity: number;
  integral: number;
};

export interface Exports {
  parseMzML: ParseMzML;
  calculateEic: typeof calculateEic;

  getPeak: (
    x: Float64Array,
    y: Float32Array,
    rt: number,
    range: number,
    options?: FindPeaksOptions
  ) => Peak;

  getPeaksFromChrom: (
    bin: Uint8Array,
    items: { idx: number; rt: number; window: number }[],
    options?: FindPeaksOptions,
    cores?: number
  ) => ChromPeakRow[];

  getPeaksFromEic: (
    bin: Uint8Array,
    items: { id?: string; rt: number; mz: number; ranges: number }[],
    fromTo: { from: number; to: number },
    options?: FindPeaksOptions,
    cores?: number
  ) => EicPeakRow[];

  findPeaks: (
    x: Float64Array,
    y: Float32Array,
    options?: FindPeaksOptions
  ) => Peak[];
  scanForPeaks: (
    x: Float64Array,
    y: Float32Array,
    options?: { epsilon?: number; windowSize?: number }
  ) => Float64Array;
  findNoiseLevel: (y: Float32Array) => number;

  __debug: {
    memory: WebAssembly.Memory;
    exports: Record<string, any>;
    heapBytes: () => number;
  };
}

type ExportsLike = Record<string, any>;
const getExports = (obj: any): ExportsLike => {
  if (obj?.exports?.memory instanceof WebAssembly.Memory) return obj.exports;
  if (obj?.asm?.memory instanceof WebAssembly.Memory) return obj.asm;
  if (obj?.memory instanceof WebAssembly.Memory) return obj;
  throw new Error("makeApi: could not find WebAssembly exports");
};
function pickFn<T extends Function>(ex: ExportsLike, names: string[]): T {
  for (let i = 0; i < names.length; i++) {
    const f = ex[names[i]];
    if (typeof f === "function") return f as unknown as T;
  }
  throw new Error("makeApi: missing exports: " + names.join("|"));
}

const BUF_PAIR_BYTES = 8;
const SIZE_COPTS = 48;

export function makeApi(instanceOrModule: any): Exports {
  const ex = getExports(instanceOrModule);
  const memory: WebAssembly.Memory = ex.memory;

  const alloc: (n: number) => number = pickFn(ex, ["alloc"]);
  const free: (p: number, n: number) => void = pickFn(ex, ["free_", "free"]);

  const parse_mzml: (
    p: number,
    n: number,
    slim: number,
    outBinBuf: number
  ) => number = pickFn(ex, ["parse_mzml"]);

  const parse_mzml_to_json: (
    p: number,
    n: number,
    slim: number,
    outJsonBuf: number,
    outBlobBuf: number
  ) => number = pickFn(ex, ["parse_mzml_to_json"]);

  const find_peaks: (
    xPtr: number,
    yPtr: number,
    len: number,
    optionsPtr: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["find_peaks"]);

  const scan_for_peaks: (
    xPtr: number,
    yPtr: number,
    len: number,
    epsilon: number,
    windowSize: number,
    outBuf: number
  ) => number = pickFn(ex, ["scan_for_peaks"]);

  const get_peak: (
    xPtr: number,
    yPtr: number,
    len: number,
    rt: number,
    range: number,
    optionsPtr: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["get_peak"]);

  const C_get_peaks_from_chrom: (
    binPtr: number,
    binLen: number,
    idxsPtr: number,
    rtsPtr: number,
    rangesPtr: number,
    nItems: number,
    optionsPtr: number,
    cores: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["C_get_peaks_from_chrom"]);

  const get_peaks_from_eic: (
    binPtr: number,
    binLen: number,
    rtsPtr: number,
    mzsPtr: number,
    rangesPtr: number,
    idsOffPtr: number,
    idsLenPtr: number,
    idsBufPtr: number,
    idsBufLen: number,
    nItems: number,
    fromLeft: number,
    toRight: number,
    optionsPtr: number,
    cores: number,
    outJsonBuf: number
  ) => number = pickFn(ex, ["C_get_peaks_from_eic"]);

  const find_noise_level: (yPtr: number, len: number) => number = pickFn(ex, [
    "find_noise_level",
  ]);

  let HEAPU8 = new Uint8Array(memory.buffer);
  let HEAPDV = new DataView(memory.buffer);
  const refreshViews = () => {
    if (HEAPU8.buffer !== memory.buffer) {
      HEAPU8 = new Uint8Array(memory.buffer);
      HEAPDV = new DataView(memory.buffer);
    }
  };
  const readBuf = (bufPtr: number) => {
    refreshViews();
    const ptr = HEAPDV.getUint32(bufPtr + 0, true);
    const len = HEAPDV.getUint32(bufPtr + 4, true);
    return { ptr, len };
  };
  const heapWrite = (dstPtr: number, src: Uint8Array) => {
    refreshViews();
    HEAPU8.set(src, dstPtr);
  };
  const heapSlice = (ptr: number, len: number): Uint8Array => {
    refreshViews();
    const out = new Uint8Array(len);
    out.set(HEAPU8.subarray(ptr, ptr + len));
    return out;
  };

  const SCRATCH_A = alloc(BUF_PAIR_BYTES);
  const SCRATCH_JSON = alloc(BUF_PAIR_BYTES);
  const SCRATCH_BLOB = alloc(BUF_PAIR_BYTES);
  const SCRATCH_PEAKS = alloc(BUF_PAIR_BYTES);
  const SCRATCH_OPTS = alloc(SIZE_COPTS);

  refreshViews();

  const td: TextDecoder =
    typeof TextDecoder !== "undefined"
      ? new TextDecoder("utf-8")
      : new (require("node:util").TextDecoder)("utf-8");
  const te: TextEncoder =
    typeof TextEncoder !== "undefined"
      ? new TextEncoder()
      : new (require("node:util").TextEncoder)();

  const writeOptions = (ptr: number, o?: FindPeaksOptions) => {
    refreshViews();
    if (!o) {
      HEAPDV.setFloat64(ptr + 0, Number.NaN, true);
      HEAPDV.setFloat64(ptr + 8, Number.NaN, true);
      HEAPDV.setInt32(ptr + 16, 0, true);
      HEAPDV.setFloat64(ptr + 24, Number.NaN, true);
      HEAPDV.setInt32(ptr + 32, 0, true);
      HEAPDV.setInt32(ptr + 36, 0, true);
      HEAPDV.setInt32(ptr + 40, 0, true);
      HEAPDV.setInt32(ptr + 44, 0, true);
      return;
    }
    HEAPDV.setFloat64(
      ptr + 0,
      typeof o.integralThreshold === "number"
        ? o.integralThreshold
        : Number.NaN,
      true
    );
    HEAPDV.setFloat64(
      ptr + 8,
      typeof o.intensityThreshold === "number"
        ? o.intensityThreshold
        : Number.NaN,
      true
    );
    HEAPDV.setInt32(
      ptr + 16,
      typeof o.widthThreshold === "number" ? o.widthThreshold | 0 : 0,
      true
    );
    HEAPDV.setFloat64(
      ptr + 24,
      typeof o.noise === "number" ? o.noise : Number.NaN,
      true
    );
    HEAPDV.setInt32(ptr + 32, o.autoNoise ? 1 : 0, true);
    HEAPDV.setInt32(ptr + 36, o.allowOverlap ? 1 : 0, true);
    HEAPDV.setInt32(
      ptr + 40,
      typeof o.windowSize === "number" ? o.windowSize | 0 : 0,
      true
    );
    HEAPDV.setInt32(
      ptr + 44,
      typeof o.snRatio === "number" ? o.snRatio | 0 : 0,
      true
    );
  };

  const parseMzML = makeParseMzML({
    alloc,
    free,
    refreshViews,
    heapWrite,
    readBuf,
    heapSlice,
    parse_mzml,
    SCRATCH_A,
    parse_mzml_to_json,
    SCRATCH_JSON,
    SCRATCH_BLOB,
    td,
  });

  const getPeak = (
    x: Float64Array,
    y: Float32Array,
    rt: number,
    range: number,
    options?: FindPeaksOptions
  ): Peak => {
    if (!(x instanceof Float64Array))
      throw new Error("getPeak: x must be Float64Array");
    if (!(y instanceof Float32Array))
      throw new Error("getPeak: y must be Float32Array");
    if (x.length !== y.length || x.length < 3)
      throw new Error("getPeak: x,y length mismatch (>=3)");

    const xb = new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
    const yb = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const xp = alloc(xb.length);
    const yp = alloc(yb.length);
    heapWrite(xp, xb);
    heapWrite(yp, yb);

    let pOpts = 0;
    if (options) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }

    const rc = get_peak(xp, yp, x.length, rt, range, pOpts, SCRATCH_JSON);
    free(xp, xb.length);
    free(yp, yb.length);
    if (rc !== 0) throw new Error(`get_peak failed: ${rc}`);

    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const getPeaksFromChrom = (
    bin: Uint8Array,
    items: { idx: number; rt: number; window: number }[],
    options?: FindPeaksOptions,
    cores = 1
  ): ChromPeakRow[] => {
    const n = items.length;
    if (n === 0) return [];

    const idxs = new Uint32Array(n);
    const rts = new Float64Array(n);
    const wins = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const it = items[i];
      const idx = it.idx;
      idxs[i] = idx == null || idx < 0 ? 0xffffffff : idx >>> 0;
      rts[i] = it.rt;
      wins[i] = it.window;
    }

    const binPtr = alloc(bin.length);
    heapWrite(binPtr, bin);

    const idxsU8 = new Uint8Array(
      idxs.buffer,
      idxs.byteOffset,
      idxs.byteLength
    );
    const rtsU8 = new Uint8Array(rts.buffer, rts.byteOffset, rts.byteLength);
    const winU8 = new Uint8Array(wins.buffer, wins.byteOffset, wins.byteLength);

    const idxsPtr = alloc(idxsU8.length);
    const rtsPtr = alloc(rtsU8.length);
    const winPtr = alloc(winU8.length);
    heapWrite(idxsPtr, idxsU8);
    heapWrite(rtsPtr, rtsU8);
    heapWrite(winPtr, winU8);

    let pOpts = 0;
    if (options) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }

    const rc = C_get_peaks_from_chrom(
      binPtr,
      bin.length,
      idxsPtr,
      rtsPtr,
      winPtr,
      n,
      pOpts,
      Math.max(1, cores | 0),
      SCRATCH_JSON
    );

    free(binPtr, bin.length);
    free(idxsPtr, idxsU8.length);
    free(rtsPtr, rtsU8.length);
    free(winPtr, winU8.length);

    if (rc !== 0) throw new Error(`C_get_peaks_from_chrom failed: ${rc}`);
    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const buildIdTable = (ids: (string | undefined | null)[]) => {
    const offs = new Uint32Array(ids.length);
    const lens = new Uint32Array(ids.length);
    const chunks: Uint8Array[] = [];
    let total = 0;
    for (let i = 0; i < ids.length; i++) {
      const s = (ids[i] ?? "") + "";
      const b = te.encode(s);
      offs[i] = total >>> 0;
      lens[i] = b.length >>> 0;
      chunks.push(b);
      total += b.length;
    }
    const buf = new Uint8Array(total);
    let pos = 0;
    for (const c of chunks) {
      buf.set(c, pos);
      pos += c.length;
    }
    return { offs, lens, buf };
  };

  type EicItem = {
    id?: string;
    rt: number;
    mz: number;
    range?: number; // preferred name
    sd?: number; // accepted alias
    window?: number; // accepted alias
  };

  const TE: TextEncoder =
    typeof TextEncoder !== "undefined"
      ? new TextEncoder()
      : new (require("node:util").TextEncoder)();

  const getPeaksFromEic = (
    bin: Uint8Array,
    items: EicItem[],
    window: { from: number; to: number },
    options?: FindPeaksOptions,
    cores: number = 1
  ) => {
    const n = items.length;
    if (n === 0) return [];

    // Build numeric arrays
    const rts = new Float64Array(n);
    const mzs = new Float64Array(n);
    const ranges = new Float64Array(n);

    // Build id buffer + offsets/lengths
    const idStrings = new Array<string>(n);
    for (let i = 0; i < n; i++) {
      const it = items[i] ?? ({} as any);
      rts[i] = Number(it.rt) || 0;
      mzs[i] = Number(it.mz) || 0;
      const rng =
        (typeof it.range === "number" && it.range > 0 ? it.range : undefined) ??
        (typeof it.sd === "number" && it.sd > 0 ? it.sd : undefined) ??
        (typeof it.window === "number" && it.window > 0
          ? it.window
          : undefined) ??
        0.25; // default fallback
      ranges[i] = rng;
      idStrings[i] = typeof it.id === "string" ? it.id : "";
    }
    // UTF-8 pack ids
    const encodedIds = idStrings.map((s) => TE.encode(s));
    let idsBufLen = 0;
    for (const b of encodedIds) idsBufLen += b.length;
    const idsBufU8 = new Uint8Array(idsBufLen);
    const idsOff = new Uint32Array(n);
    const idsLen = new Uint32Array(n);
    {
      let cur = 0;
      for (let i = 0; i < n; i++) {
        const b = encodedIds[i];
        idsOff[i] = cur >>> 0;
        idsLen[i] = b.length >>> 0;
        idsBufU8.set(b, cur);
        cur += b.length;
      }
    }

    // Write inputs to wasm heap
    const binPtr = alloc(bin.length);
    heapWrite(binPtr, bin);

    const rtsU8 = new Uint8Array(rts.buffer, rts.byteOffset, rts.byteLength);
    const mzsU8 = new Uint8Array(mzs.buffer, mzs.byteOffset, mzs.byteLength);
    const rngU8 = new Uint8Array(
      ranges.buffer,
      ranges.byteOffset,
      ranges.byteLength
    );
    const offU8 = new Uint8Array(
      idsOff.buffer,
      idsOff.byteOffset,
      idsOff.byteLength
    );
    const lenU8 = new Uint8Array(
      idsLen.buffer,
      idsLen.byteOffset,
      idsLen.byteLength
    );

    const rtsPtr = alloc(rtsU8.length);
    const mzsPtr = alloc(mzsU8.length);
    const rngPtr = alloc(rngU8.length);
    const offPtr = alloc(offU8.length);
    const lenPtr = alloc(lenU8.length);
    const idbPtr = alloc(idsBufU8.length);

    heapWrite(rtsPtr, rtsU8);
    heapWrite(mzsPtr, mzsU8);
    heapWrite(rngPtr, rngU8);
    heapWrite(offPtr, offU8);
    heapWrite(lenPtr, lenU8);
    heapWrite(idbPtr, idsBufU8);

    let pOpts = 0;
    if (options && Object.keys(options).length) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }

    const rc = get_peaks_from_eic(
      binPtr,
      bin.length,
      rtsPtr,
      mzsPtr,
      rngPtr,
      offPtr,
      lenPtr,
      idbPtr,
      idsBufU8.length,
      n,
      window.from,
      window.to,
      pOpts,
      cores,
      SCRATCH_JSON
    );

    free(binPtr, bin.length);
    free(rtsPtr, rtsU8.length);
    free(mzsPtr, mzsU8.length);
    free(rngPtr, rngU8.length);
    free(offPtr, offU8.length);
    free(lenPtr, lenU8.length);
    free(idbPtr, idsBufU8.length);

    if (rc !== 0) throw new Error(`get_peaks_from_eic failed: ${rc}`);

    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);

    const out = JSON.parse(td.decode(bytes));
    return out;
  };

  const findPeaks = (
    x: Float64Array,
    y: Float32Array,
    options: FindPeaksOptions = {}
  ) => {
    if (!(x instanceof Float64Array))
      throw new Error("findPeaks: x must be Float64Array");
    if (!(y instanceof Float32Array))
      throw new Error("findPeaks: y must be Float32Array");
    if (x.length !== y.length)
      throw new Error("findPeaks: x,y length mismatch");

    const xb = new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
    const yb = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const xp = alloc(xb.length);
    const yp = alloc(yb.length);
    heapWrite(xp, xb);
    heapWrite(yp, yb);

    let pOpts = 0;
    if (options && Object.keys(options).length > 0) {
      writeOptions(SCRATCH_OPTS, options);
      pOpts = SCRATCH_OPTS;
    }

    const rc = find_peaks(xp, yp, x.length, pOpts, SCRATCH_JSON);
    free(xp, xb.length);
    free(yp, yb.length);
    if (rc !== 0) throw new Error(`find_peaks failed: ${rc}`);

    const { ptr, len } = readBuf(SCRATCH_JSON);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return JSON.parse(td.decode(bytes));
  };

  const scanForPeaks = (
    x: Float64Array,
    y: Float32Array,
    options: { epsilon?: number; windowSize?: number } = {}
  ) => {
    const { epsilon = 1e-5, windowSize = 15 } = options;
    if (!(x instanceof Float64Array))
      throw new Error("scanForPeaks: x must be Float64Array");
    if (!(y instanceof Float32Array))
      throw new Error("scanForPeaks: y must be Float32Array");
    if (x.length !== y.length)
      throw new Error("scanForPeaks: x,y length mismatch");

    const xb = new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
    const yb = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const xp = alloc(xb.length);
    const yp = alloc(yb.length);
    heapWrite(xp, xb);
    heapWrite(yp, yb);

    const rc = scan_for_peaks(
      xp,
      yp,
      x.length,
      epsilon,
      windowSize | 0,
      SCRATCH_PEAKS
    );
    free(xp, xb.length);
    free(yp, yb.length);
    if (rc !== 0) throw new Error(`scan_for_peaks failed: ${rc}`);

    const { ptr, len } = readBuf(SCRATCH_PEAKS);
    const bytes = heapSlice(ptr, len);
    free(ptr, len);
    return new Float64Array(bytes.buffer, bytes.byteOffset, len / 8).slice();
  };

  const findNoiseLevel = (y: Float32Array) => {
    const bytes = new Uint8Array(y.buffer, y.byteOffset, y.byteLength);
    const p = alloc(bytes.length);
    heapWrite(p, bytes);
    const noise = find_noise_level(p, y.length);
    free(p, bytes.length);
    return noise;
  };

  return {
    parseMzML,
    calculateEic,
    getPeak,
    getPeaksFromChrom,
    getPeaksFromEic,
    findPeaks,
    scanForPeaks,
    findNoiseLevel,
    __debug: { memory, exports: ex, heapBytes: () => memory.buffer.byteLength },
  };
}
