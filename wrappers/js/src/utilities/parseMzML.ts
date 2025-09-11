import type { MzML } from "../types/mzml";

export type ParseMzML = {
  (
    data: Uint8Array | ArrayBuffer,
    options?: { slim?: boolean; json?: false }
  ): Uint8Array;
  (
    data: Uint8Array | ArrayBuffer,
    options: { slim?: boolean; json: true; pretty?: boolean }
  ): MzML;
};

export function makeParseMzML(deps: {
  alloc: (n: number) => number;
  free: (p: number, n: number) => void;
  refreshViews: () => void;
  heapWrite: (dstPtr: number, src: Uint8Array) => void;
  readBuf: (bufPtr: number) => { ptr: number; len: number };
  heapSlice: (ptr: number, len: number) => Uint8Array;
  parse_mzml: (p: number, n: number, slim: number, outBinBuf: number) => number;
  SCRATCH_A: number;
  parse_mzml_to_json: (
    p: number,
    n: number,
    slim: number,
    outJsonBuf: number,
    outBlobBuf: number
  ) => number;
  SCRATCH_JSON: number;
  SCRATCH_BLOB: number;
  td: TextDecoder;
}): ParseMzML {
  const {
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
  } = deps;

  const impl = (
    data: Uint8Array | ArrayBuffer,
    options: { slim?: boolean; json?: boolean; pretty?: boolean } = {}
  ): Uint8Array | MzML => {
    const { slim = true, json = true } = options;
    const bytes = data instanceof Uint8Array ? data : new Uint8Array(data);

    const inPtr = alloc(bytes.byteLength);
    refreshViews();
    heapWrite(inPtr, bytes);

    try {
      if (json) {
        // Parse XML -> (meta JSON + BIN blob). Your Rust returns a BIN1 blob here.
        const rc = parse_mzml_to_json(
          inPtr,
          bytes.byteLength,
          slim ? 1 : 0,
          SCRATCH_JSON,
          SCRATCH_BLOB
        );
        refreshViews();
        if (rc !== 0) throw new Error(`parse_mzml_to_json failed (rc=${rc})`);

        const { ptr: jsonPtr, len: jsonLen } = readBuf(SCRATCH_JSON);
        const { ptr: blobPtr, len: blobLen } = readBuf(SCRATCH_BLOB);
        if (!jsonPtr || !jsonLen)
          throw new Error("parse_mzml_to_json: empty JSON");
        if (!blobPtr || !blobLen)
          throw new Error("parse_mzml_to_json: empty BLOB");

        const jsonBytes = heapSlice(jsonPtr, jsonLen);
        const blobBytes = heapSlice(blobPtr, blobLen);

        free(jsonPtr, jsonLen);
        free(blobPtr, blobLen);
        refreshViews();

        const meta = JSON.parse(td.decode(jsonBytes)) as MzML;
        return mergeBlobIntoJson(meta, blobBytes);
      } else {
        // Request raw BIN (BIN1) only.
        const rc = parse_mzml(inPtr, bytes.byteLength, slim ? 1 : 0, SCRATCH_A);
        refreshViews();
        if (rc !== 0) throw new Error(`parse_mzml failed (rc=${rc})`);

        const { ptr: binPtr, len: binLen } = readBuf(SCRATCH_A);
        if (!binPtr || !binLen) throw new Error("parse_mzml: empty BIN output");

        const out = heapSlice(binPtr, binLen);
        free(binPtr, binLen);
        refreshViews();
        return out;
      }
    } finally {
      free(inPtr, bytes.byteLength);
      refreshViews();
    }
  };

  return impl as ParseMzML;
}

/**
 * Attach typed-array views from a BIN blob to the parsed JSON `meta`.
 * Supports BOTH formats:
 *  - "BINS" (arrays-only)
 *  - "BIN1" (full container with meta sections)
 *
 * We only need the index tables, which are identical (32 bytes per entry).
 */
export function mergeBlobIntoJson(meta: MzML, blob: Uint8Array): MzML {
  if (!meta?.run) return meta;

  const base = blob.byteOffset;
  const buf = blob.buffer;
  const end = base + blob.byteLength;
  const dv = new DataView(buf, base, blob.byteLength);

  const rdU32 = (p: number) => dv.getUint32(p, true);
  const rdU64 = (p: number) => {
    if ((dv as any).getBigUint64) return Number(dv.getBigUint64(p, true));
    const lo = dv.getUint32(p + 0, true);
    const hi = dv.getUint32(p + 4, true);
    return hi * 0x1_0000_0000 + lo;
  };

  // Accept both BINS and BIN1
  const m0 = dv.getUint8(0),
    m1 = dv.getUint8(1),
    m2 = dv.getUint8(2),
    m3 = dv.getUint8(3);
  const MAGIC = String.fromCharCode(m0, m1, m2, m3);
  if (MAGIC !== "BINS" && MAGIC !== "BIN1") {
    throw new Error(`unexpected blob magic: ${MAGIC}`);
  }

  // Header fields (shared layout)
  const nSpec = rdU32(4);
  const nCh = rdU32(8);

  const chromXFmt = dv.getUint8(12);
  const chromYFmt = dv.getUint8(13);
  const spectXFmt = dv.getUint8(14);
  const spectYFmt = dv.getUint8(15);

  const specIndexOff = rdU64(16);
  const chromIndexOff = rdU64(24);
  // For BIN1 there are extra meta offsets at 32/40; we don't need them to read arrays.

  const oob = (byteOff: number, bytes: number) => {
    const need = base + byteOff + bytes;
    if (need > end) throw new RangeError("typed view OOB");
  };
  const viewF32 = (off: number, len: number) => {
    oob(off, len * 4);
    return new Float32Array(buf, base + off, len);
  };
  const viewF64 = (off: number, len: number) => {
    oob(off, len * 8);
    return new Float64Array(buf, base + off, len);
  };

  // Spectra arrays
  for (let i = 0; i < nSpec && i < meta.run.spectra.length; i++) {
    const idxBase = specIndexOff + i * 32;
    const mzOff = rdU64(idxBase + 0);
    const mzLen = rdU32(idxBase + 8);
    const inOff = rdU64(idxBase + 12);
    const inLen = rdU32(idxBase + 20);
    const s = meta.run.spectra[i];

    if (mzOff && mzLen) {
      s.mzArray =
        spectXFmt === 2
          ? viewF64(mzOff, mzLen)
          : spectXFmt === 1
          ? viewF32(mzOff, mzLen)
          : (() => {
              throw new Error(`unknown spect_x_fmt=${spectXFmt}`);
            })();
    }
    if (inOff && inLen) {
      s.intensityArray =
        spectYFmt === 1
          ? viewF32(inOff, inLen)
          : spectYFmt === 2
          ? viewF64(inOff, inLen)
          : (() => {
              throw new Error(`unknown spect_y_fmt=${spectYFmt}`);
            })();
    }
  }

  // Chromatogram arrays
  for (let i = 0; i < nCh && i < meta.run.chromatograms.length; i++) {
    const idxBase = chromIndexOff + i * 32;
    const tOff = rdU64(idxBase + 0);
    const tLen = rdU32(idxBase + 8);
    const yOff = rdU64(idxBase + 12);
    const yLen = rdU32(idxBase + 20);
    const c = meta.run.chromatograms[i];

    if (tOff && tLen) {
      c.timeArray =
        chromXFmt === 2
          ? viewF64(tOff, tLen)
          : chromXFmt === 1
          ? viewF32(tOff, tLen)
          : (() => {
              throw new Error(`unknown chrom_x_fmt=${chromXFmt}`);
            })();
    }
    if (yOff && yLen) {
      c.intensityArray =
        chromYFmt === 1
          ? viewF32(yOff, yLen)
          : chromYFmt === 2
          ? viewF64(yOff, yLen)
          : (() => {
              throw new Error(`unknown chrom_y_fmt=${chromYFmt}`);
            })();
    }
  }

  return meta;
}
