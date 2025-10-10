import {
  makeApi,
  type Exports,
  type FindPeaksOptions,
  type Peak,
} from "./makeApi.js";
import type { MzML } from "./types/mzml.js";

declare const __INLINE__: boolean;
const INLINE = typeof __INLINE__ !== "undefined" && __INLINE__;
declare const __WASM_DATA_URL__: string | undefined;
declare var require: any;

let MOD: WebAssembly.Module | null = null;

const freshInstance = (userImports: WebAssembly.Imports = {}): Exports => {
  const td: TextDecoder =
    typeof TextDecoder !== "undefined"
      ? new TextDecoder("utf-8")
      : new (require("node:util").TextDecoder)("utf-8");
  let mem: WebAssembly.Memory | null = null;

  const userEnv = (userImports as any).env ?? {};
  const env = {
    ...userEnv,
    js_log:
      userEnv.js_log ??
      ((ptr: number, len: number) => {
        if (!mem) return;
        const bytes = new Uint8Array(mem.buffer, ptr, len);
        console.log(td.decode(bytes));
      }),
  };

  const imports: WebAssembly.Imports = { ...userImports, env };
  const instance = new WebAssembly.Instance(MOD!, imports);
  mem = (instance.exports as any).memory as WebAssembly.Memory;
  return makeApi(instance);
};

export function parseMzML(
  data: Uint8Array | ArrayBuffer,
  options: { slim?: boolean; json?: boolean; pretty?: boolean } = {}
): MzML | Uint8Array {
  if (!MOD) throw new Error("ms-utils WASM not initialized");
  return freshInstance({}).parseMzML(data, options as any) as any;
}

export type EicResult = { x: number[]; y: number[] };

export const calculateEic = (
  source: Uint8Array | ArrayBuffer,
  targetMass: number | string,
  fromTo: { from: number; to: number },
  options: { ppmTolerance?: number; mzTolerance?: number } = {}
): EicResult => {
  if (!MOD) throw new Error("ms-utils WASM not initialized");
  const bin =
    source instanceof Uint8Array
      ? source
      : new Uint8Array(source as ArrayBuffer);
  const mass =
    typeof targetMass === "string" ? Number(targetMass) : +targetMass;
  const { x, y } = freshInstance({}).calculateEic(
    bin,
    mass,
    fromTo.from,
    fromTo.to,
    options.ppmTolerance ?? 20,
    options.mzTolerance ?? 0.005
  );
  return { x: Array.from(x), y: Array.from(y) };
};

export const getPeak = (
  x: Float64Array | ArrayLike<number>,
  y: Float32Array | ArrayLike<number>,
  rt: number,
  range: number,
  options: FindPeaksOptions = {}
): Peak => {
  if (!MOD) throw new Error("ms-utils WASM not initialized");
  const x64 =
    x instanceof Float64Array ? x : new Float64Array(x as ArrayLike<number>);
  const y32 =
    y instanceof Float32Array ? y : new Float32Array(y as ArrayLike<number>);
  return freshInstance({}).getPeak(x64, y32, rt, range, options) as Peak;
};

export const findPeaks = (
  x: Float64Array | ArrayLike<number>,
  y: Float32Array | ArrayLike<number>,
  options: FindPeaksOptions = {}
): Peak[] => {
  if (!MOD) throw new Error("ms-utils WASM not initialized");
  const x64 =
    x instanceof Float64Array ? x : new Float64Array(x as ArrayLike<number>);
  const y32 =
    y instanceof Float32Array ? y : new Float32Array(y as ArrayLike<number>);
  if (x64.length !== y32.length)
    throw new Error(
      `findPeaks: x.length (${x64.length}) != y.length (${y32.length}).`
    );
  return freshInstance({}).findPeaks(x64, y32, options) as Peak[];
};

export const findNoiseLevel = (y: Float32Array): number => {
  if (!MOD) throw new Error("ms-utils WASM not initialized");
  return freshInstance({}).findNoiseLevel(y);
};

const dataUrlToBytes = (dataUrl: string): Uint8Array => {
  const b64 = dataUrl.split(",")[1] ?? "";
  if (typeof atob === "function") {
    const bin = atob(b64);
    const out = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
    return out;
  }
  return Uint8Array.from(Buffer.from(b64, "base64"));
};

const compileOnce = (bytes: Uint8Array) => {
  if (MOD) return;
  const ab = new ArrayBuffer(bytes.byteLength);
  new Uint8Array(ab).set(bytes);
  MOD = new WebAssembly.Module(ab);
};

(function init() {
  if (INLINE) {
    if (!__WASM_DATA_URL__)
      throw new Error("INLINE build missing __WASM_DATA_URL__");
    compileOnce(dataUrlToBytes(__WASM_DATA_URL__));
    return;
  }
  if (typeof process !== "undefined" && (process as any).versions?.node) {
    const fs = require("node:fs") as typeof import("node:fs");
    const path = require("node:path") as typeof import("node:path");
    const url = require("node:url") as typeof import("node:url");
    const here =
      typeof __dirname !== "undefined"
        ? __dirname
        : path.dirname(url.fileURLToPath((0, eval)("import.meta").url));
    compileOnce(new Uint8Array(fs.readFileSync(path.join(here, "msut.wasm"))));
    return;
  }
  throw new Error("Browser ESM requires the UMD inline build.");
})();
