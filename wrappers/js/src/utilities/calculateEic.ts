import { xFindClosestIndex } from "ml-spectra-processing";
import { MzML, Run } from "../types/mzml";

export function calculateEic(
  data: MzML,
  targetMass: number | string,
  fromTo: { from: number; to: number },
  options: { ppmTolerance?: number; mzTolerance?: number } = {}
) {
  const { from, to } = fromTo;
  const { ppmTolerance = 20, mzTolerance = 0.005 } = options;
  if (!targetMass) {
    throw new Error(
      "targetMass must be defined and a number or string of comma separated numbers"
    );
  }
  let targetMasses;
  if (!Number.isNaN(Number(targetMass))) {
    targetMasses = [Number(targetMass)];
  } else if (typeof targetMass === "string") {
    targetMasses = targetMass
      .split(/[ ,;\r\n\t]+/)
      .map((value) => Number(value));
  }
  for (const mass of targetMasses) {
    if (isNaN(mass)) {
      throw new Error(
        "targetMass must be defined and a number or string of comma separated numbers"
      );
    }
  }

  const { times, series } = getJson(data.run);

  const fromIndex = xFindClosestIndex(times, from);
  const toIndex = xFindClosestIndex(times, to);
  const massSpectra = series.ms.data.slice(fromIndex, toIndex);
  const x = times.slice(fromIndex, toIndex);
  const y = new Array(massSpectra.length).fill(0);
  for (const mass of targetMasses) {
    for (let i = 0; i < massSpectra.length; i++) {
      y[i] = filterSpectra(massSpectra[i], mass, { ppmTolerance, mzTolerance });
    }
  }
  return { x, y };
}

function filterSpectra(
  spectrum: [number[], number[]],
  target: number,
  options: { ppmTolerance?: number; mzTolerance?: number } = {}
) {
  const { ppmTolerance = 20, mzTolerance = 0.005 } = options;
  const tolerance = Math.max((ppmTolerance / 1e6) * target, mzTolerance);
  const minMZ = target - tolerance;
  const maxMZ = target + tolerance;
  let result = 0;
  for (let i = 0; i < spectrum[0].length; i++) {
    const value = spectrum[0][i];
    if (value >= minMZ && value <= maxMZ) {
      result = spectrum[1][i];
    }
  }
  return result;
}

function getJson(data: Run) {
  const { spectra, chromatograms } = data;
  const result = {
    series: {
      bpc: {
        dimension: 1,
        name: "bpc",
        meta: {},
        data: [],
      },
      tic: {
        dimension: 1,
        name: "tic",
        meta: {},
        data: [],
      },
      ms: {
        dimension: 1,
        name: "ms",
        meta: {},
        data: [],
      },
    },
    times: [],
  };

  for (let i = 0; i < spectra.length; i++) {
    // @ts-expect-error i know
    if (spectra[i].ms_level === 1) {
      // @ts-expect-error i know
      result.times.push(spectra[i].retention_time);
      result.series.tic.data.push(chromatograms[0].intensityArray[i]);
      result.series.bpc.data.push(chromatograms[1].intensityArray[i]);
      result.series.ms.data.push([
        spectra[i].mzArray,
        spectra[i].intensityArray,
      ]);
    }
  }

  return result;
}
