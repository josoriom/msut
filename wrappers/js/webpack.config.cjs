"use strict";
const path = require("node:path");
const fs = require("node:fs");
const webpack = require("webpack");

const wasmPath = path.resolve(__dirname, "lib-esm/msut.wasm");

let wasmDataUrl;
try {
  const bytes = fs.readFileSync(wasmPath);
  wasmDataUrl = "data:application/wasm;base64," + bytes.toString("base64");
} catch (e) {
  throw new Error(
    `Could not read ${wasmPath}. Run your TS build + copy step first. ` +
      `Original error: ${e.message}`
  );
}

module.exports = {
  entry: "./lib-esm/index-wasm.js",
  mode: "production",
  output: {
    path: path.resolve(__dirname, "dist"),
    filename: "msut.js",
    library: "msut",
    libraryTarget: "umd",
    globalObject: "this",
  },
  plugins: [
    new webpack.DefinePlugin({
      __INLINE__: "true",
      __WASM_DATA_URL__: JSON.stringify(wasmDataUrl),
    }),
    new webpack.IgnorePlugin({ resourceRegExp: /^node:fs\/promises$/ }),
  ],
  performance: { maxAssetSize: 2_000_000, maxEntrypointSize: 2_000_000 },
};
