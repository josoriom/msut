// src/types/assets.d.ts
declare module "*.wasm?asset" {
  const assetUrl: string; // data URL when using asset/inline
  export default assetUrl;
}

declare module "*.wasm" {
  const assetUrl: string; // file URL/path when using asset/resource
  export default assetUrl;
}
