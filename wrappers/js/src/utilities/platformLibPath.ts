/// <reference types="node" />
import * as path from "path";
import * as fs from "fs";

function firstExisting(...candidates: string[]) {
  for (const p of candidates) if (fs.existsSync(p)) return p;
  return candidates[0];
}

export function platformLibPath(proc: NodeJS.Process): string {
  const base = path.join(__dirname, "..", "native");
  const { platform, arch } = proc;
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
