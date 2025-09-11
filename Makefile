CRATE          := msut
CRATE_MANIFEST := core/Cargo.toml
ARTIFACTS      := artifacts
DOCKER_BULLSEYE_IMAGE   := rust:1-bullseye
WINDOWS_ARM64_IMAGE := dockcross/windows-arm64:latest


.PHONY: all macos-arm64 macos-x86_64 linux-amd64 linux-arm64 wasm rstage windows-amd64 windows-arm64

all: mac linux windows wasm rstage

mac: macos-arm64 macos-x86_64
linux: linux-amd64 linux-arm64
windows: windows-amd64

macos-arm64:
	rustup target add aarch64-apple-darwin
	cargo build --manifest-path $(CRATE_MANIFEST) --release --target aarch64-apple-darwin
	mkdir -p $(ARTIFACTS)/macos-arm64
	cp core/target/aarch64-apple-darwin/release/lib$(CRATE).dylib $(ARTIFACTS)/macos-arm64/

macos-x86_64:
	rustup target add x86_64-apple-darwin
	cargo build --manifest-path $(CRATE_MANIFEST) --release --target x86_64-apple-darwin
	mkdir -p $(ARTIFACTS)/macos-x86_64
	cp core/target/x86_64-apple-darwin/release/lib$(CRATE).dylib $(ARTIFACTS)/macos-x86_64/

linux-amd64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-linux-amd64 \
	  -v $$PWD:/work -w /work \
	  --entrypoint /usr/local/cargo/bin/cargo $(DOCKER_BULLSEYE_IMAGE) \
	  build --manifest-path $(CRATE_MANIFEST) --release
	mkdir -p $(ARTIFACTS)/linux-x86_64
	cp core/target-linux-amd64/release/lib$(CRATE).so $(ARTIFACTS)/linux-x86_64/

linux-arm64:
	docker run --rm --platform=linux/arm64 \
	  -e CARGO_TARGET_DIR=/work/core/target-linux-arm64 \
	  -v $$PWD:/work -w /work \
	  --entrypoint /usr/local/cargo/bin/cargo $(DOCKER_BULLSEYE_IMAGE) \
	  build --manifest-path $(CRATE_MANIFEST) --release
	mkdir -p $(ARTIFACTS)/linux-arm64
	cp core/target-linux-arm64/release/lib$(CRATE).so $(ARTIFACTS)/linux-arm64/

wasm:
	rustup target add wasm32-unknown-unknown
	cargo build --manifest-path $(CRATE_MANIFEST) --release --target wasm32-unknown-unknown
	mkdir -p $(ARTIFACTS)/wasm
	cp core/target/wasm32-unknown-unknown/release/$(CRATE).wasm $(ARTIFACTS)/wasm/
	cp core/target/wasm32-unknown-unknown/release/$(CRATE).wasm wrappers/js/src/ 2>/dev/null || true

rstage:
	mkdir -p wrappers/r/inst/libs
	@set -e; \
	for d in $(ARTIFACTS)/macos-* $(ARTIFACTS)/linux-* $(ARTIFACTS)/windows-* ; do \
	  [ -d "$$d" ] || continue; \
	  base=$$(basename "$$d"); \
	  mkdir -p "wrappers/r/inst/libs/$$base"; \
	  cp -f "$$d"/* "wrappers/r/inst/libs/$$base/"; \
	done

clean:
	cargo clean --manifest-path $(CRATE_MANIFEST)
	rm -rf $(ARTIFACTS) wrappers/r/inst/libs core/target-linux-amd64 core/target-linux-arm64

cross-check:
	@command -v cross >/dev/null || (echo "Please install cross: cargo install cross" && false)

windows-amd64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-windows-amd64 \
	  -v $$PWD:/work -w /work $(DOCKER_BULLSEYE_IMAGE) bash -lc '\
	    set -euo pipefail; \
	    export PATH=/usr/local/cargo/bin:$$PATH; \
	    apt-get update && apt-get install -y --no-install-recommends \
	      curl ca-certificates gcc-mingw-w64-x86-64 g++-mingw-w64-x86-64 mingw-w64; \
	    if ! command -v rustup >/dev/null 2>&1; then \
	      curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; \
	    fi; \
	    rustup target add x86_64-pc-windows-gnu; \
	    RUSTFLAGS="-C linker=x86_64-w64-mingw32-gcc" \
	      cargo build --manifest-path $(CRATE_MANIFEST) --release --target x86_64-pc-windows-gnu'
	mkdir -p $(ARTIFACTS)/windows-x86_64
	cp core/target-windows-amd64/x86_64-pc-windows-gnu/release/$(CRATE).dll $(ARTIFACTS)/windows-x86_64/

WINDOWS_ARM64_IMAGE := dockcross/windows-arm64:latest
WINDOWS_ARM64_TARGET := aarch64-pc-windows-gnullvm

windows-arm64:
	docker run --rm --platform=linux/amd64 \
	  -e CARGO_TARGET_DIR=/work/core/target-windows-arm64 \
	  -v $$PWD:/work -w /work $(WINDOWS_ARM64_IMAGE) \
	  bash -lc 'set -euo pipefail; \
	    export PATH="$$HOME/.cargo/bin:/usr/local/cargo/bin:$$PATH"; \
	    command -v rustup >/dev/null 2>&1 || { curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; . "$$HOME/.cargo/env"; }; \
	    rustup target add $(WINDOWS_ARM64_TARGET); \
	    cargo build --manifest-path $(CRATE_MANIFEST) --release --target $(WINDOWS_ARM64_TARGET)'
	mkdir -p $(ARTIFACTS)/windows-arm64
	cp core/target-windows-arm64/$(WINDOWS_ARM64_TARGET)/release/$(CRATE).dll $(ARTIFACTS)/windows-arm64/
