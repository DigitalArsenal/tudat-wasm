# tudat-wasm

WebAssembly build of the [Tudat](https://docs.tudat.space/) astrodynamics library, enabling orbital mechanics simulations directly in web browsers and Node.js.

## Overview

This repository provides:
- **WASM bindings** for the Tudat C++ astrodynamics library via Emscripten/Embind
- **NPM package** (`@tudat/tudatpy-wasm`) for JavaScript/TypeScript integration
- **Interactive demos** hosted on GitHub Pages
- **Test suite** for both Node.js and browser environments

## Quick Start

### Prerequisites

- CMake 3.20+
- Ninja build system
- Git
- Python 3.x
- Node.js (for testing)

The Emscripten SDK is downloaded automatically during the first build.

### Building

```bash
# Clone the repository
git clone https://github.com/tudat-team/tudat-wasm.git
cd tudat-wasm

# Build visualization demo (faster)
cmake -P wasm.cmake

# Build full API with tests
cmake -DFULL_BUILD=1 -P wasm.cmake
```

The build script will:
1. Download and install Emscripten SDK (v3.1.51) if needed
2. Configure the WASM build
3. Compile the WASM modules
4. Deploy to `docs/` for GitHub Pages
5. Start a local web server at http://localhost:8832

### Manual Build

```bash
# Configure
cmake -B build-wasm \
    -DCMAKE_TOOLCHAIN_FILE=cmake_modules/toolchain-emscripten.cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -G Ninja

# Build visualization demo
cmake --build build-wasm --target tudat_wasm_web

# Build full Embind API
cmake --build build-wasm --target tudatpy_wasm
```

## Project Structure

```
tudat-wasm/
├── cmake_modules/
│   └── toolchain-emscripten.cmake  # Emscripten toolchain
├── cmake/
│   └── patches/
│       └── cspice-wasm.patch       # SPICE library WASM patch
├── docs/                           # GitHub Pages deployment
├── src/
│   └── tudatpy_wasm/              # WASM bindings source
│       ├── npm/                    # NPM package
│       └── [binding modules]
├── tests/
│   └── wasm/                       # Test suite
│       ├── src/                    # C++ tests
│       ├── web/                    # Browser test UI
│       └── test_wasm_*.js          # Node.js tests
├── wasm.cmake                      # Build orchestration script
└── CMakeLists.txt                  # Main build configuration
```

## NPM Package

The `@tudat/tudatpy-wasm` package provides TypeScript-friendly bindings:

```javascript
import createTudatModule from '@tudat/tudatpy-wasm';

const tudat = await createTudatModule();

// Access astrodynamics functions
const mu = tudat.constants.GRAVITATIONAL_CONSTANT;
```

## Testing

### Node.js Tests

```bash
cd tests/wasm
node test_wasm_api.js
node test_wasm_integration.js
```

### Browser Tests

The build script starts a web server automatically. Alternatively:

```bash
cd tests/wasm/web
python start_server.py
```

Then open http://localhost:8832 in your browser.

## GitHub Pages

The `docs/` directory is deployed to GitHub Pages, providing:
- Interactive orbit visualization demos
- API documentation (TypeDoc)
- Test suite accessible in-browser

## Dependencies

- **Emscripten SDK**: v3.1.51 (auto-downloaded)
- **Tudat**: Fetched via CMake FetchContent
- **Boost**: Headers via Emscripten port
- **Eigen3**: Linear algebra library
- **CSpice**: Planetary ephemerides (with WASM patch)
- **SOFA**: Fundamental astronomy standards
- **CesiumJS**: 3D globe visualization (for demos)

## License

This project follows the same license as Tudat. See [LICENSE](LICENSE) for details.

## Links

- [Tudat Documentation](https://docs.tudat.space/)
- [Tudat GitHub](https://github.com/tudat-team/tudat)
- [TudatPy GitHub](https://github.com/tudat-team/tudatpy)
