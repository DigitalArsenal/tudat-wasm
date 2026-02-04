import { defineConfig } from 'vite';
import { spawn } from 'child_process';
import path from 'path';

// Plugin to watch C++ files and rebuild WASM
function wasmRebuildPlugin() {
  let isBuilding = false;

  return {
    name: 'wasm-rebuild',
    configureServer(server) {
      const srcDir = path.resolve(__dirname, 'src');
      const includeDir = path.resolve(__dirname, 'include');

      // Watch C++ source files
      server.watcher.add([`${srcDir}/**/*.cpp`, `${srcDir}/**/*.h`, `${includeDir}/**/*.h`]);

      server.watcher.on('change', async (file) => {
        if (isBuilding) return;
        if (!file.endsWith('.cpp') && !file.endsWith('.h') && !file.endsWith('.hpp')) return;

        console.log(`\n[wasm] C++ file changed: ${path.basename(file)}`);
        console.log('[wasm] Rebuilding...');
        isBuilding = true;

        const startTime = Date.now();

        const build = spawn('cmake', ['--build', 'build-wasm', '--target', 'deploy', '-j', '8'], {
          cwd: __dirname,
          stdio: 'inherit'
        });

        build.on('close', (code) => {
          isBuilding = false;
          const duration = ((Date.now() - startTime) / 1000).toFixed(1);

          if (code === 0) {
            console.log(`[wasm] Build complete in ${duration}s`);
            // Trigger full page reload
            server.ws.send({ type: 'full-reload' });
          } else {
            console.error(`[wasm] Build failed with code ${code}`);
          }
        });
      });
    }
  };
}

export default defineConfig({
  root: 'docs',
  plugins: [wasmRebuildPlugin()],
  server: {
    port: 8832,
    open: true,
    headers: {
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'credentialless'
    },
    watch: {
      // Also watch the build output
      ignored: ['!**/docs/**']
    }
  },
  preview: {
    port: 8832,
    headers: {
      'Cross-Origin-Opener-Policy': 'same-origin',
      'Cross-Origin-Embedder-Policy': 'credentialless'
    }
  },
  build: {
    outDir: '../dist',
    emptyOutDir: true
  },
  optimizeDeps: {
    exclude: ['*.wasm']
  },
  assetsInclude: ['**/*.wasm']
});
