import { terser } from "rollup-plugin-terser";

export default [
    {
      input: "src/index.js",
      output: {
        file: "dist/miniape.js",
        format: "es",
        sourcemap: true,
      },
      external: [
        './script/*',
        './test/*',
      ]
    },
    {
      input: "src/index.js",
      output: {
        file: "dist/miniape.min.js",
        format: "es",
        sourcemap: true
      },
      plugins: [terser()]
    },
    {
      input: "src/index.js",
      output: {
        file: "dist/miniape.cjs",
        format: "cjs",
        sourcemap: true
      }
    }
  ];
  