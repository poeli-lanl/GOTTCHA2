# Coverage browser assets

These pinned upstream distributions are embedded in each coverage report. No CDN,
network connection, or application server is required to open the generated HTML.

| Library | Version | License |
| --- | --- | --- |
| D3 | 7.9.0 | ISC |
| Vue | 3.5.13 | MIT |
| PrimeVue | 4.3.2 | MIT |
| PrimeUIX themes (Aura) | 1.0.0 | MIT |
| Bootstrap | 5.3.3 | MIT |
| PrimeIcons | 7.0.0 | MIT |

`sources.json` records the npm registry tarball URLs and verified integrity hashes.
The upstream license texts are in `licenses/`; retain them when updating these files.
Only the runtime JavaScript, CSS, and WOFF2 font needed by the renderer are included.
