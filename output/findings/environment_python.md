# Python environment

Recorded 2026-09-15 by `analyses/00_setup/02_record_environment.py`.

| | |
|---|---|
| platform | `macOS-27.0-arm64-arm-64bit` |
| machine | `arm64` |
| python | 3.12.2 |
| pymc | 6.1.0 |
| pytensor | 3.1.3 |
| arviz | 1.2.0 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| pandas | 3.0.3 |
| matplotlib | 3.10.9 |
| networkx | 3.6.1 |
| geopandas | 1.1.3 |
| shapely | 2.1.2 |
| pyproj | 3.7.2 |
| cartopy | 0.25.0 |
| openpyxl | 3.1.5 |
| xlrd | 2.0.2 |
| h5py | 3.16.0 |
| h5netcdf | 1.8.1 |
| pytest | 9.0.3 |

## Why this is recorded

`pyproject.toml` pins lower bounds only, so a rebuild does not reproduce the original environment. Any discrepancy against a previously committed number has to be attributed to a version change or to a defect before it is interpreted. Measured instance: the hierarchical convergence model's divergence count at the library default moved from at-or-below 5 to 7 across this rebuild (`output/findings/convergence_model_geometry.md`), while the Balding-Nichols fit reproduced its committed values exactly.

<details><summary>Full pip freeze</summary>

```
anyio==4.13.0
appnope==0.1.4
argon2-cffi==25.1.0
argon2-cffi-bindings==25.1.0
arrow==1.4.0
arviz==1.2.0
arviz-base==1.2.0
arviz-plots==1.2.0
arviz-stats==1.2.0
asttokens==3.0.1
async-lru==2.3.0
attrs==26.1.0
babel==2.18.0
beautifulsoup4==4.15.0
bleach==6.4.0
cachetools==6.2.6
Cartopy==0.25.0
certifi==2026.5.20
cffi==2.0.0
charset-normalizer==3.4.7
cloudpickle==3.1.2
comm==0.2.3
cons==0.4.7
contourpy==1.3.3
coverage==7.14.1
cycler==0.12.1
debugpy==1.8.21
decorator==5.3.1
defusedxml==0.7.1
et_xmlfile==2.0.0
etuples==0.3.10
executing==2.2.1
fastjsonschema==2.21.2
filelock==3.29.7
fonttools==4.63.0
fqdn==1.5.1
geopandas==1.1.3
h11==0.16.0
h5netcdf==1.8.1
h5py==3.16.0
httpcore==1.0.9
httpx==0.28.1
idna==3.18
iniconfig==2.3.0
ipykernel==7.2.0
ipython==9.14.1
ipython_pygments_lexers==1.1.1
ipywidgets==8.1.8
isoduration==20.11.0
jedi==0.20.0
Jinja2==3.1.6
json5==0.14.0
jsonpointer==3.1.1
jsonschema==4.26.0
jsonschema-specifications==2025.9.1
jupyter==1.1.1
jupyter-console==6.6.3
jupyter-events==0.12.1
jupyter-lsp==2.3.1
jupyter_client==8.9.1
jupyter_core==5.9.1
jupyter_server==2.19.0
jupyter_server_terminals==0.5.4
jupyterlab==4.5.8
jupyterlab_pygments==0.3.0
jupyterlab_server==2.28.0
jupyterlab_widgets==3.0.16
kiwisolver==1.5.0
lark==1.3.1
lazy-loader==0.5
llvmlite==0.47.0
logical-unification==0.4.7
markdown-it-py==4.2.0
MarkupSafe==3.0.3
matplotlib==3.10.9
matplotlib-inline==0.2.2
mdurl==0.1.2
miniKanren==1.0.5
mistune==3.2.1
-e git+https://github.com/clipo/mls-emergence.git@e1134c661583293979b9741b7c61f0675222ee62#egg=mls_emergence
-e git+https://github.com/rdinapoli/monument-mls@299dac76812e103d2aa0c39abc04f84ea0e7ec58#egg=monument_mls
mpmath==1.3.0
multipledispatch==1.0.0
nbclient==0.11.0
nbconvert==7.17.1
nbformat==5.10.4
nest-asyncio==1.6.0
networkx==3.6.1
notebook==7.5.7
notebook_shim==0.2.4
numba==0.65.1
numpy==2.4.6
openpyxl==3.1.5
packaging==26.2
pandas==3.0.3
pandocfilters==1.5.1
parso==0.8.7
pexpect==4.9.0
pillow==12.2.0
platformdirs==4.10.0
pluggy==1.6.0
prometheus_client==0.25.0
prompt_toolkit==3.0.52
psutil==7.2.2
ptyprocess==0.7.0
pure_eval==0.2.3
pycparser==3.0
Pygments==2.20.0
pymc==6.1.0
pyogrio==0.12.1
pyparsing==3.3.2
pypdf==6.13.1
pyproj==3.7.2
pyshp==3.0.12
pytensor==3.1.3
pytest==9.0.3
pytest-cov==7.1.0
python-dateutil==2.9.0.post0
python-json-logger==4.1.0
PyYAML==6.0.3
pyzmq==27.1.0
referencing==0.37.0
requests==2.34.2
rfc3339-validator==0.1.4
rfc3986-validator==0.1.1
rfc3987-syntax==1.1.0
rich==15.0.0
rpds-py==2026.5.1
scipy==1.17.1
Send2Trash==2.1.0
setuptools==82.0.1
shapely==2.1.2
six==1.17.0
soupsieve==2.8.4
stack-data==0.6.3
sympy==1.14.0
terminado==0.18.1
threadpoolctl==3.6.0
tinycss2==1.5.1
toolz==1.1.0
tornado==6.5.7
traitlets==5.15.1
typing_extensions==4.15.0
tzdata==2026.2
uri-template==1.3.0
urllib3==2.7.0
wcwidth==0.8.1
webcolors==25.10.0
webencodings==0.5.1
websocket-client==1.9.0
widgetsnbextension==4.0.15
xarray==2026.7.0
xarray-einstats==0.10.0
xlrd==2.0.2
```

</details>
