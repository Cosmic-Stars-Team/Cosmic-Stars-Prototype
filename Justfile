set shell := ["pwsh", "-NoLogo", "-NoProfile", "-Command"]

venv := ".venv"
python := ".venv\\Scripts\\python.exe"
pycache := "__pycache__"
py_path := "src"

default:
    just integrator-check

venv:
    if (-not (Test-Path "{{python}}")) { if (Get-Command py -ErrorAction SilentlyContinue) { py -3 -m venv {{venv}} } elseif (Get-Command python -ErrorAction SilentlyContinue) { python -m venv {{venv}} } else { throw "Python not found. Please install Python first." } }

install: venv
    {{python}} -m pip install --upgrade pip setuptools wheel
    {{python}} -m pip install -e .

integrator-check args="": install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; if ("{{args}}" -eq "") { {{python}} -m integrator_check --years 10 --steps 20000 --sample-stride 200 --sample-bodies 'Mercury,Earth' --output-csv 'data/gen/integrator_samples.csv' } else { {{python}} -m integrator_check {{args}} }

simulate: install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; {{python}} -m main_simulation

clean-pycache:
    if (Test-Path "__pycache__") { Remove-Item -Recurse -Force "__pycache__" }
