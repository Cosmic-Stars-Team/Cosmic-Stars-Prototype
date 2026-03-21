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

# 启用 REBOUNDx 的水星近日点进动测试（固定 2000 年，开箱即用）
mercury-perihelion: install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; {{python}} -m integrator_check --integrator mercurius --use-reboundx --report-mercury-perihelion --years 2000 --steps 400000 --sample-bodies 'Mercury' --sample-stride 2000 --output-csv 'data/gen/mercury_reboundx_2000y.csv'

# 开/关 REBOUNDx 对照测试（固定 2000 年，输出两个 CSV）
mercury-perihelion-compare: install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; {{python}} -m integrator_check --integrator mercurius --compare-reboundx --report-mercury-perihelion --years 2000 --steps 400000 --sample-bodies 'Mercury' --sample-stride 2000 --output-csv 'data/gen/mercury_perihelion_2000y.csv'

# 启用 REBOUNDx 的水星近日点进动测试（可自定义参数）
mercury-perihelion-custom args="": install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; if ("{{args}}" -eq "") { {{python}} -m integrator_check --integrator mercurius --use-reboundx --report-mercury-perihelion --years 2000 --steps 400000 --sample-bodies 'Mercury' --sample-stride 2000 --output-csv 'data/gen/mercury_reboundx_custom.csv' } else { {{python}} -m integrator_check --integrator mercurius --use-reboundx --report-mercury-perihelion {{args}} }

# 更短的别名（更顺手）
mp:
    just mercury-perihelion

mpc:
    just mercury-perihelion-compare

mpx args="":
    just mercury-perihelion-custom args="{{args}}"

# Water Mercury perihelion quick help (prints to terminal)
mp-help:
    @Write-Host "Mercury perihelion tests (REBOUNDx GR):"
    @Write-Host ""
    @Write-Host "  just mp"
    @Write-Host "    - Run 2000y with REBOUNDx enabled, prints drift + perihelion rate."
    @Write-Host "    - Writes: data/gen/mercury_reboundx_2000y.csv"
    @Write-Host ""
    @Write-Host "  just mpc"
    @Write-Host "    - Compare WITHOUT vs WITH REBOUNDx (same params), prints both rates + delta."
    @Write-Host "    - Writes: data/gen/mercury_perihelion_2000y_no_reboundx.csv"
    @Write-Host "             data/gen/mercury_perihelion_2000y_with_reboundx.csv"
    @Write-Host ""
    @Write-Host "Customize years/steps (dt = years/steps):"
    @Write-Host "  just mpx args=\"--years 5000 --steps 1000000 --sample-stride 5000 --output-csv data/gen/mercury_5000y.csv\""
    @Write-Host ""
    @Write-Host "Tip: smaller dt (more steps) is slower but usually more stable/accurate."

simulate: install
    New-Item -ItemType Directory -Path data\\gen -Force | Out-Null
    $env:PYTHONPYCACHEPREFIX='{{pycache}}'; $env:PYTHONPATH='{{py_path}}'; {{python}} -m main_simulation

clean-pycache:
    if (Test-Path "__pycache__") { Remove-Item -Recurse -Force "__pycache__" }
