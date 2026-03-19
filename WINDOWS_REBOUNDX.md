# Windows REBOUNDx Notes

This document records practical setup and troubleshooting guidance for running `reboundx` on Windows in this project.

## 1. Recommended Environment

- OS: Windows 10/11
- Shell: PowerShell
- Python: 3.10+
- Project venv: `.venv`

## 2. Standard Setup

From project root:

```powershell
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install --upgrade pip setuptools wheel
pip install rebound reboundx
```

Quick verify:

```powershell
.\.venv\Scripts\python.exe -c "import rebound, reboundx; print(rebound.__version__)"
```

## 3. Project-Specific Behavior

- Main prototype enables REBOUNDx GR force in `src/main_simulation.py` by default.
- Temporary local patch/build workspace may use `.tmp_reboundx/` (already gitignored).

## 4. Common Issues and Fixes

## 4.1 `pip install reboundx` fails (build error)

Symptoms:

- No wheel available for your Python version
- Build fails with compiler/toolchain errors

Actions:

1. Ensure Python version has wheel support first (recommended)
2. Update build tools:
   - `pip install --upgrade pip setuptools wheel`
3. Install MSVC C/C++ build tools if source build is required

## 4.2 Import error: DLL load failed / module not found

Symptoms:

- `ImportError` when importing `reboundx`
- `.pyd` fails to resolve dependent DLLs

Actions:

1. Confirm same Python interpreter is used for install and run
2. Reinstall into active venv:
   - `.\.venv\Scripts\python.exe -m pip install --force-reinstall rebound reboundx`
3. Avoid mixing multiple environments (`base`, `.venv`, `.venv_msys`)

## 4.3 ABI / architecture mismatch

Symptoms:

- Package installs but crashes/imports fail

Actions:

1. Ensure Python architecture and wheel architecture match (x64 vs x86)
2. Avoid cross-using binaries from MSYS/MinGW environment in CPython venv

## 4.4 Runtime behavior differs after upgrade

Symptoms:

- Same code, different output after package update

Actions:

1. Pin versions in your environment
2. Run `integrator_check` and compare drift metrics before/after upgrade

## 5. Safe Recovery Procedure

If environment is broken:

```powershell
deactivate
Remove-Item -Recurse -Force .venv
python -m venv .venv
.\.venv\Scripts\activate
python -m pip install --upgrade pip setuptools wheel
pip install -e .
```

Then re-verify:

```powershell
just integrator-check
```

## 6. Notes for Team

- Keep local experimental fixes inside `.tmp_reboundx/`.
- Do not commit generated binaries or local build artifacts.
- If machine-specific workaround is required, document exact steps in this file.

