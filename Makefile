# Makefile -- works on Windows (local) and Linux (Codespaces grader).
# Recipe lines must be TAB-indented (not spaces).
#
# Targets are separate commands:  make setup   then   make pipeline   then   make dashboard

.PHONY: setup pipeline dashboard clean

# Pick interpreter + venv paths per OS. On Windows $(OS) is "Windows_NT".
ifeq ($(OS),Windows_NT)
    VENV_PY := py
    PY      := .venv\Scripts\python.exe
    PIP     := .venv\Scripts\pip.exe
    STREAMLIT := .venv\Scripts\streamlit.exe
else
    VENV_PY := python3
    PY      := .venv/bin/python
    PIP     := .venv/bin/pip
    STREAMLIT := .venv/bin/streamlit
endif

# Create the virtual environment and install dependencies.
setup:
	$(VENV_PY) -m venv .venv
	$(PY) -m pip install --upgrade pip
	$(PIP) install -r requirements.txt

# Build the DB (Part 1) then run all analyses (Parts 2-4).
pipeline:
	$(PY) load_data.py
	$(PY) run_pipeline.py

# Start the dashboard. headless=true skips the first-run email prompt; 0.0.0.0
# lets Codespaces forward the port.
dashboard:
	$(STREAMLIT) run dashboard.py --server.headless true --server.address 0.0.0.0

# Remove generated artifacts (cross-platform via Python).
clean:
	$(PY) -c "import shutil, pathlib; pathlib.Path('cell_count.db').unlink(missing_ok=True); [shutil.rmtree(d, ignore_errors=True) for d in ['outputs', '__pycache__', 'analysis/__pycache__']]"
