
# Compressible Flow Calculator (Streamlit GUI)

A compact, course-ready GUI for **Isentropic flows**, **Normal shocks**, and **Oblique shocks (weak/strong)** under the perfect-gas model. Built with Streamlit.

## 📁 Folder layout
```
compressible_flow_calculator/
├── app.py                  # Streamlit GUI (main entrypoint)
├── requirements.txt        # Minimal dependencies (streamlit, pandas)
├── isentropic_calculator.py        # Optional CLI helper (standalone)
├── normal_shock_calculator.py      # Optional CLI helper (standalone)
└── oblique_shock_calculator.py     # Optional CLI helper (standalone)
```

> Only `app.py` and `requirements.txt` are required for the GUI.  
> The three `*_calculator.py` scripts are optional if you also want terminal CLIs.

---

## ▶️ Run locally

### Prerequisites
- Python **3.9–3.12** recommended
- (Windows) Use **PowerShell** or **CMD**; (macOS/Linux) use your shell

### 1) Create and activate a virtual environment

**Windows (PowerShell):**
```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

**Windows (cmd.exe):**
```cmd
python -m venv .venv
.venv\Scripts\activate.bat
```

**macOS / Linux:**
```bash
python3 -m venv .venv
source .venv/bin/activate
```

> If PowerShell blocks activation, you can temporarily allow it:
> ```powershell
> Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
> .venv\Scripts\Activate.ps1
> ```

### 2) Install dependencies
```bash
pip install -r requirements.txt
```

### 3) Launch the app
```bash
python -m streamlit run app.py
```
Your browser will open at `http://localhost:8501`.

---

## 🖥️ Using the app
- Set **γ** in the sidebar.
- Choose a **module**: *Isentropic*, *Normal Shock*, or *Oblique Shock*.
- Enter the required inputs (units: angles in degrees).
- Results are shown in a clean table with the exact variables you expect for coursework/design checks.

---

## 🧪 Run the optional CLI tools
From the same folder (after activating the venv):
```bash
# Isentropic
python isentropic_calculator.py
# Normal Shock
python normal_shock_calculator.py
# Oblique Shock
python oblique_shock_calculator.py
```

---

## ☁️ Deploy to Streamlit Community Cloud

> **Note:** The app’s code must live in a **public GitHub repository**.

### 1) Create a GitHub repo and push your code
```bash
# Inside the project folder
git init
git add .
git commit -m "Initial commit: compressible flow calculator"
git branch -M main
git remote add origin https://github.com/<your-username>/<your-repo>.git
git push -u origin main
```

### 2) Deploy on Streamlit Community Cloud
1. Go to https://streamlit.io/cloud and sign in (GitHub).
2. Click **New app**.
3. Select your repo, **branch = `main`**, and **file = `app.py`**.
4. Click **Deploy**.

Streamlit Cloud will install packages from `requirements.txt` automatically and start your app.
The first build can take a few minutes.

### (Optional) Pin Python version
Community Cloud generally uses a compatible Python version automatically.  
If you want to pin it, add a `runtime.txt` file:
```
python-3.11
```

Place it in the project root and push to GitHub.

### (Optional) App settings
You can add an optional `.streamlit/config.toml` to tweak the theme:
```toml
# .streamlit/config.toml
[theme]
base="light"
```
This file is not required.

---

## ❗ Troubleshooting

- **`streamlit: command not found`**  
  Streamlit isn’t installed in your venv. Run:
  ```bash
  pip install -r requirements.txt
  python -m streamlit run app.py
  ```

- **PowerShell says `source` not recognized**  
  Use the Windows activation command:
  ```powershell
  .venv\Scripts\Activate.ps1
  ```

- **Port already in use**  
  Run on a different port:
  ```bash
  python -m streamlit run app.py --server.port 8502
  ```

- **White screen / ModuleNotFoundError on Cloud**  
  Ensure `requirements.txt` is committed and in the repository root.
  Reboot the app from Streamlit Cloud.

---

## 📜 License
This project is for educational use. 
