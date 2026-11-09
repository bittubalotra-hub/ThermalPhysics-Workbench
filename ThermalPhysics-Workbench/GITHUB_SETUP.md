# Publish ThermalPhysics Workbench on GitHub (Windows PowerShell)

This is a **source-code repository**, not a GitHub Pages website. It becomes visible on your GitHub profile after you push it; the automated numerical tests run under the **Actions** tab.

## 1. Extract and test

Extract the ZIP. Open PowerShell in the folder containing `README.md` and `pyproject.toml`:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest
.\.venv\Scripts\python.exe -m thermalphysics run --config examples/steel_slab_cooling.json --output results/cooling
```

The `results/` and `.venv/` directories are ignored by Git. Avoid uploading local virtual environments.

## 2. Create a new GitHub repository

On GitHub, create **`ThermalPhysics-Workbench`** under your account. Make it public if you want it visible on your profile. Choose **no initial README, license, or .gitignore** because the folder already contains those files.

## 3. Commit and push

In the same project folder, run:

```powershell
git init
git add .
git commit -m "Initial release: ThermalPhysics Workbench"
git branch -M main
git remote add origin https://github.com/bittubalotra-hub/ThermalPhysics-Workbench.git
git push -u origin main
```

**Important:** `git add` and `git commit` must come before `git push`, otherwise `src refspec main does not match any` may occur.

If Git asks for a user identity, configure `git config --global user.name "YOUR_USERNAME"` and `git config --global user.email "YOUR_EMAIL"` then run `git commit` again. If `origin` already exists, check `git remote -v` instead of adding it twice.

## 4. Check and pin

- Look under **Actions** for tests on supported Python versions.
- Open the README and confirm the three figures display.
- In your GitHub profile, select **Customize your pins** to pin the repository.
- Add a link to this repository on your existing portfolio website, and label it an **independent computational research demonstration**.

## 5. Before publishing

Check `docs/DATA_AND_DISCLOSURE.md`. Do not attribute this source code to someone who did not actually write or review it. Do not add confidential employer information, restricted datasets or unapproved technical documents.
