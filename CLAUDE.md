# UVQgit — MATLAB Scripts Repo

## What this is
Brent's MATLAB scripts for photonics/waveguide analysis at Uviquity. Scripts cover SHG efficiency modeling, overhead loss analysis, image processing, spectrometer data, COMSOL sweep processing, and more.

## Setup status
- **VS Code extension**: MathWorks MATLAB extension (v1.3.12) installed. `matlab.installPath` set to `C:\Program Files\MATLAB\R2018b` and `matlabConnectionTiming` set to `"never"` (R2018b doesn't support the language server — syntax highlighting works, run button does not).
- **Git**: Initialized. Two commits made (initial + .gitignore cleanup). Remote is `https://github.com/uviquityBF/UVQgit_matlab`.
- **GitHub push**: PENDING — was blocked by plane WiFi during setup. Token is embedded in the remote URL. Once on solid WiFi, just run:
  ```
  git push -u origin main
  ```
  After a successful push, clean the token out of the remote URL with:
  ```
  git remote set-url origin https://github.com/uviquityBF/UVQgit_matlab
  ```

## Running scripts from VS Code terminal
R2018b doesn't support `-batch`. Use this form instead:
```
matlab -nosplash -nodesktop -r "run('scriptname.m'); exit"
```
Or without auto-exit (keeps figures open):
```
matlab -nosplash -r "run('scriptname.m')"
```

## .gitignore covers
- `*.asv` (MATLAB autosave)
- `*.mex*` / `*.mexw64` etc. (compiled MEX)
- `desktop.ini`, `Thumbs.db` (Windows metadata — these were corrupting the repo when Google Drive scattered them into .git/)
- `temp.m`

## Google Drive note
This folder lives on Google Drive (`G:\My Drive\...`). Avoid running git commands while Drive is actively syncing — it can create lock file conflicts. The previous .git was destroyed by Drive scattering desktop.ini files into .git/objects/ and .git/refs/.

## User's goal
Brent wants Claude's help writing new and improved MATLAB scripts — refactoring for reusability, vectorizing loops, better figure formatting, new analysis scripts from scratch.
