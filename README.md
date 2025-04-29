# IVV ETF Performance & Risk Modeling: Metrics, Efficient Frontier & Monte Carlo VaR (2010–2023)

![R](https://img.shields.io/badge/R-%3E%3D4.0-blue) ![License: MIT](https://img.shields.io/badge/license-MIT-green)

**FE535 Project 1**  
IVV ETF performance analysis, mean-variance efficient frontier construction & stochastic risk modeling (GBM & Monte Carlo) over Jan 2010–Sep 2023.

---

## Table of Contents

1. [Files](#files)  
2. [Data](#data)  
3. [Usage](#usage)  
4. [Results & Report](#results--report)  
5. [Project Structure](#project-structure)  
6. [Future Work](#future-work)  
7. [License](#license)  

---

## Files

- `data/raw/etf_prices.csv`  — Daily adjusted closing prices for 30 iShares ETFs  
- `r/main.R`                 — R script (quantmod + ggplot2) that loads data, computes metrics, plots, and prints tables  
- `Final Project.pdf`        — Written report with narrative & figures  
- `.gitignore`               — Excludes R artifacts, raw data & results  

---

## Data

1. **Source:** iShares / Yahoo! Finance  
2. **Path:** `data/raw/etf_prices.csv`  
3. **Format:** CSV with `Date` and 30 ETF tickers

---

## Usage

```bash
# Install (in R console):
install.packages(c("quantmod","ggplot2","PerformanceAnalytics"))

# Run the analysis:
Rscript r/main.R
```

## Project Structure

fe535-etf-risk-modeling/
├── data/
│   └── raw/etf_prices.csv
├── r/main.R
├── Final Project.pdf
├── .gitignore
└── README.md

---

## Future Work

- Python port of the analysis
- Unit tests for metric functions
- CI/CD via GitHub Actions
(Python version coming soon!)

---

## License
MIT © Mrinal Gupta
