# TF-NRD v1.0 Webpage Deployment & Requirements

This directory contains the production-ready web application suite for the **TF-NRD (Transcription Factors Non-Redundant Dataset)** database, hosted at:

🔗 **[http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html](http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/tfnrd.html)**

---

## ⚡ Key Architecture & Zero-Dependency Design

The TF-NRD webpages are engineered as **pure static, zero-dependency client-side web applications**:

- **No Backend Server Required**: No PHP, Python CGI, Node.js, or Ruby runtime needed on the production webserver.
- **No Database Server Required**: No MySQL, PostgreSQL, or MongoDB setup required. Datasets are pre-compiled and embedded directly as client-side JSON arrays within each page for instantaneous (0 ms latency) querying.
- **No Build Tools or Bundlers**: No Webpack, Vite, npm, or Babel compilation steps required. Everything runs natively in standard web browsers.
- **Offline / Local File Protocol Compatible**: Can be opened directly in any web browser via `file:///` without triggering CORS errors.

---

## 🖥️ Web Server Hosting Requirements

Any standard HTTP static file server can host this directory:

- **Supported Web Servers**:
  - Apache HTTP Server (`httpd`)
  - Nginx
  - Caddy
  - LiteSpeed
  - Built-in Python HTTP server (for local testing)

### Recommended Server Configurations

#### 1. Directory Index (Default Landing Page)
To load `tfnrd.html` automatically when users visit `http://www.csb.iitkgp.ac.in/databases/TFNRDv1.0/`:

- **Apache (`.htaccess` or virtual host configuration)**:
  ```apache
  DirectoryIndex tfnrd.html index.html
  ```

- **Nginx (`nginx.conf`)**:
  ```nginx
  location /databases/TFNRDv1.0/ {
      index tfnrd.html index.html;
  }
  ```

#### 2. Compression & Caching (Recommended for High Performance)
Because the HTML files contain embedded structured dataset records, enabling `gzip` or `brotli` compression on the web server reduces transmission size by ~80%:

- **Apache `mod_deflate`**:
  ```apache
  <IfModule mod_deflate.c>
      AddOutputFilterByType DEFLATE text/html text/css application/javascript application/json
  </IfModule>
  ```

- **Nginx**:
  ```nginx
  gzip on;
  gzip_types text/html text/css application/javascript application/json;
  ```

---

## 🌐 Client Browser Requirements

- **Supported Modern Browsers**:
  - Google Chrome / Chromium 80+
  - Mozilla Firefox 75+
  - Apple Safari 13.1+
  - Microsoft Edge 80+
  - Mobile browsers (Safari iOS, Chrome Android)
- **JavaScript**: Enabled (required for client-side search, category filtering, column sorting, pagination, and data export).
- **Fonts**: Automatically loads Google Fonts (`Inter`, `Plus Jakarta Sans`, `JetBrains Mono`) with standard system fallbacks (`-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif, monospace`) for offline environments.

---

## 📁 File Manifest

```text
TFNRDv1.0/
├── README.md               # Server requirements and deployment guide
├── tfnrd.html              # Main Landing / Home page with dataset statistics
├── TFNRDv1.0.html          # Table S1.A: Biological assemblies & complex interfaces (377 entries)
├── TFNRDv1.0_PNA.html      # Table S1.B: Unique TF-NA interfaces dataset (509 entries)
├── TFNRDv1.0_PP.html       # Table S1.C: Protein-Protein contact interfaces (582 entries)
└── css/
    └── style.css           # Unified modern responsive stylesheet
```

---

## 🧪 Local Testing & Preview

To preview the webpages locally before uploading to the server:

```bash
# Option 1: Direct browser launch
xdg-open tfnrd.html   # Linux
open tfnrd.html       # macOS

# Option 2: Lightweight Python HTTP server
cd /home/labuser/Projects/PhD_projects/TF-NRD/TFNRDv1.0
python3 -m http.server 8000
# Then visit http://localhost:8000/tfnrd.html in your browser
```

---

## 🔄 Rebuilding / Updating Webpages

If the source supplementary data in `../Supplementary/Supplementary.xlsx` is modified or updated, the webpages can be regenerated in seconds:

```bash
# From repository root:
python script/generate_webpages.py
```

*Requirements for regeneration: Python 3.9+ with `pandas` and `openpyxl` (see root `requirements.txt`).*
