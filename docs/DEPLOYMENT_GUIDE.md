# Deployment Guide: Rare Earth Ligand Database Dashboard

This guide covers deploying the Streamlit dashboard on a new host (server, VM, or local machine).

---

## Prerequisites

### System Requirements
- **OS**: Linux, macOS, or Windows
- **Python**: 3.8 or higher (3.10+ recommended)
- **RAM**: 2GB minimum, 4GB recommended
- **Disk**: ~500MB for application + dependencies
- **Network**: Port 8501 (default) or custom port for web access

---

## Quick Start (5 Minutes)

### Option 1: Clone from Git Repository

```bash
# 1. Clone the repository
git clone <your-git-repo-url>
cd Database

# 2. Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Launch dashboard
streamlit run app.py
```

Dashboard will be available at: **http://localhost:8501**

---

## Detailed Setup Instructions

### Step 1: Copy Required Files

**Essential files and directories:**

```
Database/
├── app.py                              # Main dashboard application (required)
├── requirements.txt                    # Python dependencies (required)
├── data/
│   └── sx_ligands_merged_v14_flat.jsonl     # Main database (~1.1MB, required)
├── rel_quest/                             # Python package (~217KB, required)
│   ├── __init__.py
│   ├── core.py
│   ├── smarts_patterns.py
│   ├── solvent_mapper.py
│   └── utils.py
└── vocab/                              # Vocabulary files (~24KB, required)
    ├── acids.json
    ├── elements.json
    ├── mechanisms.json
    └── solvents.json
```

**Optional files (for testing and documentation):**
- `tests/` - Automated tests
- `docs/` and the other `*.md` documentation files

**NOT required for deployment:**
- `scripts/` - Build/maintenance scripts
- `*.pptx`, `*.txt`, backup files

**Minimum deployment size**: ~2MB (excluding dependencies)

---

### Step 2: Install Python Dependencies

**Option A: Using pip (standard)**

```bash
pip install -r requirements.txt
```

**Dependencies installed:**
- `streamlit>=1.30.0` - Web framework
- `pandas>=2.0.0` - Data manipulation
- `plotly>=5.18.0` - Interactive plots
- `scipy>=1.11.0` - Statistical functions
- `numpy>=1.24.0` - Numerical operations
- `rdkit>=2023.9.1` - Molecular structure rendering (optional)
- `pillow>=10.0.0` - Image processing (optional)

**Option B: Using conda (recommended for RDKit)**

```bash
# Create conda environment
conda create -n ligand-db python=3.10

# Activate environment
conda activate ligand-db

# Install RDKit from conda-forge (recommended method)
conda install -c conda-forge rdkit

# Install other dependencies
pip install streamlit>=1.30.0 pandas>=2.0.0 plotly>=5.18.0 scipy>=1.11.0 numpy>=1.24.0 pillow>=10.0.0
```

**Note on RDKit**: If RDKit installation fails, the dashboard will still work but won't display molecular structures (will show SMILES text instead).

---

### Step 3: Configure Application (Optional)

#### Change Database Path

If you want to use a different database file location, edit `app.py` line 136:

```python
# app.py line 136
db_path = 'data/sx_ligands_merged_v14_flat.jsonl'  # Change this path if needed
```

#### Change Default Port

```bash
# Run on custom port
streamlit run app.py --server.port 8080
```

#### Enable External Access

```bash
# Allow access from other machines
streamlit run app.py --server.address 0.0.0.0
```

---

## Deployment Scenarios

### Scenario 1: Local Development

**Use case**: Testing on your laptop/desktop

```bash
cd Database
source venv/bin/activate
streamlit run app.py
```

Access at: `http://localhost:8501`

---

### Scenario 2: Internal Server (LAN)

**Use case**: Lab/organization internal server

```bash
# On the server
streamlit run app.py --server.address 0.0.0.0 --server.port 8501

# Access from other machines on the network
# http://<server-ip>:8501
```

**Firewall configuration** (if needed):
```bash
# Ubuntu/Debian
sudo ufw allow 8501/tcp

# CentOS/RHEL
sudo firewall-cmd --permanent --add-port=8501/tcp
sudo firewall-cmd --reload
```

---

### Scenario 3: Production Server (Internet-Facing)

**Use case**: Public or authenticated web service

#### Option A: Using systemd (Linux)

Create service file: `/etc/systemd/system/ligand-dashboard.service`

```ini
[Unit]
Description=Ligand Database Dashboard
After=network.target

[Service]
Type=simple
User=ligand-user
WorkingDirectory=/opt/ligand-db
Environment="PATH=/opt/ligand-db/venv/bin"
ExecStart=/opt/ligand-db/venv/bin/streamlit run app.py --server.port 8501 --server.address 0.0.0.0
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

**Start service:**
```bash
sudo systemctl daemon-reload
sudo systemctl enable ligand-dashboard
sudo systemctl start ligand-dashboard
sudo systemctl status ligand-dashboard
```

#### Option B: Using Docker

Create `Dockerfile`:

```dockerfile
FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy application files
COPY requirements.txt .
# Copy the whole checkout: LigandDatabase resolves data/ and vocab/ relative
# to the package, and pyproject.toml is needed for `pip install -e .`.
# The repo's own Dockerfile does exactly this.
COPY . .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Expose port
EXPOSE 8501

# Run application
CMD ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
```

**Build and run:**
```bash
# Build image
docker build -t ligand-dashboard .

# Run container
docker run -d -p 8501:8501 --name ligand-db ligand-dashboard

# View logs
docker logs -f ligand-db
```

#### Option C: Using Nginx Reverse Proxy

**Setup Nginx** (for HTTPS and domain mapping):

```nginx
# /etc/nginx/sites-available/ligand-dashboard
server {
    listen 80;
    server_name ligands.example.com;

    location / {
        proxy_pass http://localhost:8501;
        proxy_http_version 1.1;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header Host $host;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_read_timeout 86400;
    }
}
```

**Enable site:**
```bash
sudo ln -s /etc/nginx/sites-available/ligand-dashboard /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl reload nginx
```

**Add HTTPS with Let's Encrypt:**
```bash
sudo certbot --nginx -d ligands.example.com
```

---

### Scenario 4: Cloud Hosting

#### Streamlit Community Cloud (Easiest)

This is how the dashboard is shared with external users. Everything the app
reads at runtime (`data/sx_ligands_merged_v14_flat.jsonl`, `vocab/*.json`) is
tracked in git and resolved relative to the `rel_quest` package, so the clone
is self-sufficient — no data upload step.

1. Push the branch to GitHub. For a private repository, authorize the
   Streamlit GitHub app for it when prompted.
2. Go to https://share.streamlit.io → **New app**.
3. Choose the repository and its default branch, main file **`app.py`**.
4. **Advanced settings → Python 3.11** — matches the CI matrix and guarantees
   prebuilt rdkit/scipy wheels.
5. **Settings → Secrets**: paste the contents of
   [`.streamlit/secrets.toml.example`](../.streamlit/secrets.toml.example) with
   a real key. See "AI assistant on Streamlit Cloud" below.
6. Deploy. Python dependencies install from `requirements.txt`; `packages.txt`
   adds the X11 shared libraries (`libxrender1`, `libxext6`, `libx11-6`) that
   RDKit's drawing code links against. Without them RDKit raises `ImportError`
   on a slim image and molecular structures silently stop rendering.

**Pros**: Free, automatic HTTPS, redeploys on push
**Cons**: Public by default, resource limits, no local Ollama

##### AI assistant on Streamlit Cloud

The assistant cannot reach a local Ollama server from Cloud, so it offers two
cloud modes instead. The local-Ollama mode auto-hides: it is shown only when an
Ollama server actually answers (or when `RELQUEST_ENABLE_OLLAMA` is set), so
cloud visitors never see a dead option.

| Mode | Key | Models |
|---|---|---|
| 🎁 Free models | the app owner's `OPENROUTER_FREE_API_KEY`, never shown to visitors | zero-cost OpenRouter IDs only — a `:free` suffix or the `openrouter/free` router; the app refuses anything else, so the shared key cannot be spent on paid usage |
| 🔑 Your own API key | entered by the visitor, session-only | any OpenAI-compatible provider (Anthropic, OpenRouter, OpenAI, DeepSeek, Groq) |

The free-mode default is `openrouter/free`, OpenRouter's "Free Models Router",
which dispatches to whichever free models are currently up — unlike a pinned
model ID, it cannot go stale when a provider retires one.

Claude models work in the your-own-key mode through Anthropic's OpenAI-SDK
compatibility layer: base URL `https://api.anthropic.com/v1/` (keep the
trailing slash) with a model such as `claude-opus-5`. Streaming is fully
supported there; prompt caching and strict tool schemas are not, neither of
which this app uses.

Set the shared key in **Settings → Secrets**:

```toml
OPENROUTER_FREE_API_KEY = "sk-or-v1-..."
```

Omit it and the free mode simply is not offered — visitors then supply their
own key. Config is read by `rel_quest.llm_backends.get_config`, which checks
`st.secrets` first and falls back to environment variables, so the same names
work locally, in Docker, and on Cloud.

**Free-tier quota**: OpenRouter allows roughly 50 requests/day per account
(1000/day once the account holds ≥$10 of credit), and that quota is shared by
*every* visitor using free mode. The UI says so up front, and a rate-limit error
points the user at ⚙️ Settings to add their own key.

#### AWS EC2

```bash
# Launch EC2 instance (Ubuntu 22.04, t2.small or larger)
# SSH into instance
ssh -i your-key.pem ubuntu@<ec2-ip>

# Install dependencies
sudo apt update
sudo apt install python3-pip python3-venv

# Clone repo and setup
git clone <repo-url>
cd Database
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Run with systemd (see Scenario 3)
```

#### Google Cloud Run (Serverless)

```bash
# Build and deploy
gcloud builds submit --tag gcr.io/PROJECT-ID/ligand-dashboard
gcloud run deploy --image gcr.io/PROJECT-ID/ligand-dashboard --platform managed
```

#### Heroku

Create `Procfile`:
```
web: streamlit run app.py --server.port=$PORT --server.address=0.0.0.0
```

Deploy:
```bash
heroku create ligand-dashboard
git push heroku master
```

---

## Configuration Options

### Streamlit Configuration File

Create `.streamlit/config.toml` for custom settings:

```toml
[server]
port = 8501
address = "0.0.0.0"
maxUploadSize = 200
enableXsrfProtection = true
enableCORS = false

[browser]
serverAddress = "ligands.example.com"
serverPort = 443
gatherUsageStats = false

[theme]
primaryColor = "#0068c9"
backgroundColor = "#FFFFFF"
secondaryBackgroundColor = "#F0F2F6"
textColor = "#262730"
font = "sans serif"

[logger]
level = "info"
```

### Environment Variables

```bash
# Set custom port
export STREAMLIT_SERVER_PORT=8080

# Set custom database path
export LIGAND_DB_PATH=/path/to/custom/database.jsonl

# Enable debug mode
export STREAMLIT_LOGGER_LEVEL=debug
```

---

## Verification & Testing

### 1. Check Installation

```bash
# Test Python import
python3 -c "from rel_quest import LigandDatabase; print('✓ rel_quest package works')"

# Test database loading
python3 -c "from rel_quest import LigandDatabase; db = LigandDatabase('data/sx_ligands_merged_v14_flat.jsonl'); print(f'✓ Loaded {len(db.get_all_ligands())} ligands')"

# Test RDKit (optional)
python3 -c "from rdkit import Chem; print('✓ RDKit works')"
```

**Expected output:**
```
✓ rel_quest package works
✓ Loaded 170 ligands
✓ RDKit works
```

### 2. Run Test Suite

```bash
python3 -m pytest -q
```

**Expected output:**
```
All tests passed! ✓
```

### 3. Launch Dashboard

```bash
streamlit run app.py
```

**Check for errors in terminal** - should see:
```
You can now view your Streamlit app in your browser.

Local URL: http://localhost:8501
Network URL: http://192.168.x.x:8501
```

### 4. Access Dashboard

Open browser to `http://localhost:8501` and verify:
- ✅ Page loads without errors
- ✅ Sidebar shows filters
- ✅ Results count badge appears (🟢 170 ligands found)
- ✅ Database Overview plots render
- ✅ Ligand Browser table loads
- ✅ "Jump to Results" button works

---

## Troubleshooting

### Issue 1: ModuleNotFoundError: No module named 'rel_quest'

**Cause**: rel_quest package not in Python path

**Solution**:
```bash
# Check if rel_quest/ directory exists
ls -la rel_quest/

# Make sure you're running from the Database/ directory
cd /path/to/Database
streamlit run app.py
```

### Issue 2: FileNotFoundError: data/sx_ligands_merged_v14_flat.jsonl

**Cause**: Database file missing or wrong path

**Solution**:
```bash
# Check if file exists
ls -lh data/sx_ligands_merged_v14_flat.jsonl

# If missing, you may need to build it or copy from source
# For deployment, copy the file from the original installation
```

### Issue 3: RDKit Import Error

**Cause**: RDKit not installed properly

**Solution**:
```bash
# Option 1: Install via conda (recommended)
conda install -c conda-forge rdkit

# Option 2: Skip RDKit (dashboard still works, just no molecular structures)
# Edit requirements.txt and comment out:
# rdkit>=2023.9.1
```

**Note**: If RDKit fails to import, the dashboard automatically falls back to showing SMILES text instead of molecular structures.

### Issue 4: Port Already in Use

**Cause**: Port 8501 already taken

**Solution**:
```bash
# Option 1: Kill existing process
lsof -ti:8501 | xargs kill -9

# Option 2: Use different port
streamlit run app.py --server.port 8080
```

### Issue 5: "External URL" Not Working

**Cause**: Firewall blocking port or wrong address binding

**Solution**:
```bash
# Check if port is listening
netstat -tuln | grep 8501

# Make sure you're binding to 0.0.0.0
streamlit run app.py --server.address 0.0.0.0

# Check firewall (Ubuntu/Debian)
sudo ufw status
sudo ufw allow 8501/tcp
```

### Issue 6: Slow Performance

**Cause**: Large database or insufficient resources

**Solution**:
```bash
# Enable caching (already enabled in app.py)
# Increase memory if on cloud platform
# Use smaller database subset for testing
```

### Issue 7: CSS/HTML Not Rendering (Button Shows Raw HTML)

**Cause**: Streamlit version issue or HTML sanitizer blocking

**Solution**:
```bash
# Upgrade Streamlit
pip install --upgrade streamlit>=1.30.0

# Clear Streamlit cache
rm -rf ~/.streamlit/cache
```

---

## Performance Optimization

### For Production Deployment

1. **Enable Caching** (already enabled in `app.py`):
   ```python
   @st.cache_resource
   def load_database():
       # ... cached database loading
   ```

2. **Use CDN for Static Assets** (if self-hosting):
   - Serve plots as static images for faster loading
   - Cache Plotly.js library

3. **Limit Concurrent Users**:
   ```toml
   # .streamlit/config.toml
   [server]
   maxConcurrentUsers = 10
   ```

4. **Optimize Database**:
   - Currently at 170 ligands (~1.5 MB) - already optimized
   - If database grows significantly, consider SQLite or PostgreSQL

5. **Use Load Balancer** (for high traffic):
   - Deploy multiple instances behind Nginx/HAProxy
   - Session affinity for Streamlit's websockets

---

## Security Considerations

### Authentication

**Option 1: Basic HTTP Auth (Nginx)**
```nginx
location / {
    auth_basic "Ligand Database";
    auth_basic_user_file /etc/nginx/.htpasswd;
    proxy_pass http://localhost:8501;
}
```

**Option 2: Streamlit Cloud Authentication**
- Built-in Google/GitHub OAuth
- Free tier includes authentication

**Option 3: Custom Authentication**
- Implement in `app.py` using `st.session_state`
- Or use Streamlit-Authenticator package

### HTTPS

**Always use HTTPS for production**:
- Use Nginx reverse proxy with SSL certificate
- Use Let's Encrypt for free SSL certificates
- Or use cloud platform's built-in HTTPS (Streamlit Cloud, AWS, etc.)

### API Keys / Secrets

Secrets live in `.streamlit/secrets.toml` locally (gitignored) or in the
Streamlit Cloud **Settings → Secrets** box. Copy
[`.streamlit/secrets.toml.example`](../.streamlit/secrets.toml.example) and
fill it in — never commit the filled-in file.

```toml
OPENROUTER_FREE_API_KEY = "sk-or-v1-..."   # powers the assistant's free mode
```

Read them through the helper rather than `st.secrets` directly, so the same
name also works as an environment variable (Docker, local runs, tests):

```python
from rel_quest.llm_backends import get_config
api_key = get_config("OPENROUTER_FREE_API_KEY")
```

---

## Backup & Maintenance

### Regular Backups

**Database file**:
```bash
# Backup database
cp data/sx_ligands_merged_v14_flat.jsonl data/sx_ligands_merged_v6_backup_$(date +%Y%m%d).jsonl

# Automated backup (cron)
0 2 * * * cd /opt/ligand-db && cp data/sx_ligands_merged_v14_flat.jsonl backups/backup_$(date +\%Y\%m\%d).jsonl
```

### Updates

**To update application**:
```bash
# Pull latest code
git pull origin master

# Update dependencies
pip install --upgrade -r requirements.txt

# Restart service
sudo systemctl restart ligand-dashboard
```

### Monitoring

**Check service status**:
```bash
# Systemd
sudo systemctl status ligand-dashboard

# Docker
docker logs -f ligand-db

# Check resource usage
htop
```

**Logs location**:
- Streamlit logs: `~/.streamlit/logs/`
- System logs: `/var/log/syslog` or journalctl

---

## Minimal Deployment Checklist

- [ ] Python 3.8+ installed
- [ ] Create virtual environment
- [ ] Copy essential files:
  - [ ] `app.py`
  - [ ] `requirements.txt`
  - [ ] `data/sx_ligands_merged_v14_flat.jsonl`
  - [ ] `rel_quest/` directory
  - [ ] `vocab/` directory
- [ ] Install dependencies: `pip install -r requirements.txt`
- [ ] Test database loading: `python3 -m pytest -q`
- [ ] Launch dashboard: `streamlit run app.py`
- [ ] Verify in browser: `http://localhost:8501`
- [ ] Configure firewall (if remote access needed)
- [ ] Setup systemd service (for production)
- [ ] Configure Nginx reverse proxy (for HTTPS)
- [ ] Setup automated backups

---

## Quick Reference Commands

```bash
# Development
streamlit run app.py

# Production (systemd)
sudo systemctl start ligand-dashboard
sudo systemctl status ligand-dashboard
sudo journalctl -u ligand-dashboard -f

# Production (Docker)
docker-compose up -d
docker logs -f ligand-db

# Testing
python3 -m pytest -q
streamlit run app.py --server.headless true

# Debugging
streamlit run app.py --logger.level=debug
```

---

## Support & Documentation

- **Project overview**: `README.md`
- **Troubleshooting a self-hosted server**: `docs/DEPLOYMENT_CHECKLIST.md`

**Questions?** Check the documentation files or review git commit history for implementation details.

---

## Summary

**Minimum requirements to deploy**:
1. Python 3.8+ with pip
2. ~500MB disk space (app + dependencies)
3. Files: `app.py`, `requirements.txt`, `data/sx_ligands_merged_v14_flat.jsonl`, `rel_quest/`, `vocab/`
4. Run: `pip install -r requirements.txt && streamlit run app.py`

**That's it!** The dashboard will be accessible at `http://localhost:8501`
