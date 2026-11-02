FROM python:3.11-slim

WORKDIR /app

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies, fully pinned (requirements.lock is verified against the
# test suite and the derived-field reproducibility check; see its header).
COPY requirements.lock .
RUN pip install --no-cache-dir -r requirements.lock

# Copy project
COPY . .

EXPOSE 8501

CMD ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501"]
