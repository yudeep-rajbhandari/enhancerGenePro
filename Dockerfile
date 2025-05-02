FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    bedtools \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Install Python packages directly
RUN pip install --no-cache-dir \
    pandas \
    pybedtools \
    seaborn \
    matplotlib \
    swifter \
    matplotlib-venn \
    Flask \
    waitress \
    flask-crontab \
    Flask-Mail

# Copy your app code
ENV FLASK_APP=waitress_server.py

COPY . .

CMD ["sh", "-c", "flask crontab add && python3 waitress_server.py"]

