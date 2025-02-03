# Use the uv image with Python 3.10 based on debian slim buster 
FROM python:3.10-slim-buster
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Set environment variables
ENV UV_COMPILE_BYTECODE=1 \
    BGDATA_LOCAL="/bgdatacache" \
    BGDATA_OFFLINE=TRUE \
    BBGLAB_HOME="/home/user/.config/bbglab/" \
    LC_ALL=C.UTF-8

# Set the working directory to /boostdm
WORKDIR /boostdm

# Stage necessary files into the container
COPY uv.lock uv.lock
COPY pyproject.toml pyproject.toml
COPY boostdm boostdm
COPY README.md README.md

# Install dependencies and build the project
RUN mkdir -p /root/.cache/uv \
    && cd /boostdm \
    && uv sync --frozen --no-dev \
    && uv build \
    && pip install dist/*.tar.gz \
    # Clean up unnecessary files
    && rm -rf /root/.cache/uv /boostdm/build /boostdm/.venv

RUN apt-get update &&
    apt-get -y install locales-all zlib1g-dev tabix bedtools