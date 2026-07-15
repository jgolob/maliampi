FROM python:3.11-slim-bookworm AS builder

WORKDIR /app
ENV UV_LINK_MODE=copy
RUN pip install --no-cache-dir uv==0.11.7
COPY pyproject.toml uv.lock README.md LICENSE.txt ./
COPY src ./src
RUN uv sync --frozen --no-dev --no-editable

FROM python:3.11-slim-bookworm

WORKDIR /app
ENV PATH="/app/.venv/bin:$PATH" \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1
COPY --from=builder /app/.venv /app/.venv
