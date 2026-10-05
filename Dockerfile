FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
COPY web ./web
RUN pip install --no-cache-dir '.[web]' && useradd --system --uid 10001 intpath
USER intpath
ENV INTPATH_RELEASE_ROOT=/data
EXPOSE 8000
# workers: $WEB_CONCURRENCY (uvicorn default); each worker holds every organism's library in memory
CMD ["uvicorn", "web.app:app", "--host", "0.0.0.0", "--port", "8000", "--proxy-headers", "--forwarded-allow-ips", "*"]
