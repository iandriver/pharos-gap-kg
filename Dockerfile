# The question server over MCP's streamable HTTP, with the pharos319 bundle inside: read-only, no
# secrets, no outside calls. The bundle is unpacked from the release archive, each file checked against
# the digests committed in docs/data/bundle_archive_pharos319.json, so the image refuses to build on
# any other bundle. It serves plain HTTP on 8080 and expects a host that terminates TLS in front.
#
# From the repository root, with the archive at data/release/ (python -m pharos_gap_kg.release get pharos319):
#   docker build -t pharos-gap-kg .
#   docker run --rm -p 8080:8080 pharos-gap-kg                              # /mcp and /health, any Host
#   docker run --rm -p 8080:8080 pharos-gap-kg --public-host kg.example.org # only that Host
# For an x86-64 host, build with: docker buildx build --platform linux/amd64 -t pharos-gap-kg .

FROM python:3.12-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704 AS build
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY pyproject.toml constraints.txt README.md LICENSE ./
COPY src ./src
RUN python -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir ".[server]" -c constraints.txt

FROM python:3.12-slim@sha256:a0f57064cbdcd6285fa45f51237ab0b35ea8ebaed06a1f5973e89e19dfb86704
COPY --from=build /opt/venv /opt/venv
ENV PATH=/opt/venv/bin:$PATH PHAROS_GAP_KG_DATA=/data
COPY docs/data/bundle_archive_pharos319.json /data/docs/
COPY data/release/pharos-gap-kg-bundle-pharos319.tar.gz /tmp/
RUN python -m pharos_gap_kg.release unpack pharos319 /tmp/pharos-gap-kg-bundle-pharos319.tar.gz \
    && rm /tmp/pharos-gap-kg-bundle-pharos319.tar.gz
USER nobody
EXPOSE 8080
ENTRYPOINT ["python", "-m", "pharos_gap_kg.server", "--http", "--host", "0.0.0.0", "--port", "8080", "--behind-proxy", "--bundle", "/data/bundle/pharos319"]
CMD ["--public-host", "*"]
