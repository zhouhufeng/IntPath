# Deploying intpath.genohub.org

```bash
# 1. build releases (one per organism; the KEGG KGML download is rate-limited, so the first build is slowest)
pip install -e '.[web]'
for org in sapiens musculus cerevisiae tuberculosis; do intpath build $org; done

# 2. run the service
INTPATH_RELEASE_ROOT=Data/intpathv2/release uvicorn web.app:app --host 0.0.0.0 --port 8000 --workers 2
#   or: docker build -t intpath . && docker run -p 8000:8000 -v $PWD/Data/intpathv2/release:/data intpath
```

Put a reverse proxy in front for TLS, for example Caddy:

```
intpath.genohub.org {
    encode gzip
    reverse_proxy 127.0.0.1:8000
}
```

Point a DNS A/AAAA record for `intpath.genohub.org` at the server. Each worker keeps every
loaded organism release in memory: about 1–2 GB for human with GO and PPIs. Interactive API
docs are served at `/docs`.

Before making the server public, check licensing (see `src/intpath/sources.py`). KEGG- and
BioCyc-derived files must not be offered for download without a licence. The web analysis can
use them only under the terms you hold. For a fully open deployment, build with
`intpath build <org> --public`.
