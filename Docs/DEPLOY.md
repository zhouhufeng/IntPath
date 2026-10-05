# Deploying intpath.genohub.org

IntPathV2 runs on the `igvfagent-prod` Jetstream2 VM beside igvfagent.genohub.org and
igvfkg.genohub.org, the same way igvfkg does:

```
internet ──▶ Cloudflare (TLS, DNS) ──▶ igvfagent cloudflared tunnel (outbound only)
                                             │  intpath.genohub.org → http://gateway:80
                                             ▼
                          igvfagent-gateway (nginx) ── sites/intpath.conf (public, rate-limited)
                                             │
                                             ▼
                          intpath-web (compose project "intpath", network igvfagent_default)
                                             │  read-only
                                             ▼
                          /mnt/igvf-data/intpath/release/<organism>/intpath.sqlite
```

No port is published. IntPath is a public database, so its server block doesn't use the IGVF
sign-in gate; instead it rate-limits per client (10 req/s for the API, 30 enrichment runs per
minute). The container is read-only, drops all capabilities and is limited to 8 GB / 4 CPUs.

## Release and deploy

```bash
# 1. build the open-tier releases (no KEGG/BioCyc/restricted MSigDB: redistributable)
for org in sapiens musculus cerevisiae tuberculosis; do
  intpath build $org --public --out Data/intpathv2/release-open
done

# 2. ship code + releases to the VM and run Deploy/deploy.sh there
INTPATH_VM=ubuntu@<vm-floating-ip> INTPATH_SSH_KEY=<deploy key> Deploy/ship.sh            # all organisms
INTPATH_VM=... INTPATH_SSH_KEY=... Deploy/ship.sh sapiens                                  # one organism

# 3. once only: route the hostname through the igvfagent tunnel (idempotent; --check reports only)
CF_API_TOKEN=... CF_ACCOUNT_ID=... CF_TUNNEL_ID=<igvfagent tunnel id> python3 Deploy/cloudflare.py
```

The VM address, deploy key, tunnel id and Cloudflare token are in the IGVFagent operator notes
(`IGVFagent/Docs/Secretes/`), never in this repository.

`Deploy/deploy.sh` (run on the VM):
1. checks that release databases exist;
2. builds `intpath:latest` and starts `intpath-web`;
3. waits until the container is healthy;
4. installs `nginx/intpath.conf` into the gateway's `sites/`;
5. runs `nginx -t`, reloading only if it passes (the block is removed again on failure);
6. checks intpath, igvfagent and igvfkg health through the gateway.

## Operating

```bash
ssh ubuntu@<vm> 'docker logs -f intpath-web'                      # API log
ssh ubuntu@<vm> 'bash /srv/intpath/Deploy/deploy.sh --check'       # status only
curl https://intpath.genohub.org/healthz
```

Each uvicorn worker (`INTPATH_WORKERS`, default 2) holds every organism's compact library in
memory: about 1.1 GB per worker for human, mouse, yeast and *M. tuberculosis*. New releases are
picked up on container restart (`deploy.sh` restarts it).
