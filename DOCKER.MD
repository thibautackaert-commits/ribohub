# Docker Setup

Runs RiboHub as a Galaxy tool with automatic hub serving via Apache.
Galaxy generates track hubs, Apache serves them — one `docker compose up` and you're live.

## Architecture

```
User → Galaxy (port 8080)  → runs RiboHub → writes hub to shared volume
GWIPS → Apache (port 80)   → serves hub files + bigWig/bigBed data
                               ↑
                     shared volume (hub-output)
```

## Quick Start

### 1. Configure

```bash
cp .env.example .env
```

Edit `.env`:

```env
SERVER_IP=143.239.109.114          # public IP (GWIPS needs to reach this)
DATA_DIR=/path/to/data             # directory with sorted bigWig/bigBed files
METADATA_CSV=/path/to/metadata.csv # RiboSeqOrg metadata CSV
```

### 2. Start

```bash
docker compose up -d --build
```

Galaxy takes 2-3 minutes on first boot. Check progress:

```bash
docker logs ribohub-galaxy --tail 10 -f
```

### 3. Use

1. Open Galaxy at `http://localhost:8080`
2. Log in with `admin@galaxy.org` / `password`
3. Search "RiboHub" in the tool panel
4. Upload a sample list, type SRR IDs or use a metadata filter
5. Click **Execute**
6. Copy the hub URL from the output
7. Paste into GWIPS-viz → My Data → Track Hubs → My Hubs

Each run creates a unique URL under `/genomebrowsertrack/<random>/` so hubs never overwrite each other.

## File Layout

```
docker-compose.yml            # Galaxy + Apache orchestration
.env.example                  # Configuration template
docker/
├── Dockerfile.galaxy          # Extends Galaxy, pre-installs RiboHub
├── apache.conf                # Apache config with CORS for GWIPS
└── register_tool.sh           # Auto-registers tool on Galaxy boot
```

## Environment Variables

| Variable | Required | Description |
|---|---|---|
| `SERVER_IP` | Yes | Public IP/hostname (GWIPS needs to reach this) |
| `DATA_DIR` | Yes | Host path to sorted bigWig/bigBed directory |
| `METADATA_CSV` | Yes | Host path to RiboSeqOrg metadata CSV |
| `APACHE_PORT` | No | Apache port (default: 80) |
| `GALAXY_PORT` | No | Galaxy port (default: 8080) |

## How It Works

- **Galaxy** runs inside a container with RiboHub pre-installed via `pip install`.
  The tool is auto-registered on every boot via `register_tool.sh`.
- **Apache** serves the hub output directory and the bigWig/bigBed data with CORS
  headers so GWIPS can fetch files cross-origin.
- A **shared Docker volume** bridges the two containers: Galaxy writes to it,
  Apache serves from it. New hubs appear instantly.
- Each run gets a **unique random folder** (`genomebrowsertrack/<15-char-id>/`),
  so multiple users and runs never overwrite each other.

## Hub URL Structure

```
http://<SERVER_IP>/genomebrowsertrack/<random_id>/RiboSeqHub.hub.txt
```

The `bigDataUrl` entries inside the hub point to:

```
http://<SERVER_IP>/data/<SRR_path>/<filename>.bw
```

Both served by the same Apache container.

## Cleanup

Remove old hubs:

```bash
docker exec ribohub-apache find /usr/local/apache2/htdocs/genomebrowsertrack -maxdepth 1 -type d -mtime +30 -exec rm -rf {} +
```

## Troubleshooting

**Galaxy shows 502**: Still booting. Wait 2-3 minutes.

**Tool not in Galaxy**: Check registration:
```bash
docker exec ribohub-galaxy grep ribohub /etc/galaxy/tool_conf.xml
```

**Permission denied on hub_output**: Fix manually:
```bash
docker exec -u root ribohub-galaxy chmod 777 /export/galaxy/hub_output
```
