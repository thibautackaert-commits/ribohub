#!/bin/bash
# Register RiboHub in Galaxy's tool_conf.xml (idempotent).
# Runs on every container start before Galaxy boots.

TOOL_CONF="/etc/galaxy/tool_conf.xml"
TOOL_LINE='  <tool file="/etc/galaxy/ribohub.xml" />'

# Wait for tool_conf.xml to exist (bgruening image may generate it)
for i in $(seq 1 30); do
    [ -f "$TOOL_CONF" ] && break
    sleep 1
done

if [ ! -f "$TOOL_CONF" ]; then
    echo "[RiboHub] WARNING: $TOOL_CONF not found. Tool not registered."
    exit 0
fi

# Only add if not already present
if ! grep -q "ribohub.xml" "$TOOL_CONF"; then
    sed -i "s|</toolbox>|${TOOL_LINE}\n</toolbox>|" "$TOOL_CONF"
    echo "[RiboHub] Registered in $TOOL_CONF"
else
    echo "[RiboHub] Already registered in $TOOL_CONF"
fi

# Fix permissions on mounted data so Galaxy's job runner can read them
chmod -R 755 /data/bigwig 2>/dev/null || true
chmod 644 /data/metadata.csv 2>/dev/null || true
chmod 777 /export/galaxy/hub_output 2>/dev/null || truechmod 777 /export/galaxy/hub_output 2>/dev/null || true
