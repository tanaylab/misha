#!/bin/bash
# Install repo-tracked git hooks into .git/hooks/ as relative symlinks.
# Re-run safely; existing hooks are backed up unless they already point
# at the tracked source.
#
# Usage:
#   tools/install-hooks.sh

set -e

repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"

hooks_src_rel="tools/githooks"
hooks_dst=".git/hooks"

if [ ! -d "$hooks_src_rel" ]; then
    echo "ERROR: $hooks_src_rel not found (run from a misha checkout)." >&2
    exit 1
fi

# Symlinks live at .git/hooks/<name>, so the relative target is ../../tools/...
link_target_prefix="../../$hooks_src_rel"

ts=$(date +%s)
installed=()

for src in "$hooks_src_rel"/*; do
    [ -f "$src" ] || continue
    name=$(basename "$src")
    dst="$hooks_dst/$name"
    target="$link_target_prefix/$name"

    if [ -L "$dst" ] && [ "$(readlink "$dst")" = "$target" ]; then
        echo "ok        $dst -> $target (already linked)"
        continue
    fi

    if [ -e "$dst" ] || [ -L "$dst" ]; then
        backup="$dst.bak.$ts"
        mv "$dst" "$backup"
        echo "backup    $dst -> $backup"
    fi

    ln -s "$target" "$dst"
    chmod +x "$src"
    installed+=("$name")
    echo "installed $dst -> $target"
done

if [ ${#installed[@]} -gt 0 ]; then
    echo ""
    echo "Installed: ${installed[*]}"
fi
