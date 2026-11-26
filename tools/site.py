#!/usr/bin/env python3
"""Helper for writing posts. Needs only Python 3.8+, nothing to install.

    python3 tools/site.py new "Title of my post"   start a post from a template
    python3 tools/site.py serve                    live preview at http://localhost:8000
    python3 tools/site.py build                    refresh posts/index.json
    python3 tools/site.py build --out _site        also copy the publishable site there
                                                   (GitHub Actions runs this on every push)
"""

import argparse
import datetime as dt
import http.server
import json
import re
import shutil
import sys
import unicodedata
import webbrowser
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
POSTS = ROOT / "posts"
INDEX = POSTS / "index.json"

# Files and folders that are not part of the published site.
EXCLUDE = {
    ".git", ".github", ".gitignore", ".DS_Store", "tools", "_site",
    "README.md", "node_modules", "__pycache__",
}

FRONT_MATTER = re.compile(r"\A---[ \t]*\n(.*?)\n---[ \t]*(?:\n|\Z)", re.S)
ISO_DATE = re.compile(r"^\d{4}-\d{2}(-\d{2})?$")
NICE_SLUG = re.compile(r"^[A-Za-z0-9_-]+$")


# ---------------------------------------------------------------------------
# Front matter (same rules as parseFrontMatter in assets/blog.js)
# ---------------------------------------------------------------------------

def parse_scalar(raw):
    v = raw.strip()
    if len(v) > 1 and v[0] == v[-1] and v[0] in "\"'":
        return v[1:-1]
    if v.startswith("[") and v.endswith("]"):
        return [x for x in (parse_scalar(p) for p in v[1:-1].split(",")) if x != ""]
    if v.lower() in ("true", "yes"):
        return True
    if v.lower() in ("false", "no"):
        return False
    return v


def parse_front_matter(text):
    text = text.lstrip("﻿").replace("\r\n", "\n").replace("\r", "\n")
    m = FRONT_MATTER.match(text)
    if not m:
        return {}, text
    meta, list_key = {}, None
    for line in m.group(1).split("\n"):
        if not line.strip() or line.strip().startswith("#"):
            continue
        item = re.match(r"^\s*-\s+(.*)$", line)
        if item and list_key:
            meta[list_key].append(parse_scalar(item.group(1)))
            continue
        kv = re.match(r"^([A-Za-z_][\w-]*)\s*:\s*(.*)$", line)
        if not kv:
            continue
        key = kv.group(1).lower()
        if kv.group(2).strip() == "":
            meta[key] = []
            list_key = key
        else:
            meta[key] = parse_scalar(kv.group(2))
            list_key = None
    return meta, text[m.end():]


def as_text(v):
    if v is None or v is True or v is False:
        return ""
    if isinstance(v, list):
        return ", ".join(str(x) for x in v)
    return str(v).strip()


def as_list(v):
    if isinstance(v, list):
        return [str(x) for x in v if str(x).strip()]
    if isinstance(v, str) and v.strip():
        return [t.strip() for t in v.split(",") if t.strip()]
    return []


# ---------------------------------------------------------------------------
# Posts
# ---------------------------------------------------------------------------

def post_files():
    """Every posts/*.md except ones starting with _ or . (handy for templates)."""
    if not POSTS.is_dir():
        return []
    return [p for p in sorted(POSTS.glob("*.md")) if not p.name.startswith(("_", "."))]


def read_post(path):
    meta, _ = parse_front_matter(path.read_text(encoding="utf-8"))
    return meta


def collect_posts(include_drafts=False):
    posts, warnings = [], []
    for path in post_files():
        meta = read_post(path)
        draft = meta.get("draft") is True
        if draft and not include_drafts:
            continue
        slug = path.stem
        date = as_text(meta.get("date"))
        if not as_text(meta.get("title")):
            warnings.append(f"{path.name}: no title in the front matter, using '{slug}'")
        if not date:
            warnings.append(f"{path.name}: no date (add e.g. 'date: {dt.date.today().isoformat()}')")
        elif not ISO_DATE.match(date):
            warnings.append(f"{path.name}: date should look like 2026-01-31 or 2026-01, got '{date}'")
        if not NICE_SLUG.match(slug):
            warnings.append(f"{path.name}: works, but a name like 'my-post.md' gives nicer links")
        entry = {
            "slug": slug,
            "title": as_text(meta.get("title")) or slug,
            "date": date,
            "summary": as_text(meta.get("summary")),
            "tags": as_list(meta.get("tags")),
        }
        if as_text(meta.get("status")):
            entry["status"] = as_text(meta.get("status"))
        if draft:
            entry["draft"] = True
        posts.append(entry)
    posts.sort(key=lambda p: (p["date"], p["slug"]), reverse=True)
    return posts, warnings


def index_json(posts):
    return json.dumps(posts, indent=2, ensure_ascii=False) + "\n"


def write_if_changed(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return False
    path.write_text(text, encoding="utf-8")
    return True


def update_index(verbose=True):
    """Rewrite posts/index.json (published posts only)."""
    posts, warnings = collect_posts(include_drafts=False)
    changed = write_if_changed(INDEX, index_json(posts))
    if verbose:
        for w in warnings:
            print(f"  warning: {w}")
    return posts, changed


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def slugify(title):
    s = unicodedata.normalize("NFKD", title).encode("ascii", "ignore").decode()
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


TEMPLATE = """---
title: {title}
date: {date}
summary: One or two sentences for the list on the home page.
tags: []
draft: true
---

Start writing here. Inline math looks like $e^{i\\pi} + 1 = 0$ and display math like

$$
\\int_{-\\infty}^{\\infty} e^{-x^2}\\,dx = \\sqrt{\\pi}.
$$

::: theorem
Theorems, lemmas, definitions and proofs go in blocks like this one.
:::

The full list of what you can write is in posts/formatting-guide.md.
"""


def cmd_new(args):
    title = " ".join(args.title).strip()
    slug = args.slug or slugify(title)
    if not slug:
        sys.exit("Please give the post a title, e.g.  python3 tools/site.py new \"Galois theory\"")
    path = POSTS / f"{slug}.md"
    if path.exists():
        sys.exit(f"{path.relative_to(ROOT)} already exists.")
    shown_title = title
    if title[:1] in "\"'[{":
        shown_title = '"' + title.replace('"', "'") + '"'
    text = (TEMPLATE
            .replace("{title}", shown_title)
            .replace("{date}", dt.date.today().isoformat()))
    POSTS.mkdir(exist_ok=True)
    path.write_text(text, encoding="utf-8")
    print(f"Created {path.relative_to(ROOT)}")
    print("It's a draft: it shows up in the local preview but not on the live site")
    print("until you delete the 'draft: true' line (or change it to 'draft: false').")
    print()
    print("Preview it with:  python3 tools/site.py serve")
    print(f"           then open  http://localhost:8000/post.html?p={slug}")


class PreviewHandler(http.server.SimpleHTTPRequestHandler):
    extensions_map = {
        **http.server.SimpleHTTPRequestHandler.extensions_map,
        ".md": "text/markdown; charset=utf-8",
        ".json": "application/json; charset=utf-8",
        ".js": "text/javascript; charset=utf-8",
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self):
        path = self.path.split("?", 1)[0].split("#", 1)[0]
        if path == "/posts/index.json":
            # Built fresh on every request, drafts included, so new posts show up at once.
            update_index(verbose=False)
            posts, _ = collect_posts(include_drafts=True)
            body = index_json(posts).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()

    def log_message(self, fmt, *args):
        # The post page polls for changes every second; only report problems.
        if len(args) > 1 and str(args[1])[:1] in "45":
            sys.stderr.write(f"  {args[1]}  {self.path}\n")


def cmd_serve(args):
    posts, _ = update_index()
    httpd = None
    for port in range(args.port, args.port + 20):
        try:
            httpd = http.server.ThreadingHTTPServer(("127.0.0.1", port), PreviewHandler)
            break
        except OSError:
            continue
    if httpd is None:
        sys.exit(f"Couldn't find a free port from {args.port} upward.")
    port = httpd.server_address[1]
    url = f"http://localhost:{port}/"
    drafts = [p for p in collect_posts(include_drafts=True)[0] if p.get("draft")]
    print(f"Previewing the site at {url}   (Ctrl+C to stop)")
    print("Post pages reload by themselves when you save the .md file.")
    if drafts:
        print("Drafts (visible here only):")
        for d in drafts:
            print(f"  {url}post.html?p={d['slug']}")
    if args.open:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        httpd.server_close()


def cmd_build(args):
    posts, changed = update_index()
    print(f"posts/index.json: {len(posts)} published post(s){' (updated)' if changed else ''}")
    if not args.out:
        return

    out = Path(args.out).resolve()
    if out == ROOT or out in ROOT.parents or out == Path.home():
        sys.exit(f"Refusing to use {out} as the output folder.")
    if out.exists():
        shutil.rmtree(out)

    skip = set(EXCLUDE)
    if out.parent == ROOT:
        skip.add(out.name)
    unpublished = {p.name for p in POSTS.glob("*.md")} - {p.name for p in post_files()}
    unpublished |= {p.name for p in post_files() if read_post(p).get("draft") is True}

    def ignore(folder, names):
        here = Path(folder).resolve()
        ignored = set()
        for name in names:
            if here == ROOT and (name in skip or (name.startswith(".") and name != ".nojekyll")):
                ignored.add(name)
            elif here == POSTS and name in unpublished:
                ignored.add(name)
            elif name in ("__pycache__", ".DS_Store"):
                ignored.add(name)
        return ignored

    shutil.copytree(ROOT, out, ignore=ignore)
    print(f"Site copied to {out}")


def main():
    parser = argparse.ArgumentParser(
        description="Write, preview and build the website's posts.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command")

    p_new = sub.add_parser("new", help="start a new post")
    p_new.add_argument("title", nargs="+", help="the post's title")
    p_new.add_argument("--slug", help="file name to use instead of one made from the title")
    p_new.set_defaults(func=cmd_new)

    p_serve = sub.add_parser("serve", help="preview the site locally")
    p_serve.add_argument("--port", type=int, default=8000)
    p_serve.add_argument("--open", action="store_true", help="open the browser too")
    p_serve.set_defaults(func=cmd_serve)

    p_build = sub.add_parser("build", help="refresh posts/index.json (and copy the site with --out)")
    p_build.add_argument("--out", help="folder to copy the publishable site into, e.g. _site")
    p_build.set_defaults(func=cmd_build)

    args = parser.parse_args()
    if not getattr(args, "func", None):
        parser.print_help()
        return
    args.func(args)


if __name__ == "__main__":
    main()
