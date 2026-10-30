"""Sections of event_report.md written by the analysis commands (pb, sources, spectra)."""
from pathlib import Path


def update(home, title, lines):
    """Write (or replace) the section `## title` of event_report.md when the report exists, before the
    section Figures; returns the path or None."""
    f = Path(home) / "event_report.md"
    if not f.exists():
        return None
    text = f.read_text(encoding="utf-8")
    head = f"## {title}"
    sec = "\n".join([head, ""] + list(lines)) + "\n"
    if head + "\n" in text:
        before, rest = text.split(head + "\n", 1)
        nxt = rest.find("\n## ")
        text = before + sec + (rest[nxt + 1:] if nxt >= 0 else "")
    else:
        marker = "## Figures"
        text = text.replace(marker, sec + "\n" + marker, 1) if marker in text else text.rstrip("\n") + "\n\n" + sec
    f.write_text(text, encoding="utf-8")
    return f
