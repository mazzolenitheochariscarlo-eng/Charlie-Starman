"""Inline the engine scripts into the app page: python3 webapp/build.py -> webapp/dist/biothing-sound-lab.html"""
from pathlib import Path

here = Path(__file__).parent
page = (here / "app.html").read_text(encoding="utf-8")
for marker, name in (("/*ENGINE*/", "engine.js"), ("/*ROOMS*/", "rooms.js")):
    assert page.count(marker) == 1, marker
    page = page.replace(marker, (here / name).read_text(encoding="utf-8"))
out = here / "dist" / "biothing-sound-lab.html"
out.parent.mkdir(exist_ok=True)
out.write_text(page, encoding="utf-8")
print(out, f"{len(page) / 1024:.0f} KB")
